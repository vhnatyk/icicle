use icicle_cuda_runtime::memory::HostOrDeviceSlice;

use crate::{
    curve::Curve,
    ntt::{FieldImpl, IcicleResult, NTTConfig, NTTDir},
};

pub use crate::curve::Projective;

#[cfg(feature = "arkworks")]
#[doc(hidden)]
pub mod tests;

#[doc(hidden)]
pub trait ECNTT<C: Curve>: ECNTTUnchecked<Projective<C>, C::ScalarField> {}

pub trait ECNTTUnchecked<T, F: FieldImpl> {
    fn ntt_unchecked(
        input: &(impl HostOrDeviceSlice<T> + ?Sized),
        dir: NTTDir,
        cfg: &NTTConfig<F>,
        output: &mut (impl HostOrDeviceSlice<T> + ?Sized),
    ) -> IcicleResult<()>;

    fn ntt_inplace_unchecked(
        inout: &mut (impl HostOrDeviceSlice<T> + ?Sized),
        dir: NTTDir,
        cfg: &NTTConfig<F>,
    ) -> IcicleResult<()>;
}

#[macro_export]
macro_rules! impl_ecntt {
    (
        $field_prefix:literal,
        $field_prefix_ident:ident,
        $field:ident,
        $field_config:ident,
        $curve:ident
    ) => {
        mod $field_prefix_ident {

            use crate::curve;
            use crate::curve::BaseCfg;
            use crate::ecntt::IcicleResult;
            use crate::ecntt::Projective;
            use crate::ecntt::{
                $curve, $field, $field_config, CudaError, DeviceContext, NTTConfig, NTTDir, DEFAULT_DEVICE_ID,
            };
            use icicle_core::ecntt::ECNTT;
            use icicle_core::ecntt::ECNTTUnchecked;
            use icicle_core::impl_ntt_without_domain;
            use icicle_core::ntt::NTT;
            use icicle_core::traits::IcicleResultWrap;
            use icicle_cuda_runtime::memory::HostOrDeviceSlice;

            pub type ProjectiveC = Projective<$curve>;
            impl_ntt_without_domain!($field_prefix, $field, $field_config, ECNTTUnchecked, "ECNTT", ProjectiveC);

            fn ntt_unchecked(
                input: &(impl HostOrDeviceSlice<Projective<$curve>> + ?Sized),
                dir: NTTDir,
                cfg: &NTTConfig<$field>,
                output: &mut (impl HostOrDeviceSlice<Projective<$curve>> + ?Sized),
            ) -> IcicleResult<()> {
                <curve::ScalarCfg as ECNTTUnchecked<Projective<$curve>, $field>>::ntt_unchecked(input, dir, cfg, output)
            }

            fn ntt_inplace_unchecked(
                inout: &mut (impl HostOrDeviceSlice<Projective<$curve>> + ?Sized),
                dir: NTTDir,
                cfg: &NTTConfig<$field>,
            ) -> IcicleResult<()> {
                <curve::ScalarCfg as ECNTTUnchecked<Projective<$curve>, $field>>::ntt_inplace_unchecked(inout, dir, cfg)
            }

            impl ECNTT<$curve> for $field_config {}
        }
    };
}

/// Computes the ECNTT, or a batch of several NTTs.
///
/// # Arguments
///
/// * `input` - inputs of the NTT.
///
/// * `dir` - whether to compute forward of inverse NTT.
///
/// * `cfg` - config used to specify extra arguments of the NTT.
///
/// * `output` - buffer to write the NTT outputs into. Must be of the same size as `input`.
pub fn ecntt<C: Curve>(
    input: &(impl HostOrDeviceSlice<Projective<C>> + ?Sized),
    dir: NTTDir,
    cfg: &NTTConfig<C::ScalarField>,
    output: &mut (impl HostOrDeviceSlice<Projective<C>> + ?Sized),
) -> IcicleResult<()>
where
    C::ScalarField: FieldImpl,
    <C::ScalarField as FieldImpl>::Config: ECNTT<C>,
{
    <<C::ScalarField as FieldImpl>::Config as ECNTTUnchecked<Projective<C>, C::ScalarField>>::ntt_unchecked(
        input, dir, &cfg, output,
    )
}

/// Computes the ECNTT, or a batch of several NTTs inplace.
///
/// # Arguments
///
/// * `inout` - buffer with inputs to also write the NTT outputs into.
///
/// * `dir` - whether to compute forward of inverse NTT.
///
/// * `cfg` - config used to specify extra arguments of the NTT.
pub fn ecntt_inplace<C: Curve>(
    inout: &mut (impl HostOrDeviceSlice<Projective<C>> + ?Sized),
    dir: NTTDir,
    cfg: &NTTConfig<C::ScalarField>,
) -> IcicleResult<()>
where
    C::ScalarField: FieldImpl,
    <C::ScalarField as FieldImpl>::Config: ECNTT<C>,
{
    <<C::ScalarField as FieldImpl>::Config as ECNTTUnchecked<Projective<C>, C::ScalarField>>::ntt_inplace_unchecked(
        inout, dir, &cfg,
    )
}

#[macro_export]
macro_rules! impl_ecntt_tests {
    (
      $field:ident,
      $curve:ident
    ) => {
        use icicle_core::ntt::tests::init_domain;
        use icicle_cuda_runtime::device_context::DEFAULT_DEVICE_ID;
        const MAX_SIZE: u64 = 1 << 18;
        static INIT: OnceLock<()> = OnceLock::new();
        const FAST_TWIDDLES_MODE: bool = false;

        #[test]
        fn test_ecntt() {
            INIT.get_or_init(move || init_domain::<$field>(MAX_SIZE, DEFAULT_DEVICE_ID, FAST_TWIDDLES_MODE));
            check_ecntt::<$curve>()
        }

        #[test]
        fn test_ecntt_batch() {
            INIT.get_or_init(move || init_domain::<$field>(MAX_SIZE, DEFAULT_DEVICE_ID, FAST_TWIDDLES_MODE));
            check_ecntt_batch::<$curve>()
        }

        // #[test]
        // fn test_ntt_device_async() {
        //     // init_domain is in this test is performed per-device
        //     check_ntt_device_async::<$field>()
        // }
    };
}

#[macro_export]
macro_rules! impl_ecntt_bench {
    (
      $field_prefix:literal,
      $field:ident,
      $base_field:ident,
      $curve:ident
    ) => {
        use std::sync::OnceLock;

        #[cfg(feature = "arkworks")]
        use ark_ff::FftField;

        use criterion::{black_box, criterion_group, criterion_main, Criterion};
        use icicle_core::{
            traits::ArkConvertible,
            curve::Curve,
            ecntt::{ecntt, Projective},
            ntt::{FieldImpl, HostOrDeviceSlice, NTTConfig, NTTDir, NttAlgorithm, Ordering},
        };

        use icicle_core::ecntt::ECNTT;
        use icicle_core::ntt::NTT;

        fn ecntt_for_bench<F: FieldImpl + ArkConvertible, C: Curve>(
            points: &HostOrDeviceSlice<Projective<C>>,
            mut batch_ntt_result: &mut HostOrDeviceSlice<Projective<C>>,
            test_sizes: usize,
            batch_size: usize,
            is_inverse: NTTDir,
            ordering: Ordering,
            config: &mut NTTConfig<C::ScalarField>,
            _seed: u32,
        ) where
            <C::BaseField as FieldImpl>::Config: NTT<Projective<C>, F>,
        {
            ntt(&points, is_inverse, config, &mut batch_ntt_result).unwrap();
        }

        static INIT: OnceLock<()> = OnceLock::new();

        fn benchmark_ecntt<F: FieldImpl + ArkConvertible, C: Curve>(c: &mut Criterion)
        where
            F::ArkEquivalent: FftField,
            <F as FieldImpl>::Config: NTT<Projective<C>, F>,
            <F as FieldImpl>::Config: NTTDomain<F>
        {
            use icicle_core::ntt::tests::init_domain;
            use icicle_cuda_runtime::device_context::DEFAULT_DEVICE_ID;
            use criterion::SamplingMode;

            let group_id = format!("{} EC NTT", $field_prefix);
            let mut group = c.benchmark_group(&group_id);
            group.sampling_mode(SamplingMode::Flat);
            group.sample_size(10);

            const MAX_SIZE: u64 = 1 << 18;
            const FAST_TWIDDLES_MODE: bool = false;
            INIT.get_or_init(move || init_domain::<F>(MAX_SIZE, DEFAULT_DEVICE_ID, FAST_TWIDDLES_MODE));

            let test_sizes = [1 << 4, 1 << 8];
            let batch_sizes = [1, 1 << 4, 128];
            for test_size in test_sizes {
                for batch_size in batch_sizes {
                    let points = HostOrDeviceSlice::on_host(C::generate_random_projective_points(test_size * batch_size));
                    let mut batch_ntt_result = HostOrDeviceSlice::on_host(vec![Projective::zero(); batch_size * test_size]);
                    let mut config = NTTConfig::default();
                    for is_inverse in [NTTDir::kInverse, NTTDir::kForward] {
                        for ordering in [
                            Ordering::kNN,
                            // Ordering::kNR, // times are ~ same as kNN
                            // Ordering::kRN,
                            // Ordering::kRR,
                            // Ordering::kNM, // no mixed radix ecntt
                            // Ordering::kMN,
                        ] {
                            config.ordering = ordering;
                            for alg in [NttAlgorithm::Radix2] {
                                config.batch_size = batch_size as i32;
                                config.ntt_algorithm = alg;
                                let bench_descr = format!("{:?} {:?} {:?} {} x {}", alg, ordering, is_inverse, test_size, batch_size);
                                group.bench_function(&bench_descr, |b| {
                                    b.iter(|| {
                                        ecntt_for_bench::<F, C>(
                                            &points,
                                            &mut batch_ntt_result,
                                            test_size,
                                            batch_size,
                                            is_inverse,
                                            ordering,
                                            &mut config,
                                            black_box(1),
                                        )
                                    })
                                });
                            }
                        }
                    }
                }
            }

            group.finish();
        }

        criterion_group!(benches, benchmark_ecntt<$field, $curve>);
        criterion_main!(benches);
    };
}