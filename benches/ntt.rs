extern crate criterion;

use std::os::unix::raw::dev_t;

use criterion::{
    criterion_group, criterion_main, measurement::WallTime, BenchmarkGroup, Criterion,
};

use icicle_utils::{
    curves::bls12_381::{Point_BLS12_381, ScalarField_BLS12_381},
    test_bls12_381::{
        bailey_ntt_bls12_381, bench_add_fr, bench_mul_fr, evaluate_points_batch_bls12_381,
        evaluate_scalars_batch_bls12_381, fast_ntt_batch_bls12_381,
        interpolate_points_batch_bls12_381, interpolate_scalars_batch_bls12_381,
        set_up_points_bls12_381, set_up_scalars_bls12_381,
    },
};
use rustacuda::prelude::DeviceBuffer;
type BenchGroup<'a> = BenchmarkGroup<'a, WallTime>;

const LOG_NTT_SIZES: [usize; 1] = [10]; //, 23, 9, 10, 11, 12, 18];
                                        // const BATCH_SIZES: [usize; 2] = [1<<17, 256];
const BATCH_SIZES: [usize; 1] = [1 << 10];

const MAX_POINTS_LOG2: usize = 18;
const MAX_SCALARS_LOG2: usize = 25;

fn bench_ntt(c: &mut Criterion) {
    for log_ntt_size in LOG_NTT_SIZES {
        for batch_size in BATCH_SIZES {
            let ntt_size = 1 << log_ntt_size;
            let group = &mut c.benchmark_group("NTT");

            fn fast_ntt(
                d_evaluations: &mut DeviceBuffer<ScalarField_BLS12_381>,
                d_domain: &mut DeviceBuffer<ScalarField_BLS12_381>,
                _d_domain_full: &mut DeviceBuffer<ScalarField_BLS12_381>,
                batch_size: usize,
            ) -> DeviceBuffer<ScalarField_BLS12_381> {
                //bailey_ntt_bls12_381(d_evaluations, d_domain, batch_size);
                //println!("domain: {} {}", d_domain.len(), batch_size);
                fast_ntt_batch_bls12_381(d_evaluations, d_domain, batch_size);

                unsafe { DeviceBuffer::uninitialized(d_domain.len()).unwrap() }
            }

            fn bailey_ntt(
                d_evaluations: &mut DeviceBuffer<ScalarField_BLS12_381>,
                d_domain: &mut DeviceBuffer<ScalarField_BLS12_381>,
                d_domain_full: &mut DeviceBuffer<ScalarField_BLS12_381>,
                _batch_size: usize,
            ) -> DeviceBuffer<ScalarField_BLS12_381> {
                bailey_ntt_bls12_381(
                    d_evaluations,
                    d_domain,
                    d_domain_full,
                    d_domain_full.len() / d_domain.len(),
                );
                unsafe { DeviceBuffer::uninitialized(d_domain.len()).unwrap() }
            }

            pub fn set_up_bailey_scalars_bls12_381(
                test_size: usize,
                log_domain_size: usize,
                inverse: bool,
            ) -> (
                Vec<ScalarField_BLS12_381>,
                DeviceBuffer<ScalarField_BLS12_381>,
                DeviceBuffer<ScalarField_BLS12_381>,
                Option<DeviceBuffer<ScalarField_BLS12_381>>,
            ) {
                let (a, b, c) = set_up_scalars_bls12_381(test_size, log_domain_size / 2, inverse);
                let (_, _, d) = set_up_scalars_bls12_381(0, log_domain_size, inverse);
                (a, b, c, Some(d))
            }

            pub fn set_up_ntt_scalars_bls12_381(
                test_size: usize,
                log_domain_size: usize,
                inverse: bool,
            ) -> (
                Vec<ScalarField_BLS12_381>,
                DeviceBuffer<ScalarField_BLS12_381>,
                DeviceBuffer<ScalarField_BLS12_381>,
                Option<DeviceBuffer<ScalarField_BLS12_381>>,
            ) {
                let (a, b, c) = set_up_scalars_bls12_381(test_size, log_domain_size, inverse);
                (a, b, c, None)
            }

            pub fn set_up_ecntt_point_bls12_381(
                test_size: usize,
                log_domain_size: usize,
                inverse: bool,
            ) -> (
                Vec<Point_BLS12_381>,
                DeviceBuffer<Point_BLS12_381>,
                DeviceBuffer<ScalarField_BLS12_381>,
                Option<&'static mut DeviceBuffer<ScalarField_BLS12_381>>,
            ) {
                let (a, b, c) = set_up_points_bls12_381(test_size, log_domain_size, inverse);
                (a, b, c, None)
            }

            bench_template(
                MAX_SCALARS_LOG2,
                ntt_size,
                batch_size,
                log_ntt_size,
                set_up_ntt_scalars_bls12_381,
                fast_ntt,
                group,
                "fast NTT",
                false,
                100,
            );

            bench_template(
                MAX_SCALARS_LOG2,
                ntt_size * batch_size,
                1,
                log_ntt_size * 2,
                set_up_bailey_scalars_bls12_381,
                bailey_ntt,
                group,
                "Bailey NTT",
                false,
                100,
            );

            // bench_template(
            //     MAX_SCALARS_LOG2,
            //     ntt_size,
            //     batch_size,
            //     log_ntt_size,
            //     set_up_ntt_scalars_bls12_381,
            //     interpolate_scalars_batch_bls12_381,
            //     group,
            //     "iNTT",
            //     true,
            //     100,
            // );
        }
    }
}

fn bench_template<E, S>(
    log_max_size: usize,
    ntt_size: usize,
    batch_size: usize,
    log_ntt_size: usize,
    set_data: fn(
        test_size: usize,
        log_domain_size: usize,
        inverse: bool,
    ) -> (
        Vec<E>,
        DeviceBuffer<E>,
        DeviceBuffer<S>,
        Option<DeviceBuffer<S>>,
    ),
    bench_fn: fn(
        d_evaluations: &mut DeviceBuffer<E>,
        d_domain: &mut DeviceBuffer<S>,
        d_domain_full: &mut DeviceBuffer<S>,
        batch_size: usize,
    ) -> DeviceBuffer<E>,
    g: &mut BenchmarkGroup<WallTime>,
    id: &str,
    inverse: bool,
    samples: usize,
) -> Option<(Vec<E>, Option<DeviceBuffer<E>>)> {
    let count = ntt_size * batch_size;

    let bench_id = format!("{} of size 2^{} in batch {}", id, log_ntt_size, batch_size);

    if count > 1 << log_max_size {
        println!("Bench size exceeded: {}", bench_id);
        return None;
    }

    let (input, mut d_evals, mut d_domain, d_domain_full) =
        set_data(ntt_size * batch_size, log_ntt_size, inverse);
    //let (_, _, mut d_domain_b) = set_data(ntt_size * batch_size, log_ntt_size/2, inverse);

    //let first = bench_fn(&mut d_evals, &mut d_domain, d_domain_full, batch_size);
    let mut d_domain_full = d_domain_full
        .unwrap_or_else(|| unsafe { DeviceBuffer::uninitialized(d_domain.len()).unwrap() });
    g.sample_size(samples).bench_function(&bench_id, |b| {
        b.iter(|| bench_fn(&mut d_evals, &mut d_domain, &mut d_domain_full, batch_size))
    });

    Some((input, None))
}

fn bench_fr_arith(c: &mut Criterion) {
    use std::str::FromStr;
    let bench_npow = std::env::var("ARITH_BENCH_NPOW").unwrap_or("8".to_string());
    let npoints_npow = usize::from_str(&bench_npow).unwrap();

    let mut group: BenchGroup = c.benchmark_group("FR CUDA");
    group.sample_size(10);

    for blocks in [128, 256, 1024, 2048] {
        for threads in [128] {
            for lg_domain_size in 2..=npoints_npow {
                let domain_size = 10_usize.pow(lg_domain_size as u32) as usize;

                let name = format!("FR ADD 10**{}*{}*{}", lg_domain_size, blocks, threads);
                group.bench_function(name, |b| {
                    b.iter(|| {
                        bench_add_fr(domain_size, blocks, threads);
                    })
                });

                let name = format!("FR MUL 10**{}*{}*{}", lg_domain_size, blocks, threads);
                group.bench_function(name, |b| {
                    b.iter(|| {
                        bench_mul_fr(domain_size, blocks, threads);
                    })
                });
            }
        }
    }
    group.finish();
}

criterion_group!(ntt_benches, bench_ntt, bench_fr_arith);
criterion_main!(ntt_benches);
