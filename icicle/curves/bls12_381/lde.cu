#ifndef _BLS12_381_LDE
#define _BLS12_381_LDE
#include <cuda.h>
#include "../../appUtils/ntt/lde.cu"
#include "../../appUtils/ntt/ntt.cuh"
#include "../../appUtils/vector_manipulation/ve_mod_mult.cuh"
#include "curve_config.cuh"

extern "C" BLS12_381::scalar_t *build_domain_cuda_bls12_381(uint32_t domain_size, uint32_t logn, bool inverse, size_t device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        cudaStreamCreate(&stream);
        if (inverse)
        {
            return fill_twiddle_factors_array(domain_size, BLS12_381::scalar_t::omega_inv(logn), stream);
        }
        else
        {
            return fill_twiddle_factors_array(domain_size, BLS12_381::scalar_t::omega(logn), stream);
        }
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());
        return nullptr;
    }
}

extern "C" int ntt_cuda_bls12_381(BLS12_381::scalar_t *arr, uint32_t n, bool inverse, Decimation decimation, size_t device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        cudaStreamCreate(&stream);
        return ntt_end2end_template<BLS12_381::scalar_t, BLS12_381::scalar_t>(arr, n, inverse, stream); // TODO: pass device_id
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());

        return -1;
    }
}

extern "C" int ecntt_cuda_bls12_381(BLS12_381::projective_t *arr, uint32_t n, bool inverse, Decimation decimation, size_t device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        cudaStreamCreate(&stream);
        return ntt_end2end_template<BLS12_381::projective_t, BLS12_381::scalar_t>(arr, n, inverse, stream); // TODO: pass device_id
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());
        return -1;
    }
}

extern "C" int ntt_batch_cuda_bls12_381(BLS12_381::scalar_t *arr, uint32_t arr_size, uint32_t batch_size, bool inverse, size_t device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        cudaStreamCreate(&stream);
        return ntt_end2end_batch_template<BLS12_381::scalar_t, BLS12_381::scalar_t>(arr, arr_size, batch_size, inverse, stream); // TODO: pass device_id
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());
        return -1;
    }
}

// extern "C" int ecntt_batch_cuda_bls12_381(BLS12_381::projective_t *arr, uint32_t arr_size, uint32_t batch_size, bool inverse, size_t device_id = 0)
// {
//     try
//     {
//         return ntt_end2end_batch_template<BLS12_381::projective_t, BLS12_381::scalar_t>(arr, arr_size, batch_size, inverse); // TODO: pass device_id
//     }
//     catch (const std::runtime_error &ex)
//     {
//         printf("error %s", ex.what());
//         return -1;
//     }
// }

extern "C" int interpolate_scalars_cuda_bls12_381(BLS12_381::scalar_t *d_out, BLS12_381::scalar_t *d_evaluations, BLS12_381::scalar_t *d_domain, unsigned n, unsigned device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        BLS12_381::scalar_t *_null = nullptr;
        return interpolate(d_out, d_evaluations, d_domain, n, false, _null, stream);
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());
        return -1;
    }
}

extern "C" int interpolate_scalars_batch_cuda_bls12_381(BLS12_381::scalar_t *d_out, BLS12_381::scalar_t *d_evaluations, BLS12_381::scalar_t *d_domain, unsigned n,
                                                        unsigned batch_size, size_t device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        BLS12_381::scalar_t *_null = nullptr;
        cudaStreamCreate(&stream);
        return interpolate_batch(d_out, d_evaluations, d_domain, n, batch_size, false, _null, stream);
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());
        return -1;
    }
}

// extern "C" int interpolate_points_cuda_bls12_381(BLS12_381::projective_t *d_out, BLS12_381::projective_t *d_evaluations, BLS12_381::scalar_t *d_domain, unsigned n, size_t device_id = 0)
// {
//     try
//     {
//         return interpolate(d_out, d_evaluations, d_domain, n);
//     }
//     catch (const std::runtime_error &ex)
//     {
//         printf("error %s", ex.what());
//         return -1;
//     }
// }

// extern "C" int interpolate_points_batch_cuda_bls12_381(BLS12_381::projective_t *d_out, BLS12_381::projective_t *d_evaluations, BLS12_381::scalar_t *d_domain,
//                                                        unsigned n, unsigned batch_size, size_t device_id = 0)
// {
//     try
//     {
//         return interpolate_batch(d_out, d_evaluations, d_domain, n, batch_size);
//     }
//     catch (const std::runtime_error &ex)
//     {
//         printf("error %s", ex.what());
//         return -1;
//     }
// }

extern "C" int evaluate_scalars_cuda_bls12_381(BLS12_381::scalar_t *d_out, BLS12_381::scalar_t *d_coefficients, BLS12_381::scalar_t *d_domain,
                                               unsigned domain_size, unsigned n, unsigned device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        BLS12_381::scalar_t *_null = nullptr;
        cudaStreamCreate(&stream);
        return evaluate(d_out, d_coefficients, d_domain, domain_size, n, false, _null, stream);
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());
        return -1;
    }
}

extern "C" int evaluate_scalars_batch_cuda_bls12_381(BLS12_381::scalar_t *d_out, BLS12_381::scalar_t *d_coefficients, BLS12_381::scalar_t *d_domain, unsigned domain_size,
                                                     unsigned n, unsigned batch_size, size_t device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        BLS12_381::scalar_t *_null = nullptr;
        cudaStreamCreate(&stream);
        auto result_code = evaluate_batch(d_out, d_coefficients, d_domain, domain_size, n, batch_size, false, _null, 0);
        cudaStreamDestroy(stream);
        return result_code;
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());
        return -1;
    }
}

extern "C" int bailey_ntt_cuda_bls12_381(BLS12_381::scalar_t *d_inout, BLS12_381::scalar_t *d_twf, BLS12_381::scalar_t *d_full_twf, uint32_t n, uint32_t batch_size, size_t device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        cudaStreamCreate(&stream);
        auto result_code = bailey_ntt(d_inout, d_twf, d_full_twf, n, batch_size, stream); // TODO: pass device_id
        cudaStreamDestroy(stream);
        return result_code;
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());
        return -1;
    }
}

extern "C" int bench_fr_add_cuda(size_t device_id, size_t samples, size_t blocks, size_t threads)
{
    // auto &gpu = select_gpu(device_id);
    BLS12_381::scalar_t f1 = BLS12_381::scalar_t::omega(8); // TODO: any value, random
    BLS12_381::scalar_t f2 = BLS12_381::scalar_t::omega(7);

    BLS12_381::scalar_t h_answer;
    BLS12_381::scalar_t *d_answer;
    cudaMalloc(&d_answer, sizeof(BLS12_381::scalar_t));

    bench_add_kernel<<<blocks, threads>>>(f1, f2, d_answer, (size_t)(blocks * threads), samples);
    cudaDeviceSynchronize();

    cudaMemcpy(&h_answer, d_answer, sizeof(BLS12_381::scalar_t), cudaMemcpyDeviceToHost);
    cudaFree(d_answer);
    return 0;
}

extern "C" int bench_fr_sub_cuda(size_t device_id, size_t samples)
{
    // auto &gpu = select_gpu(device_id);
    // fr_t f1 = group_gen;
    // fr_t f2 = f1 * group_gen_inverse;

    // fr_t t;

    // for (int s = 0; s < samples; s++)
    // {
    //     t = f1 - f2;
    // }

    // fr_t f = t;

    return 0;
}

extern "C" int bench_fr_mul_cuda(size_t device_id, size_t samples, size_t blocks, size_t threads)
{
    // auto &gpu = select_gpu(device_id);
    BLS12_381::scalar_t f1 = BLS12_381::scalar_t::omega(8); // TODO: any value, random
    BLS12_381::scalar_t f2 = BLS12_381::scalar_t::omega(7);

    BLS12_381::scalar_t h_answer;
    BLS12_381::scalar_t *d_answer;
    cudaMalloc(&d_answer, sizeof(BLS12_381::scalar_t));

    bench_mul_kernel<<<blocks, threads>>>(f1, f2, d_answer, (size_t)(blocks * threads), samples);
    cudaDeviceSynchronize();

    cudaMemcpy(&h_answer, d_answer, sizeof(BLS12_381::scalar_t), cudaMemcpyDeviceToHost);
    cudaFree(d_answer);
    return 0;
}



extern "C" int evaluate_points_cuda_bls12_381(BLS12_381::projective_t *d_out, BLS12_381::projective_t *d_coefficients, BLS12_381::scalar_t *d_domain,
                                              unsigned domain_size, unsigned n, size_t device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        BLS12_381::scalar_t *_null = nullptr;
        cudaStreamCreate(&stream);
        return evaluate(d_out, d_coefficients, d_domain, domain_size, n, false, _null, stream);
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());
        return -1;
    }
}

extern "C" int evaluate_points_batch_cuda_bls12_381(BLS12_381::projective_t *d_out, BLS12_381::projective_t *d_coefficients, BLS12_381::scalar_t *d_domain, unsigned domain_size,
                                                    unsigned n, unsigned batch_size, size_t device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        BLS12_381::scalar_t *_null = nullptr;
        cudaStreamCreate(&stream);
        auto result_code = evaluate_batch(d_out, d_coefficients, d_domain, domain_size, n, batch_size, false, _null, stream);
        cudaStreamDestroy(stream);
        return result_code;
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());
        return -1;
    }
}

extern "C" int evaluate_scalars_on_coset_cuda_bls12_381(BLS12_381::scalar_t *d_out, BLS12_381::scalar_t *d_coefficients, BLS12_381::scalar_t *d_domain, unsigned domain_size,
                                                        unsigned n, BLS12_381::scalar_t *coset_powers, unsigned device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        cudaStreamCreate(&stream);
        return evaluate(d_out, d_coefficients, d_domain, domain_size, n, true, coset_powers, stream);
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());
        return -1;
    }
}

extern "C" int evaluate_scalars_on_coset_batch_cuda_bls12_381(BLS12_381::scalar_t *d_out, BLS12_381::scalar_t *d_coefficients, BLS12_381::scalar_t *d_domain, unsigned domain_size,
                                                              unsigned n, unsigned batch_size, BLS12_381::scalar_t *coset_powers, size_t device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        cudaStreamCreate(&stream);
        return evaluate_batch(d_out, d_coefficients, d_domain, domain_size, n, batch_size, true, coset_powers, stream);
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());
        return -1;
    }
}

extern "C" int evaluate_points_on_coset_cuda_bls12_381(BLS12_381::projective_t *d_out, BLS12_381::projective_t *d_coefficients, BLS12_381::scalar_t *d_domain, unsigned domain_size,
                                                       unsigned n, BLS12_381::scalar_t *coset_powers, size_t device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        cudaStreamCreate(&stream); // TODO: don't create if default was passed, destroy what was created, same applies to all calls
        return evaluate(d_out, d_coefficients, d_domain, domain_size, n, true, coset_powers, stream);
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());
        return -1;
    }
}

extern "C" int evaluate_points_on_coset_batch_cuda_bls12_381(BLS12_381::projective_t *d_out, BLS12_381::projective_t *d_coefficients, BLS12_381::scalar_t *d_domain, unsigned domain_size,
                                                             unsigned n, unsigned batch_size, BLS12_381::scalar_t *coset_powers, size_t device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        cudaStreamCreate(&stream);
        return evaluate_batch(d_out, d_coefficients, d_domain, domain_size, n, batch_size, true, coset_powers, stream);
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());
        return -1;
    }
}

extern "C" int ntt_inplace_batch_cuda_bls12_381(BLS12_381::scalar_t *d_inout, BLS12_381::scalar_t *d_twiddles,
                                                unsigned n, unsigned batch_size, bool inverse, size_t device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        cudaStreamCreate(&stream);
        BLS12_381::scalar_t *_null = nullptr;
        ntt_inplace_batch_template(d_inout, d_twiddles, n, batch_size, inverse, false, _null, stream, true);
        return CUDA_SUCCESS; // TODO: we should implement this https://leimao.github.io/blog/Proper-CUDA-Error-Checking/
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());
        return -1;
    }
}

extern "C" int fast_ntt_batch_cuda_bls12_381(BLS12_381::scalar_t *d_inout, BLS12_381::scalar_t *d_twf, uint32_t n, uint32_t batch_size, size_t device_id = 0, cudaStream_t stream = 0)
{
    return ntt_inplace_batch_cuda_bls12_381(d_inout, d_twf, n, batch_size, false, device_id, stream);
}

extern "C" int reverse_order_scalars_cuda_bls12_381(BLS12_381::scalar_t *arr, int n, size_t device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        uint32_t logn = uint32_t(log(n) / log(2));
        cudaStreamCreate(&stream);
        reverse_order(arr, n, logn, stream);
        return 0;
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());
        return -1;
    }
}

extern "C" int reverse_order_scalars_batch_cuda_bls12_381(BLS12_381::scalar_t *arr, int n, int batch_size, size_t device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        uint32_t logn = uint32_t(log(n) / log(2));
        cudaStreamCreate(&stream);
        reverse_order_batch(arr, n, logn, batch_size, stream);
        return 0;
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());
        return -1;
    }
}

extern "C" int reverse_order_points_cuda_bls12_381(BLS12_381::projective_t *arr, int n, size_t device_id = 0, cudaStream_t stream = 0)
{
    try
    {
        uint32_t logn = uint32_t(log(n) / log(2));
        cudaStreamCreate(&stream);
        reverse_order(arr, n, logn, stream);
        return 0;
    }
    catch (const std::runtime_error &ex)
    {
        printf("error %s", ex.what());
        return -1;
    }
}

// extern "C" int reverse_order_points_cuda_bls12_381(BLS12_381::projective_t *arr, int n, size_t device_id = 0)
// {
//     try
//     {
//         uint32_t logn = uint32_t(log(n) / log(2));
//         reverse_order(arr, n, logn);
//         return 0;
//     }
//     catch (const std::runtime_error &ex)
//     {
//         printf("error %s", ex.what());
//         return -1;
//     }
// }

// extern "C" int reverse_order_points_batch_cuda_bls12_381(BLS12_381::projective_t *arr, int n, int batch_size, size_t device_id = 0)
// {
//     try
//     {
//         uint32_t logn = uint32_t(log(n) / log(2));
//         reverse_order_batch(arr, n, logn, batch_size);
//         return 0;
//     }
//     catch (const std::runtime_error &ex)
//     {
//         printf("error %s", ex.what());
//         return -1;
//     }
// }
#endif