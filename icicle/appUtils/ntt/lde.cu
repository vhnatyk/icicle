#pragma once
#ifndef LDE
#define LDE
#include <cuda.h>
#include "ntt.cuh"
#include "lde.cuh"
#include "../vector_manipulation/ve_mod_mult.cuh"
#include "../../utils/tr.cuh"
#include <assert.h>

/**
 * Interpolate a batch of polynomials from their evaluations on the same subgroup.
 * Note: this function does not preform any bit-reverse permutations on its inputs or outputs.
 * @param d_out The variable to write coefficients of the resulting polynomials into (the coefficients are in bit-reversed order if the evaluations weren't bit-reversed and vice-versa).
 * @param d_evaluations Input array of evaluations of all polynomials of type E (elements).
 * @param d_domain Domain on which the polynomials are evaluated. Must be a subgroup.
 * @param n Length of `d_domain` array, also equal to the number of evaluations of each polynomial.
 * @param batch_size The size of the batch; the length of `d_evaluations` is `n` * `batch_size`.
 */
template <typename E, typename S>
int interpolate_batch(E *d_out, E *d_evaluations, S *d_domain, unsigned n, unsigned batch_size)
{
  uint32_t logn = uint32_t(log(n) / log(2));
  cudaMemcpy(d_out, d_evaluations, sizeof(E) * n * batch_size, cudaMemcpyDeviceToDevice);

  int NUM_THREADS = min(n / 2, MAX_THREADS_BATCH);
  int chunks = max(int((n / 2) / NUM_THREADS), 1);
  int total_tasks = batch_size * chunks;
  int NUM_BLOCKS = total_tasks;
  int max_sharedmem = 512 * sizeof(E);
  int shared_mem = 2 * NUM_THREADS * sizeof(E); // TODO: calculator, as shared mem size may be more efficient less then max to allow more concurrent blocks on SM
  uint32_t logn_shmem = uint32_t(log(2 * NUM_THREADS) / log(2));
  // ntt_template_kernel_shared<<<NUM_BLOCKS, NUM_THREADS, shared_mem, 0>>>(d_out, 1 << logn_shmem, d_domain, n, total_tasks, 0, logn_shmem, false);

  // for (uint32_t s = logn_shmem; s < logn; s++) // TODO: this loop also can be unrolled
  //for (uint32_t s = logn; s > 0; s--) // TODO: this loop also can be unrolled
  for (uint32_t s = 0; s < logn; s++) // TODO: this loop also can be unrolled
  {
    ntt_template_kernel_bc<<<NUM_BLOCKS, NUM_THREADS>>>(d_out, n, d_domain, n, total_tasks, s, false);
  }

  NUM_BLOCKS = (n * batch_size + NUM_THREADS - 1) / NUM_THREADS;
  template_normalize_kernel<<<NUM_BLOCKS, NUM_THREADS>>>(d_out, n * batch_size, S::inv_log_size(logn));
  return 0;
}

/**
 * Interpolate a polynomial from its evaluations on a subgroup.
 * Note: this function does not preform any bit-reverse permutations on its inputs or outputs.
 * @param d_out The variable to write coefficients of the resulting polynomial into (the coefficients are in bit-reversed order if the evaluations weren't bit-reversed and vice-versa).
 * @param d_evaluations Input array of evaluations that have type E (elements).
 * @param d_domain Domain on which the polynomial is evaluated. Must be a subgroup.
 * @param n Length of `d_evaluations` and the size `d_domain` arrays (they should have equal length).
 */
template <typename E, typename S>
int interpolate(E *d_out, E *d_evaluations, S *d_domain, unsigned n)
{
  return interpolate_batch<E, S>(d_out, d_evaluations, d_domain, n, 1);
}

template <typename E>
__global__ void fill_array(E *arr, E val, uint32_t n)
{
  int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (tid < n)
  {
    arr[tid] = val;
  }
}

template <typename E, typename S>
__global__ void bench_mul_kernel(E a, S b, E *r, size_t n, size_t samples)
{
  // S f1 = group_gen;
  // S f2 = f1 * group_gen_inverse;

  int tid = blockDim.x * blockIdx.x + threadIdx.x;
  if (tid < n)
  {
    // int scalar_id = tid % n_scalars;
    // element_vec[tid] = scalar_vec[scalar_id] * element_vec[tid];

    S t;

    for (int s2 = 0; s2 < samples; s2++)
    {
      t = t * b;
    }

    t = a * t;

    if (tid == 0)
    {
      *r = t;
    }
  }
}

template <typename E, typename S>
__global__ void bench_add_kernel(E a, S b, E *r, size_t n, size_t samples)
{
  // S f1 = group_gen;
  // S f2 = f1 * group_gen_inverse;

  int tid = blockDim.x * blockIdx.x + threadIdx.x;
  if (tid < n)
  {
    // int scalar_id = tid % n_scalars;
    // element_vec[tid] = scalar_vec[scalar_id] * element_vec[tid];

    S t;
    // for (int s1 = 0; s1 < samples; s1++)
    // {
    for (int s2 = 0; s2 < samples; s2++)
    {
      t = t + b;
    }
    // }

    t = a + t;

    if (tid == 0)
    {
      *r = t;
    }
  }
}

/**
 * Evaluate a batch of polynomials on the same coset.
 * @param d_out The evaluations of the polynomials on coset `u` * `d_domain`.
 * @param d_coefficients Input array of coefficients of all polynomials of type E (elements) to be evaluated in-place on a coset.
 * @param d_domain Domain on which the polynomials are evaluated (see `coset` flag). Must be a subgroup.
 * @param domain_size Length of `d_domain` array, on which the polynomial is computed.
 * @param n The number of coefficients, which might be different from `domain_size`.
 * @param batch_size The size of the batch; the length of `d_coefficients` is `n` * `batch_size`.
 * @param coset The flag that indicates whether to evaluate on a coset. If false, evaluate on a subgroup `d_domain`.
 * @param coset_powers If `coset` is true, a list of powers `[1, u, u^2, ..., u^{n-1}]` where `u` is the generator of the coset.
 */
template <typename E, typename S>
int evaluate_batch(E *d_out, E *d_coefficients, S *d_domain, unsigned domain_size, unsigned n, unsigned batch_size, bool coset, S *coset_powers)
{
  uint32_t logn = uint32_t(log(domain_size) / log(2));
  if (domain_size > n)
  {
    // allocate and initialize an array of stream handles to parallelize data copying across batches
    cudaStream_t *memcpy_streams = (cudaStream_t *)malloc(batch_size * sizeof(cudaStream_t));
    for (int i = 0; i < batch_size; i++)
    {
      cudaStreamCreate(&(memcpy_streams[i]));

      cudaMemcpyAsync(&d_out[i * domain_size], &d_coefficients[i * n], n * sizeof(E), cudaMemcpyDeviceToDevice, memcpy_streams[i]);
      int NUM_THREADS = MAX_THREADS_BATCH;
      int NUM_BLOCKS = (domain_size - n + NUM_THREADS - 1) / NUM_THREADS;
      fill_array<E><<<NUM_BLOCKS, NUM_THREADS, 0, memcpy_streams[i]>>>(&d_out[i * domain_size + n], E::zero(), domain_size - n);

      cudaStreamSynchronize(memcpy_streams[i]);
      cudaStreamDestroy(memcpy_streams[i]);
    }
  }
  else
    cudaMemcpy(d_out, d_coefficients, sizeof(E) * domain_size * batch_size, cudaMemcpyDeviceToDevice);

  if (coset)
    batch_vector_mult(coset_powers, d_out, domain_size, batch_size);

  int NUM_THREADS = min(domain_size / 2, MAX_THREADS_BATCH);
  int chunks = max(int((domain_size / 2) / NUM_THREADS), 1);
  int total_tasks = batch_size * chunks;
  int NUM_BLOCKS = total_tasks;
  int max_sharedmem = 512 * sizeof(E);
  int shared_mem = (2 * NUM_THREADS) * sizeof(E); // TODO: calculator, as shared mem size may be more efficient less then max to allow more concurrent blocks on SM
  uint32_t logn_shmem = uint32_t(log(2 * NUM_THREADS) / log(2));
  // for (uint32_t s = logn - 1; s >= logn_shmem; s--) // TODO: this loop also can be unrolled
  for (uint32_t s = logn; s > 0; s--) // TODO: this loop also can be unrolled
  // for (uint32_t s = 0; s < logn; s++) // TODO: this loop also can be unrolled
  {
    // ntt_template_kernel<<<NUM_BLOCKS, NUM_THREADS>>>(d_out, domain_size, d_domain, domain_size, total_tasks, s - 1, true);
    //  ntt_template_kernel<<<NUM_BLOCKS, NUM_THREADS>>>(d_out, domain_size, d_domain, domain_size, total_tasks, s, false);

    ntt_template_kernel_bc<<<NUM_BLOCKS, NUM_THREADS>>>(d_out, n, d_domain, n, total_tasks, s - 1, true);
  }

  // uint32_t log2_num_blocks = (log(NUM_BLOCKS) / log(2));
  // uint32_t n_div_log2_blocks = (((1 << logn_shmem) >> (log2_num_blocks + 1)) - 1);
  // uint32_t num_blocks2x = NUM_BLOCKS * 2; // TODO: ? uint32_t

  // ntt_template_kernel_shared_rev<<<NUM_BLOCKS, NUM_THREADS, shared_mem, 0>>>(d_out, 1 << logn_shmem, d_domain, n / 2, total_tasks, 0, logn_shmem - 1, n_div_log2_blocks, num_blocks2x, (1 << logn_shmem) - 1);
  // ntt_template_kernel_shared<<<NUM_BLOCKS, NUM_THREADS, shared_mem, 0>>>(d_out, 1 << logn_shmem, d_domain, n, total_tasks, 0, logn_shmem, false);

  return 0;
}

///
/**
 * Evaluate a batch of polynomials on the same coset.
 * @param d_inout Input array of type E (elements)
 * @param d_twf Twiddle factors of type S (scalars) array allocated on the device memory (must be a power of 2).
 * @param n The size of single input.
 * @param batch_size The size of the batch; the length of `d_inout` is `n` * `batch_size`.
 */
template <typename E, typename S>
int ntt_batch_template(E *d_inout, S *d_twf, unsigned n, unsigned batch_size)
{
  uint32_t logn = uint32_t(log(n) / log(2));

  int NUM_THREADS = min(n / 2, MAX_THREADS_BATCH);
  int chunks = max(int((n / 2) / NUM_THREADS), 1);
  int total_tasks = batch_size * chunks;
  int NUM_BLOCKS = total_tasks;
  int max_sharedmem = 512 * sizeof(E);
  int shared_mem = (2 * NUM_THREADS) * sizeof(E); // TODO: calculator, as shared mem size may be more efficient less then max to allow more concurrent blocks on SM
  uint32_t logn_shmem = uint32_t(log(2 * NUM_THREADS) / log(2));
  // for (uint32_t s = logn - 1; s >= logn_shmem; s--) // TODO: this loop also can be unrolled
  // for (uint32_t s = 0; s < logn; s++) // TODO: this loop also can be unrolled
  for (uint32_t s = logn; s > 0; s--) // TODO: this loop also can be unrolled
  {
    ntt_template_kernel<<<NUM_BLOCKS, NUM_THREADS>>>(d_inout, n, d_twf, n, total_tasks, s - 1, true);
  }

  // uint32_t log2_num_blocks = (log(NUM_BLOCKS) / log(2));
  // uint32_t n_div_log2_blocks = (((1 << logn_shmem) >> (log2_num_blocks + 1)) - 1);
  // uint32_t num_blocks2x = NUM_BLOCKS * 2; // TODO: ? uint32_t

  // ntt_template_kernel_shared_rev<<<NUM_BLOCKS, NUM_THREADS, shared_mem, 0>>>(d_inout, 1 << logn_shmem, d_twf, n / 2, total_tasks, 0, logn_shmem - 1, n_div_log2_blocks, num_blocks2x, (1 << logn_shmem) - 1);

  return 0;
}

template <typename E, typename S>
int ntt_batch_bc_template(E *d_inout, S *d_twf, unsigned n, unsigned batch_size, bool r, bool t)
{
  uint32_t logn = uint32_t(log(n) / log(2));

  int NUM_THREADS = min(n / 2, MAX_THREADS_BATCH);
  int chunks = max(int((n / 2) / NUM_THREADS), 1);
  int total_tasks = batch_size * chunks;
  int NUM_BLOCKS = total_tasks;
  int max_sharedmem = 512 * sizeof(E);
  int shared_mem = (2 * NUM_THREADS) * sizeof(E); // TODO: calculator, as shared mem size may be more efficient less then max to allow more concurrent blocks on SM
  uint32_t logn_shmem = uint32_t(log(2 * NUM_THREADS) / log(2));
  if (r)
  {
    // for (uint32_t s = logn - 1; s >= logn_shmem; s--) // TODO: this loop also can be unrolled
    // for (uint32_t s = 0; s < logn; s++) // TODO: this loop also can be unrolled
    for (uint32_t s = logn; s > 0; s--) // TODO: this loop also can be unrolled
    {
      ntt_template_kernel_bc<<<NUM_BLOCKS, NUM_THREADS>>>(d_inout, n, d_twf, n, total_tasks, s - 1, t);
    }
  }
  else
  {
    // for (uint32_t s = 0; s < logn; s++) // TODO: this loop also can be unrolled
    for (uint32_t s = 0; s < logn; s++) // TODO: this loop also can be unrolled
    {
      ntt_template_kernel_bc<<<NUM_BLOCKS, NUM_THREADS>>>(d_inout, n, d_twf, n, total_tasks, s, t);
    }
  }

  // uint32_t log2_num_blocks = (log(NUM_BLOCKS) / log(2));
  // uint32_t n_div_log2_blocks = (((1 << logn_shmem) >> (log2_num_blocks + 1)) - 1);
  // uint32_t num_blocks2x = NUM_BLOCKS * 2; // TODO: ? uint32_t

  // ntt_template_kernel_shared_rev<<<NUM_BLOCKS, NUM_THREADS, shared_mem, 0>>>(d_inout, 1 << logn_shmem, d_twf, n / 2, total_tasks, 0, logn_shmem - 1, n_div_log2_blocks, num_blocks2x, (1 << logn_shmem) - 1);

  return 0;
}

template <typename S>
int ntt_batch_bc(S *d_inout, S *d_twf, unsigned n, unsigned batch_size, bool r, bool t)
{
  return ntt_batch_bc_template(d_inout, d_twf, n, batch_size, r, t);
}

template <typename S>
int ntt_batch(S *d_inout, S *d_twf, unsigned n, unsigned batch_size)
{
  return ntt_batch_template(d_inout, d_twf, n, batch_size);
}

template <typename S>
int bailey_ntt(S *d_inout, S *d_twf, S *d_full_twf, unsigned n, unsigned batch_size, bool r1, bool t1, bool tt1, bool r2, bool t2, bool tt2)
{
  uint32_t logn = uint32_t(log(n) / log(2));

  dim3 threads(TILE_DIM, BLOCK_ROWS);
  dim3 blocks(batch_size / TILE_DIM, n / TILE_DIM);

  ntt_batch_bc(d_inout, d_twf, n, batch_size, r1, t1);

  transpose<<<blocks, threads>>>(d_inout);

  batch_mul_tw_ij<<<batch_size, n>>>(d_inout, d_full_twf, n, batch_size, logn, tt1, tt2);
  ntt_batch_bc(d_inout, d_twf, n, batch_size, r2, t2);

  transpose<<<blocks, threads>>>(d_inout);

  return 0;
}
///

/**
 * Evaluate a polynomial on a coset.
 * Note: this function does not preform any bit-reverse permutations on its inputs or outputs, so the order of outputs is bit-reversed.
 * @param d_out The evaluations of the polynomial on coset `u` * `d_domain`.
 * @param d_coefficients Input array of coefficients of a polynomial of type E (elements).
 * @param d_domain Domain on which the polynomial is evaluated (see `coset` flag). Must be a subgroup.
 * @param domain_size Length of `d_domain` array, on which the polynomial is computed.
 * @param n The number of coefficients, which might be different from `domain_size`.
 * @param coset The flag that indicates whether to evaluate on a coset. If false, evaluate on a subgroup `d_domain`.
 * @param coset_powers If `coset` is true, a list of powers `[1, u, u^2, ..., u^{n-1}]` where `u` is the generator of the coset.
 */
template <typename E, typename S>
int evaluate(E *d_out, E *d_coefficients, S *d_domain, unsigned domain_size, unsigned n, bool coset, S *coset_powers)
{
  return evaluate_batch<E, S>(d_out, d_coefficients, d_domain, domain_size, n, 1, coset, coset_powers);
}

template <typename S>
int interpolate_scalars(S *d_out, S *d_evaluations, S *d_domain, unsigned n)
{
  return interpolate(d_out, d_evaluations, d_domain, n);
}

template <typename S>
int interpolate_scalars_batch(S *d_out, S *d_evaluations, S *d_domain, unsigned n, unsigned batch_size)
{
  return interpolate_batch(d_out, d_evaluations, d_domain, n, batch_size);
}

template <typename E, typename S>
int interpolate_points(E *d_out, E *d_evaluations, S *d_domain, unsigned n)
{
  return interpolate(d_out, d_evaluations, d_domain, n);
}

template <typename E, typename S>
int interpolate_points_batch(E *d_out, E *d_evaluations, S *d_domain, unsigned n, unsigned batch_size)
{
  return interpolate_batch(d_out, d_evaluations, d_domain, n, batch_size);
}

template <typename S>
int evaluate_scalars(S *d_out, S *d_coefficients, S *d_domain, unsigned domain_size, unsigned n)
{
  S *_null = nullptr;
  return evaluate(d_out, d_coefficients, d_domain, domain_size, n, false, _null);
}

template <typename S>
int evaluate_scalars_batch(S *d_out, S *d_coefficients, S *d_domain, unsigned domain_size, unsigned n, unsigned batch_size)
{
  S *_null = nullptr;
  return evaluate_batch(d_out, d_coefficients, d_domain, domain_size, n, batch_size, false, _null);
}

template <typename E, typename S>
int evaluate_points(E *d_out, E *d_coefficients, S *d_domain, unsigned domain_size, unsigned n)
{
  S *_null = nullptr;
  return evaluate(d_out, d_coefficients, d_domain, domain_size, n, false, _null);
}

template <typename E, typename S>
int evaluate_points_batch(E *d_out, E *d_coefficients, S *d_domain,
                          unsigned domain_size, unsigned n, unsigned batch_size)
{
  S *_null = nullptr;
  return evaluate_batch(d_out, d_coefficients, d_domain, domain_size, n, batch_size, false, _null);
}

template <typename S>
int evaluate_scalars_on_coset(S *d_out, S *d_coefficients, S *d_domain,
                              unsigned domain_size, unsigned n, S *coset_powers)
{
  return evaluate(d_out, d_coefficients, d_domain, domain_size, n, true, coset_powers);
}

template <typename E, typename S>
int evaluate_scalars_on_coset_batch(S *d_out, S *d_coefficients, S *d_domain, unsigned domain_size,
                                    unsigned n, unsigned batch_size, S *coset_powers)
{
  return evaluate_batch(d_out, d_coefficients, d_domain, domain_size, n, batch_size, true, coset_powers);
}

template <typename E, typename S>
int evaluate_points_on_coset(E *d_out, E *d_coefficients, S *d_domain,
                             unsigned domain_size, unsigned n, S *coset_powers)
{
  return evaluate(d_out, d_coefficients, d_domain, domain_size, n, true, coset_powers);
}

template <typename E, typename S>
int evaluate_points_on_coset_batch(E *d_out, E *d_coefficients, S *d_domain, unsigned domain_size,
                                   unsigned n, unsigned batch_size, S *coset_powers)
{
  return evaluate_batch(d_out, d_coefficients, d_domain, domain_size, n, batch_size, true, coset_powers);
}

#endif