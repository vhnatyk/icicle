#pragma once
#ifndef NTT_H
#define NTT_H

#include "../../utils/sharedmem.cuh"

const uint32_t MAX_NUM_THREADS = 1024;
const uint32_t MAX_THREADS_BATCH = 512;

/**
 * Computes the twiddle factors.
 * Outputs: d_twiddles[i] = omega^i.
 * @param d_twiddles input empty array.
 * @param n_twiddles number of twiddle factors.
 * @param omega multiplying factor.
 */
template <typename S>
__global__ void twiddle_factors_kernel(S *d_twiddles, uint32_t n_twiddles, S omega)
{
  for (uint32_t i = 0; i < n_twiddles; i++)
  {
    d_twiddles[i] = S::zero();
  }
  d_twiddles[0] = S::one();
  for (uint32_t i = 0; i < n_twiddles - 1; i++)
  {
    d_twiddles[i + 1] = omega * d_twiddles[i];
  }
}

/**
 * Fills twiddles array with twiddle factors.
 * @param twiddles input empty array.
 * @param n_twiddles number of twiddle factors.
 * @param omega multiplying factor.
 */
template <typename S>
S *fill_twiddle_factors_array(uint32_t n_twiddles, S omega)
{
  size_t size_twiddles = n_twiddles * sizeof(S);
  S *d_twiddles;
  cudaMalloc(&d_twiddles, size_twiddles);
  twiddle_factors_kernel<S><<<1, 1>>>(d_twiddles, n_twiddles, omega);
  return d_twiddles;
}

/**
 * Returns the bit reversed order of a number.
 * for example: on inputs num = 6 (110 in binary) and logn = 3
 * the function should return 3 (011 in binary.)
 * @param num some number with bit representation of size logn.
 * @param logn length of bit representation of `num`.
 * @return bit reveresed order or `num`.
 */
__device__ __host__ uint32_t reverseBits(uint32_t num, uint32_t logn)
{
  unsigned int reverse_num = 0;
  for (int i = 0; i < logn; i++)
  {
    if ((num & (1 << i)))
      reverse_num |= 1 << ((logn - 1) - i);
  }
  return reverse_num;
}

/**
 * Returns the bit reversal ordering of the input array.
 * for example: on input ([a[0],a[1],a[2],a[3]], 4, 2) it returns
 * [a[0],a[3],a[2],a[1]] (elements in indices 3,1 swhich places).
 * @param arr array of some object of type T of size which is a power of 2.
 * @param n length of `arr`.
 * @param logn log(n).
 * @return A new array which is the bit reversed version of input array.
 */
template <typename T>
T *template_reverse_order(T *arr, uint32_t n, uint32_t logn)
{
  T *arrReversed = new T[n];
  for (uint32_t i = 0; i < n; i++)
  {
    uint32_t reversed = reverseBits(i, logn);
    arrReversed[i] = arr[reversed];
  }
  return arrReversed;
}

template <typename T>
__global__ void reverse_order_kernel(T *arr, T *arr_reversed, uint32_t n, uint32_t logn, uint32_t batch_size)
{
  int threadId = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (threadId < n * batch_size)
  {
    int idx = threadId % n;
    int batch_idx = threadId / n;
    int idx_reversed = __brev(idx) >> (32 - logn);
    arr_reversed[batch_idx * n + idx_reversed] = arr[batch_idx * n + idx];
  }
}

/**
 * Bit-reverses a batch of input arrays in-place inside GPU.
 * for example: on input array ([a[0],a[1],a[2],a[3]], 4, 2) it returns
 * [a[0],a[3],a[2],a[1]] (elements at indices 3 and 1 swhich places).
 * @param arr batch of arrays of some object of type T. Should be on GPU.
 * @param n length of `arr`.
 * @param logn log(n).
 * @param batch_size the size of the batch.
 */
template <typename T>
void reverse_order_batch(T *arr, uint32_t n, uint32_t logn, uint32_t batch_size)
{
  T *arr_reversed;
  cudaMalloc(&arr_reversed, n * batch_size * sizeof(T));
  int number_of_threads = MAX_THREADS_BATCH;
  int number_of_blocks = (n * batch_size + number_of_threads - 1) / number_of_threads;
  reverse_order_kernel<<<number_of_blocks, number_of_threads>>>(arr, arr_reversed, n, logn, batch_size);
  cudaMemcpy(arr, arr_reversed, n * batch_size * sizeof(T), cudaMemcpyDeviceToDevice);
  cudaFree(arr_reversed);
}

/**
 * Bit-reverses an input array in-place inside GPU.
 * for example: on array ([a[0],a[1],a[2],a[3]], 4, 2) it returns
 * [a[0],a[3],a[2],a[1]] (elements at indices 3 and 1 swhich places).
 * @param arr array of some object of type T of size which is a power of 2. Should be on GPU.
 * @param n length of `arr`.
 * @param logn log(n).
 */
template <typename T>
void reverse_order(T *arr, uint32_t n, uint32_t logn)
{
  reverse_order_batch(arr, n, logn, 1);
}

// /**
//  * Cooley-Tukey butterfly kernel.
//  * @param arr array of objects of type E (elements).
//  * @param twiddles array of twiddle factors of type S (scalars).
//  * @param n size of arr.
//  * @param n_twiddles size of omegas.
//  * @param m "pair distance" - indicate distance of butterflies inputs.
//  * @param i Cooley-Tukey FFT stage number.
//  * @param max_thread_num maximal number of threads in stage.
//  */
// template <typename E, typename S>
// __global__ void template_butterfly_kernel(E *arr, S *twiddles, uint32_t n, uint32_t n_twiddles, uint32_t m, uint32_t i, uint32_t max_thread_num)
// {
//   int j = (blockIdx.x * blockDim.x) + threadIdx.x;
//   if (j < max_thread_num)
//   {
//     uint32_t g = j * (n / m);
//     uint32_t k = i + j + (m >> 1);
//     E u = arr[i + j];
//     E v = twiddles[g * n_twiddles / n] * arr[k];
//     arr[i + j] = u + v;
//     arr[k] = u - v;
//   }
// }

/**
 * Multiply the elements of an input array by a scalar in-place.
 * @param arr input array.
 * @param n size of arr.
 * @param n_inv scalar of type S (scalar).
 */
template <typename E, typename S>
__global__ void template_normalize_kernel(E *arr, uint32_t n, S scalar)
{
  int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (tid < n)
  {
    arr[tid] = scalar * arr[tid];
  }
}

// /**
//  * Cooley-Tukey NTT.
//  * NOTE! this function assumes that d_arr and d_twiddles are located in the device memory.
//  * @param d_arr input array of type E (elements) allocated on the device memory.
//  * @param n length of d_arr.
//  * @param logn log(n).
//  * @param d_twiddles twiddle factors of type S (scalars) array allocated on the device memory (must be a power of 2).
//  * @param n_twiddles length of d_twiddles.
//  */
// template <typename E, typename S>
// void template_ntt_on_device_memory(E *d_arr, uint32_t n, uint32_t logn, S *d_twiddles, uint32_t n_twiddles)
// {
//   uint32_t m = 2;
//   for (uint32_t s = 0; s < logn; s++)
//   {
//     for (uint32_t i = 0; i < n; i += m)
//     {
//       int shifted_m = m >> 1;
//       int number_of_threads = MAX_NUM_THREADS ^ ((shifted_m ^ MAX_NUM_THREADS) & -(shifted_m < MAX_NUM_THREADS));
//       int number_of_blocks = shifted_m / MAX_NUM_THREADS + 1;
//       template_butterfly_kernel<E, S><<<number_of_threads, number_of_blocks>>>(d_arr, d_twiddles, n, n_twiddles, m, i, m >> 1);
//     }
//     m <<= 1;
//   }
// }

// /**
//  * Cooley-Tukey NTT.
//  * NOTE! this function assumes that d_twiddles are located in the device memory.
//  * @param arr input array of type E (elements).
//  * @param n length of d_arr.
//  * @param d_twiddles twiddle factors of type S (scalars) array allocated on the device memory (must be a power of 2).
//  * @param n_twiddles length of d_twiddles.
//  * @param inverse indicate if the result array should be normalized by n^(-1).
//  */
// template <typename E, typename S>
// E *ntt_template(E *arr, uint32_t n, S *d_twiddles, uint32_t n_twiddles, bool inverse)
// {
//   uint32_t logn = uint32_t(log(n) / log(2));
//   size_t size_E = n * sizeof(E);
//   E *arrReversed = template_reverse_order<E>(arr, n, logn);
//   E *d_arrReversed;
//   cudaMalloc(&d_arrReversed, size_E);
//   cudaMemcpy(d_arrReversed, arrReversed, size_E, cudaMemcpyHostToDevice);
//   template_ntt_on_device_memory<E, S>(d_arrReversed, n, logn, d_twiddles, n_twiddles);
//   if (inverse)
//   {
//     int NUM_THREADS = MAX_NUM_THREADS;
//     int NUM_BLOCKS = (n + NUM_THREADS - 1) / NUM_THREADS;
//     template_normalize_kernel<E, S><<<NUM_THREADS, NUM_BLOCKS>>>(d_arrReversed, n, S::inv_log_size(logn));
//   }
//   cudaMemcpy(arrReversed, d_arrReversed, size_E, cudaMemcpyDeviceToHost);
//   cudaFree(d_arrReversed);
//   return arrReversed;
// }

/**
 * Cooley-Tukey (scalar) NTT.
 * @param arr input array of type E (element).
 * @param n length of d_arr.
 * @param inverse indicate if the result array should be normalized by n^(-1).
 */
template <typename E, typename S>
uint32_t ntt_end2end_template(E *arr, uint32_t n, bool inverse)
{
  // uint32_t logn = uint32_t(log(n) / log(2));
  // uint32_t n_twiddles = n;
  // S *twiddles = new S[n_twiddles];
  // S *d_twiddles;
  // if (inverse)
  // {
  //   d_twiddles = fill_twiddle_factors_array(n_twiddles, S::omega_inv(logn));
  // }
  // else
  // {
  //   d_twiddles = fill_twiddle_factors_array(n_twiddles, S::omega(logn));
  // }
  // E *result = ntt_template<E, S>(arr, n, d_twiddles, n_twiddles, inverse);
  // for (int i = 0; i < n; i++)
  // {
  //   arr[i] = result[i];
  // }
  // cudaFree(d_twiddles);
  return 0; // TODO add
}

/**
 * Returens the bit reversal ordering of the input array according to the batches *in place*.
 * The assumption is that arr is divided into N tasks of size n.
 * Tasks indicates the index of the task (out of N).
 * @param arr input array of type T.
 * @param n length of arr.
 * @param logn log(n).
 * @param task log(n).
 */
template <typename T>
__device__ __host__ void reverseOrder_batch(T *arr, uint32_t n, uint32_t logn, uint32_t task)
{
  for (uint32_t i = 0; i < n; i++)
  {
    uint32_t reversed = reverseBits(i, logn);
    if (reversed > i)
    {
      T tmp = arr[task * n + i];
      arr[task * n + i] = arr[task * n + reversed];
      arr[task * n + reversed] = tmp;
    }
  }
}

/**
 * Cooley-Tukey butterfly kernel.
 * @param arr array of objects of type E (elements).
 * @param twiddles array of twiddle factors of type S (scalars).
 * @param n size of arr.
 * @param n_twiddles size of omegas.
 * @param m "pair distance" - indicate distance of butterflies inputs.
 * @param i Cooley-TUckey FFT stage number.
 * @param offset offset corr. to the specific taks (in batch).
 */
template <typename E, typename S>
__device__ __host__ void butterfly(E *arrReversed, S *omegas, uint32_t n, uint32_t n_omegas, uint32_t m, uint32_t i, uint32_t j, uint32_t offset)
{
  uint32_t g = j * (n / m);
  uint32_t k = i + j + (m >> 1);
  E u = arrReversed[offset + i + j];
  E v = omegas[g * n_omegas / n] * arrReversed[offset + k];
  arrReversed[offset + i + j] = u + v;
  arrReversed[offset + k] = u - v;
}

template <typename E, typename S>
__global__ void batch_mul_tw_ij(E *element_vec, S *scalar_vec, unsigned n2, unsigned n1, unsigned logn)
{
  int j = threadIdx.x;
  int i = blockIdx.x;
  int tid = blockDim.x * i + j;
  
  if (tid < n2 * n1)
  {
    element_vec[tid] = scalar_vec[(__brev(j) >> (32 - logn)) * (__brev(i) >> (32 - logn))] * element_vec[tid];
    //     element_vec[i * n1 + j] = scalar_vec[i * j] * element_vec[i * n1 + j];
  }

  // for (int j = 0; j < n2; j++)
  // {
  //   for (int i = 0; i < n1; i++)
  //   {
  //     element_vec[i * n1 + j] = scalar_vec[i * j] * element_vec[i * n1 + j];
  //   }
  // }
}

//************************************************************************************************
/**
 * Cooley-Tuckey NTT.
 * NOTE! this function assumes that d_twiddles are located in the device memory.
 * @param arr input array of type E (elements).
 * @param n length of d_arr.
 * @param twiddles twiddle factors of type S (scalars) array allocated on the device memory (must be a power of 2).
 * @param n_twiddles length of twiddles.
 * @param max_task max count of parallel tasks.
 * @param s log2(n) loop index.
 */
template <typename E, typename S>
__launch_bounds__(MAX_THREADS_BATCH, 3)
    __global__ void ntt_template_kernel_shared_rev(E *__restrict__ arr_g, uint32_t n, const S *__restrict__ r_twiddles, uint32_t n_div_2_twiddles, uint32_t max_task, uint8_t ss, uint8_t logn_m_1, uint32_t n_div_log2_blocks, uint32_t num_blocks2x, uint32_t n_m1)
{
  if (blockIdx.x < max_task)
  {
    // flattened loop allows parallel processing

    if (threadIdx.x < blockDim.x)
    {
      SharedMemory<E> smem;
      E *arr = smem.getPointer();

      uint16_t l = (blockIdx.x & n_div_log2_blocks) * blockDim.x + threadIdx.x; // to l from chunks to full
      // uint32_t offset = blockIdx.x * num_blocks2x;
      arr_g += blockIdx.x * (blockDim.x * 2);

      uint8_t s = logn_m_1 - ss;
      uint16_t shift_s = 1 << s;
      uint16_t j = l & (shift_s - 1); // Equivalent to: l % (1 << s)
      auto tw = r_twiddles[j * (n_div_2_twiddles >> s)];
      uint16_t oij = (((l >> s) * (shift_s << 1)) & n_m1) + j;
      j = oij + shift_s; // reuse for k

      E *uu = arr_g + oij;
      E *vv = arr_g + j;
      E *uuu = arr + oij;
      E *vvv = arr + j;

      auto u = *uu;
      auto v = *vv;
      *uuu = u + v;
      *vvv = tw * (u - v);

      ss++;
#pragma unroll 7
      for (; ss < logn_m_1 - 1; ss++)
      {
        s = logn_m_1 - ss;
        if (s > 4)
          __syncthreads();

        shift_s = 1 << s;
        j = l & (shift_s - 1); // Equivalent to: l % (1 << s)
        tw = r_twiddles[j * (n_div_2_twiddles >> s)];
        oij = (((l >> s) * (shift_s << 1)) & n_m1) + j;
        j = oij + shift_s; // reuse for k

        uu = arr + oij;
        vv = arr + j;

        uuu = uu;
        vvv = vv;

        u = *uu;
        v = *vv;
        *uuu = u + v;
        *vvv = tw * (u - v);

        //__syncthreads();
      }

      s = 1;
      shift_s = 2;
      j = l & 1; // Equivalent to: l % (1 << s)
      tw = r_twiddles[j * (n_div_2_twiddles >> 1)];
      // if (j != 0)
      //   tw = r_twiddles[n_div_2_twiddles >> 1];
      oij = (((l >> 1) * 4) & n_m1) + j;
      shift_s = oij + shift_s; // reuse for k

      uu = arr + oij;
      vv = arr + oij + 2;

      uuu = uu;
      vvv = vv;

      u = *uu;
      v = *vv;
      *uuu = u + v;
      *vvv = j == 0 ? (u - v) : tw * (u - v);
      // if (blockIdx.x != 0 && threadIdx.x != 0)
      // {
      //   tw = r_twiddles[blockIdx.x];
      //   tw2 = tw * tw2;
      //   tw = tw * r_twiddles[threadIdx.x * 2];
      //   //   // tw = r_twiddles[blockIdx.x];
      //   //   // tw2 = tw * r_twiddles[threadIdx.x * 2 + 1];
      //   //   // tw = tw * r_twiddles[threadIdx.x * 2];

      //   //   // tw = r_twiddles[blockIdx.x] * r_twiddles[threadIdx.x * 2];
      //   //   // tw2 = r_twiddles[blockIdx.x] * r_twiddles[threadIdx.x * 2 + 1];
      // }
      ////////
      // s = 0;
      // shift_s = 1;
      // j = l & (shift_s - 1); // Equivalent to: l % (1 << s)
      //  tw = r_twiddles[j * (n_twiddles >> s)];         // TODO: it all can be constant here except oij

      oij = ((l << 1) & n_m1); // but the simplification oij = (l >> 1) & n_m1
                               // actually breaks correctness and decreases performance?!!
      j = oij + 1;             // reuse for k

      uu = arr + oij;
      vv = arr + j;

      uuu = arr_g + oij;
      vvv = arr_g + j;

      u = *uu;
      v = *vv;
      // if (blockIdx.x == 0 || threadIdx.x == 0)
      // {
      *uuu = u + v;
      *vvv = u - v; // s is 0 here - so twiddle factor is 1 - so no need to multiply
      // }
      // else
      // {
      //   *uuu = tw * (u + v);
      //   *vvv = tw2 * (u - v); // s is 0 here - so twiddle factor is 1 - so no need to multiply
      // }
    }
  }
}
//************************************************************************************************

//************************************************************************************************
/**
 * Cooley-Tuckey NTT.
 * NOTE! this function assumes that d_twiddles are located in the device memory.
 * @param arr input array of type E (elements).
 * @param n length of d_arr.
 * @param twiddles twiddle factors of type S (scalars) array allocated on the device memory (must be a power of 2).
 * @param n_twiddles length of twiddles.
 * @param max_task max count of parallel tasks.
 * @param s log2(n) loop index.
 */
template <typename E, typename S>
__global__ void ntt_template_kernel_shared(E *__restrict__ arr_g, uint32_t n, const S *__restrict__ r_twiddles, uint32_t n_twiddles, uint32_t max_task, uint32_t s, uint32_t logn, bool rev)
{
  SharedMemory<E> smem;
  E *arr = smem.getPointer();

  uint32_t task = blockIdx.x;
  uint32_t loop_limit = blockDim.x;
  uint32_t chunks = n / (loop_limit * 2);
  uint32_t offset = (task / chunks) * n;
  if (task < max_task)
  {
    // flattened loop allows parallel processing
    uint32_t l = threadIdx.x;

    if (l < loop_limit)
    {
      // #pragma unroll 8
      for (; s < logn; s++) // TODO: this loop also can be unrolled
      {
        // if (s == 0) //this actually can be faster even by introducing extra read and (see below)...
        // {
        //   arr[l] = arr_g[offset + l];
        //   arr[loop_limit + l] = arr_g[offset + loop_limit + l];
        //   __syncthreads();
        // }

        uint32_t ntw_i = task % chunks;

        uint32_t n_twiddles_div = n_twiddles >> (s + 1);

        uint32_t shift_s = 1 << s;
        uint32_t shift2_s = 1 << (s + 1);

        l = ntw_i * loop_limit + l; // to l from chunks to full

        uint32_t j = l & (shift_s - 1);               // Equivalent to: l % (1 << s)
        uint32_t i = ((l >> s) * shift2_s) & (n - 1); // (..) % n (assuming n is power of 2)
        uint32_t oij = i + j;
        uint32_t k = oij + shift_s;

        E u = s == 0 ? arr_g[offset + oij] : arr[oij];
        E v = rev ? (s == 0 ? arr_g[offset + k] : arr[k]) : (r_twiddles[j * n_twiddles_div] * (s == 0 ? arr_g[offset + k] : arr[k]));
        if (s == (logn - 1))
        {
          arr_g[offset + oij] = u + v;
          if (rev)
            arr_g[offset + k] = r_twiddles[j * n_twiddles_div] * (u - v);
          else
            arr_g[offset + k] = u - v;
        }
        else
        {
          arr[oij] = u + v;
          arr[k] = u - v;
          if (rev)
            arr[k] = r_twiddles[j * n_twiddles_div] * arr[k];
        }

        __syncthreads();
      }

      // ... and extra write
      // arr_g[offset + l] = arr[l];
      // arr_g[offset + loop_limit + l] = arr[l + loop_limit];
      // __syncthreads();
    }
  }
}
//************************************************************************************************

/**
 * Cooley-Tukey NTT.
 * NOTE! this function assumes that d_twiddles are located in the device memory.
 * @param arr input array of type E (elements).
 * @param n length of d_arr.
 * @param twiddles twiddle factors of type S (scalars) array allocated on the device memory (must be a power of 2).
 * @param n_twiddles length of twiddles.
 * @param max_task max count of parallel tasks.
 * @param s log2(n) loop index.
 */
template <typename E, typename S>
__global__ void ntt_template_kernel(E *arr, uint32_t n, S *twiddles, uint32_t n_twiddles, uint32_t max_task, uint32_t s, bool rev)
{
  int task = blockIdx.x;
  int chunks = n / (blockDim.x * 2);

  if (task < max_task)
  {
    // flattened loop allows parallel processing
    uint32_t l = threadIdx.x;
    uint32_t loop_limit = blockDim.x;

    if (l < loop_limit)
    {
      uint32_t ntw_i = task % chunks;

      uint32_t shift_s = 1 << s;
      uint32_t shift2_s = 1 << (s + 1);
      uint32_t n_twiddles_div = n_twiddles >> (s + 1);

      l = ntw_i * blockDim.x + l; // to l from chunks to full

      uint32_t j = l & (shift_s - 1);               // Equivalent to: l % (1 << s)
      uint32_t i = ((l >> s) * shift2_s) & (n - 1); // (..) % n (assuming n is power of 2)
      uint32_t k = i + j + shift_s;

      uint32_t offset = (task / chunks) * n;
      E u = arr[offset + i + j];
      E v = rev ? arr[offset + k] : twiddles[j * n_twiddles_div] * arr[offset + k];
      arr[offset + i + j] = u + v;
      arr[offset + k] = u - v;
      if (rev)
        arr[offset + k] = twiddles[j * n_twiddles_div] * arr[offset + k];
    }
  }
}

template <typename E, typename S>
HOST_DEVICE_INLINE void butterfly_bc(E* u, E* v, S* w, bool rev) {
    if (rev) {
        // Input RBO, Output RBO
        E diff = *u - *v;
        *u = *u + *v;
        *v = diff * (*w);
    } else {
        // Input RBO, Output Natural
        // *v = *u - (*v * (*w));
        // *u = *u + *v;
        E vv = (*v * (*w));
        *v = *u - vv;
        *u = *u + vv;
    // } else if (!rev && output_rbo) {
    //     // // Input Natural, Output RBO
    //     // *u = *u + *v;
    //     // *v = *u - (*v * (*w));
    // } else {
    //     // // Input Natural, Output Natural
    //     // E sum = *u - *v;
    //     // *v = *u - *v;
    //     // *u = sum;
    }
}

/**
 * Cooley-Tukey NTT.
 * NOTE! this function assumes that d_twiddles are located in the device memory.
 * @param arr input array of type E (elements).
 * @param n length of d_arr.
 * @param twiddles twiddle factors of type S (scalars) array allocated on the device memory (must be a power of 2).
 * @param n_twiddles length of twiddles.
 * @param max_task max count of parallel tasks.
 * @param s log2(n) loop index.
 */
template <typename E, typename S>
__global__ void ntt_template_kernel_bc(E *arr, uint32_t n, S *twiddles, uint32_t n_twiddles, uint32_t max_task, uint32_t s, bool rev)
{
  int task = blockIdx.x;
  int chunks = n / (blockDim.x * 2);

  if (task < max_task)
  {
    // flattened loop allows parallel processing
    uint32_t l = threadIdx.x;
    uint32_t loop_limit = blockDim.x;

    if (l < loop_limit)
    {
      uint32_t ntw_i = task % chunks;

      uint32_t shift_s = 1 << s;
      uint32_t shift2_s = 1 << (s + 1);
      uint32_t n_twiddles_div = n_twiddles >> (s + 1);

      l = ntw_i * blockDim.x + l; // to l from chunks to full

      uint32_t j = l & (shift_s - 1);               // Equivalent to: l % (1 << s)
      uint32_t i = ((l >> s) * shift2_s) & (n - 1); // (..) % n (assuming n is power of 2)
      uint32_t k = i + j + shift_s;

      uint32_t offset = (task / chunks) * n;
      // E u = arr[offset + i + j];
      // E v = rev ? arr[offset + k] : twiddles[j * n_twiddles_div] * arr[offset + k];
      // arr[offset + i + j] = u + v;
      // arr[offset + k] = u - v;
      // if (rev)
      //   arr[offset + k] = twiddles[j * n_twiddles_div] * arr[offset + k];

      butterfly_bc(&arr[offset+ i + j], &arr[offset + k], &twiddles[j * n_twiddles_div], rev);
    }
  }
}

/**
 * Cooley-Tukey NTT.
 * NOTE! this function assumes that d_twiddles are located in the device memory.
 * @param arr input array of type E (elements).
 * @param n length of arr.
 * @param logn log2(n).
 * @param max_task max count of parallel tasks.
 */
template <typename E, typename S>
__global__ void ntt_template_kernel_rev_ord(E *arr, uint32_t n, uint32_t logn, uint32_t max_task)
{
  int task = (blockIdx.x * blockDim.x) + threadIdx.x;

  if (task < max_task)
  {
    reverseOrder_batch<E>(arr, n, logn, task);
  }
}

/**
 * Cooley-Tukey (scalar) NTT.
 * This is a bached version - meaning it assumes than the input array
 * consists of N arrays of size n. The function performs n-size NTT on each small array.
 * @param arr input array of type BLS12_381::scalar_t.
 * @param arr_size number of total elements = n * N.
 * @param n size of batch.
 * @param inverse indicate if the result array should be normalized by n^(-1).
 */
template <typename E, typename S>
uint32_t ntt_end2end_batch_template(E *arr, uint32_t arr_size, uint32_t n, bool inverse)
{
  // int batches = int(arr_size / n);
  // uint32_t logn = uint32_t(log(n) / log(2));
  // uint32_t n_twiddles = n; // n_twiddles is set to 4096 as BLS12_381::scalar_t::omega() is of that order.
  // size_t size_E = arr_size * sizeof(E);
  // S *d_twiddles;
  // if (inverse)
  // {
  //   d_twiddles = fill_twiddle_factors_array(n_twiddles, S::omega_inv(logn));
  // }
  // else
  // {
  //   d_twiddles = fill_twiddle_factors_array(n_twiddles, S::omega(logn));
  // }
  // E *d_arr;
  // cudaMalloc(&d_arr, size_E);
  // cudaMemcpy(d_arr, arr, size_E, cudaMemcpyHostToDevice);
  // E *arr_reversed;
  // cudaMalloc(&arr_reversed, n * batches * sizeof(E));
  // int number_of_threads = MAX_THREADS_BATCH;
  // int number_of_blocks = (arr_size + number_of_threads - 1) / number_of_threads;
  // // ntt_template_kernel_rev_ord<E, S><<<NUM_BLOCKS, NUM_THREADS>>>(d_arr, n, logn, batches);
  // reverse_order_kernel<<<number_of_blocks, number_of_threads>>>(d_arr, arr_reversed, n, logn, batches);
  // d_arr = arr_reversed;

  // int NUM_THREADS = min(n / 2, MAX_THREADS_BATCH);
  // int chunks = max(int((n / 2) / NUM_THREADS), 1);
  // int total_tasks = batches * chunks;
  // int NUM_BLOCKS = total_tasks;
  // int max_sharedmem = 512 * sizeof(E);
  // int shared_mem = 2 * NUM_THREADS * sizeof(E); // TODO: calculator, as shared mem size may be more efficient less then max to allow more concurrent blocks on SM
  // uint32_t logn_shmem = uint32_t(log(2 * NUM_THREADS) / log(2));
  // ntt_template_kernel_shared<<<NUM_BLOCKS, NUM_THREADS, shared_mem, 0>>>(d_arr, 1 << logn_shmem, d_twiddles, n, total_tasks, 0, logn_shmem, false);

  // for (uint32_t s = logn_shmem; s < logn; s++) // TODO: this loop also can be unrolled
  // {
  //   ntt_template_kernel<<<NUM_BLOCKS, NUM_THREADS>>>(d_arr, n, d_twiddles, n_twiddles, total_tasks, s, false);
  // }
  // if (inverse == true)
  // {
  //   NUM_THREADS = 64;
  //   NUM_BLOCKS = (arr_size + NUM_THREADS - 1) / NUM_THREADS;
  //   template_normalize_kernel<<<NUM_BLOCKS, NUM_THREADS>>>(d_arr, arr_size, S::inv_log_size(logn));
  // }
  // cudaMemcpy(arr, d_arr, size_E, cudaMemcpyDeviceToHost);
  // cudaFree(d_arr);
  // cudaFree(d_twiddles);
  return 0;
}

#endif