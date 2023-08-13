#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include "tr.cuh"
#define uS_PER_SEC 1000000
#define uS_PER_mS 1000
#define N 1024
#define M 1024

template <typename T>
__global__ void transposeCoalesced(T *odata, const T *idata)
{
  __shared__ T tile[TILE_DIM][TILE_DIM + 1];

  int x = blockIdx.x * TILE_DIM + threadIdx.x;
  int y = blockIdx.y * TILE_DIM + threadIdx.y;
  int width = gridDim.x * TILE_DIM;

  for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
    tile[threadIdx.y + j][threadIdx.x] = idata[(y + j) * width + x];

  __syncthreads();

  x = blockIdx.y * TILE_DIM + threadIdx.x; // transpose block offset
  y = blockIdx.x * TILE_DIM + threadIdx.y;

  for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
    odata[(y + j) * width + x] = tile[threadIdx.x][threadIdx.y + j];
}

template <typename T>
int validate(const T *mat, const T *mat_t, int n, int m)
{
  int result = 1;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      if (mat[(i * m) + j] != mat_t[(j * n) + i])
        result = 0;
  return result;
}

int main()
{

  timeval t1, t2;
  float *matrix = (float *)malloc(N * M * sizeof(float));
  for (int i = 0; i < N; i++)
    for (int j = 0; j < M; j++)
      matrix[(i * M) + j] = i;
  // Starting the timer
  gettimeofday(&t1, NULL);
  float *matrixT = (float *)malloc(N * M * sizeof(float));
  for (int i = 0; i < N; i++)
    for (int j = 0; j < M; j++)
      matrixT[(j * N) + i] = matrix[(i * M) + j]; // matrix is obviously filled
                                                  // Ending the timer
  gettimeofday(&t2, NULL);
  if (!validate(matrix, matrixT, N, M))
  {
    printf("fail!\n");
    return 1;
  }
  float et1 = (((t2.tv_sec * uS_PER_SEC) + t2.tv_usec) - ((t1.tv_sec * uS_PER_SEC) + t1.tv_usec)) / (float)uS_PER_mS;
  printf("CPU time = %fms\n", et1);

  float *h_matrixT, *d_matrixT, *d_matrix;
  h_matrixT = (float *)(malloc(N * M * sizeof(float)));
  cudaMalloc((void **)&d_matrixT, N * M * sizeof(float));
  cudaMalloc((void **)&d_matrix, N * M * sizeof(float));
  cudaMemcpy(d_matrix, matrix, N * M * sizeof(float), cudaMemcpyHostToDevice);

  // Starting the timer
  gettimeofday(&t1, NULL);

  const float alpha = 1.0;
  const float beta = 0.0;
  cublasHandle_t handle;
  // gettimeofday(&t1, NULL);
  cublasCreate(&handle);
  gettimeofday(&t1, NULL);
  cublasSgeam(handle, CUBLAS_OP_T, CUBLAS_OP_N, N, M, &alpha, d_matrix, M, &beta, d_matrix, N, d_matrixT, N);
  cudaDeviceSynchronize();
  gettimeofday(&t2, NULL);
  cublasDestroy(handle);

  // Ending the timer
  float et2 = (((t2.tv_sec * uS_PER_SEC) + t2.tv_usec) - ((t1.tv_sec * uS_PER_SEC) + t1.tv_usec)) / (float)uS_PER_mS;
  printf("GPU Sgeam time = %fms\n", et2);

  cudaMemcpy(h_matrixT, d_matrixT, N * M * sizeof(float), cudaMemcpyDeviceToHost);
  if (!validate(matrix, h_matrixT, N, M))
  {
    printf("fail!\n");
    return 1;
  }
  cudaMemset(d_matrixT, 0, N * M * sizeof(float));
  memset(h_matrixT, 0, N * M * sizeof(float));
  dim3 threads(TILE_DIM, BLOCK_ROWS);
  dim3 blocks(N / TILE_DIM, M / TILE_DIM);
  gettimeofday(&t1, NULL);
  transposeCoalesced<<<blocks, threads>>>(d_matrixT, d_matrix);
  cudaDeviceSynchronize();
  gettimeofday(&t2, NULL);
  cudaMemcpy(h_matrixT, d_matrixT, N * M * sizeof(float), cudaMemcpyDeviceToHost);
  if (!validate(matrix, h_matrixT, N, M))
  {
    printf("fail!\n");
    return 1;
  }
  float et3 = (((t2.tv_sec * uS_PER_SEC) + t2.tv_usec) - ((t1.tv_sec * uS_PER_SEC) + t1.tv_usec)) / (float)uS_PER_mS;
  printf("GPU kernel time = %fms\n", et3);

  memset(h_matrixT, 0, N * M * sizeof(float));
  gettimeofday(&t1, NULL);
  // dim3 threads(TILE_DIM, BLOCK_ROWS);
  // dim3 blocks(N / TILE_DIM, M / TILE_DIM);
  transpose<<<blocks, threads>>>(d_matrix);
  cudaDeviceSynchronize();
  gettimeofday(&t2, NULL);
  cudaMemcpy(h_matrixT, d_matrix, N * M * sizeof(float), cudaMemcpyDeviceToHost);
  if (!validate(matrix, h_matrixT, N, M))
  {
    printf("fail!\n");
    return 1;
  }
  float et4 = (((t2.tv_sec * uS_PER_SEC) + t2.tv_usec) - ((t1.tv_sec * uS_PER_SEC) + t1.tv_usec)) / (float)uS_PER_mS;
  printf("GPU in-place kernel time = %fms\n", et4);

  cudaFree(d_matrix);
  cudaFree(d_matrixT);
  return 0;
}
// $ nvcc -arch=sm_86 -o tr tr.cu -lcublas; ./tr
// $ ./t469
// CPU time = 450.095001ms
// GPU Sgeam time = 1.937000ms
// GPU kernel time = 1.694000ms
// GPU in-place kernel time = 1.839000ms
// $