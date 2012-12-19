#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <cublas_v2.h>


/* Helper function */

/* Function to call cudaSetDevice */
extern "C" int set_gpu(const int *device_id)
{
  cudaError_t cuda_status;

  cuda_status = cudaSetDevice(*device_id);

  if (cuda_status != cudaSuccess)
    {
      printf("Failed to set device!\n");
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

/* Function to call cudaDeviceSynchronize */
extern "C" int sync_gpu()
{
  cudaError_t cuda_status;  

  cuda_status = cudaDeviceSynchronize();
  if (cuda_status != cudaSuccess)
    {
      printf("Streams failed to synchronise!\n");
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
