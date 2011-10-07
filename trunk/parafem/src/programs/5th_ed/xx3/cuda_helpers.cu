#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <cublas_v2.h>

/* Kernel definitions */

extern __shared__ float work_array[];

/* 
  This kernel assumes matrix is stored column wise

  This attempts 2D thread blocks
*/
__global__ void MultiMatVecMultiply3(int n_mat,
				     int n_row,
				     int n_col,
				     double* lhs_vector,
				     double* matrix,
				     double* rhs_vector)
{
  int global_row_id;
  int matrix_id;
  int row_id;
  int j;
  double tmp;

  /* Get the global index that corresponds to a row of some matrix */
  global_row_id =  threadIdx.x + blockIdx.x * blockDim.x;

  if (global_row_id < n_mat*n_row)
    {
      /* get matrix id */
      matrix_id = global_row_id/n_row;
      
      /* Get local row id */
      row_id = global_row_id%n_row;
      
      /* Variable to store partial row value */
      tmp = 0.0;
      
      /* Each thread loop through a part of the row */
      for (j=threadIdx.y; j<n_col; j+= blockDim.y)
	{
	  tmp +=
	    lhs_vector[j+matrix_id*n_col] *
	    matrix[row_id+j*n_row+matrix_id*n_col*n_row];
	}

      /* Put values in shared memory array */
      work_array[threadIdx.x + threadIdx.y*n_row] = tmp;

      /* Sync threads */
      __syncthreads();

      /* Reduce array values */
      j = threadIdx.y;

      while ( j>1 )
	{
	  j >>= 1; /* Divide j by 2 */

	  if (threadIdx.y<j)
	    {
	      work_array[threadIdx.x + threadIdx.y*n_row]
		+= work_array[threadIdx.x + threadIdx.y*(n_row+j)];

	      /* Sync threads */
	      __syncthreads();
	    }
	}
      
      /* Copy reduced value to result vector */
      rhs_vector[global_row_id] = work_array[j];
    }
}


/* This kernel used shared memory but assumes:

   - the matrix is stored column wise 
   - the 1D thread block is of size 60
   - 60x60 matrix

*/
__global__ void MultiMatVecMultiply2(int n_mat,
				     int n_row,
				     int n_col,
				     double* lhs_vector,
				     double* matrix,
				     double* rhs_vector)
{
  int local_thread_id;
  int block_id;
  int matrix_id;
  int global_thread_id;
  int row_id;
  int j; 
  double tmp;
  
  /* Get 1D thread index */
  local_thread_id = threadIdx.x;
  
  /* Get 1D block index */
  block_id = blockIdx.x;

  /* Global thread index */
  global_thread_id = local_thread_id + block_id * blockDim.x; 

  /* get matrix id */
  matrix_id = global_thread_id/n_row;

  /* Shared memory to store copy of lhs vector */
  __shared__ double lhs_vector_shared[60];
 /*  double* lhs_vector_shared = (double*)work_array; */
  
  /* Get row id */
  row_id = global_thread_id%n_row; 
  
  /* Load values into shared vector */
  if (local_thread_id < 60)
    lhs_vector_shared[local_thread_id] = 
      lhs_vector[matrix_id*n_col+local_thread_id];

  /* Sunc threads */
  __syncthreads();

  if (global_thread_id < n_mat*n_row)
    {
      tmp = 0.0;
      for (j=0; j<n_col; j++)
	{
	  tmp +=
	    lhs_vector_shared[j] * 
	    matrix[j*n_row+row_id+matrix_id*n_col*n_row];
	}
      rhs_vector[global_thread_id] = tmp;
    }
}

/* 
  This kernel assumes matrix is stored column wise

  This is a naive implementation, no attempt at optimisation
*/
__global__ void MultiMatVecMultiply1(int n_mat,
				     int n_row,
				     int n_col,
				     double* lhs_vector,
				     double* matrix,
				     double* rhs_vector)
{
  int local_thread_id;
  int block_id;
  int matrix_id;
  int global_thread_id;
  int row_id;
  int j; 
  double tmp;
  
  /* Get 1D thread index */
  local_thread_id = threadIdx.x;
  
  /* Get 1D block index */
  block_id = blockIdx.x;

  global_thread_id = local_thread_id + block_id * blockDim.x;  
  
  if (global_thread_id < n_mat*n_row)
    {
      
      /* get matrix id */
      matrix_id = global_thread_id/n_row;
      
      /* Get row id */
      row_id = global_thread_id%n_row; 
      
      tmp = 0.0;
      for (j=0; j<n_col; j++)
	{
	  tmp +=
	    lhs_vector[j+matrix_id*n_col] * 
	    matrix[row_id+j*n_row+matrix_id*n_col*n_row];
	}
      rhs_vector[global_thread_id] = tmp;
    }
}


/* Helper functions */

/* Function to allocate memory on device */
extern "C" int allocate_memory_on_gpu(const int *n_elements, 
				      const int *element_size,
				      void **device_pointer)
{
  cudaError_t cuda_status;
  
  cuda_status = cudaMalloc(device_pointer, *element_size * *n_elements);
  
  if (cuda_status != cudaSuccess)
    {
      printf("Device memory failed to allocate!\n");
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

/* Function to free memory on device */
extern "C" int free_memory_on_gpu(void **device_pointer)
{
  cudaError_t cuda_status;
  
  cuda_status = cudaFree(*device_pointer);
  
  if (cuda_status != cudaSuccess)
    {
      printf("Device memory failed to deallocate!\n");
      return EXIT_FAILURE;
    }
  
  return EXIT_SUCCESS;
}

/* Function to copy data to the gpu */
extern "C" int copy_data_to_gpu(const int *n_elements,
				const int *element_size,
				const void *host_data,
				void **device_pointer)
{
  cudaError_t cuda_status;
  
  cuda_status = cudaMemcpy(*device_pointer,
			   host_data,
			   *n_elements * *element_size,
			   cudaMemcpyHostToDevice);
  
  if (cuda_status != cudaSuccess)
    {
      printf("Failed to copy data to device!\n");
      return EXIT_FAILURE;
    }
  
  return EXIT_SUCCESS;
 }

/* Function to copy data from the gpu */
extern "C" int copy_data_from_gpu(const int *n_elements, 
				  const int *element_size,
				  void *host_data,
				  void **device_pointer)
{
  cudaError_t cuda_status;
  
  cuda_status = cudaMemcpy(host_data,
			   *device_pointer,
			   *n_elements * *element_size,
			   cudaMemcpyDeviceToHost);
  
  if (cuda_status != cudaSuccess)
    {
      printf("Failed to copy data from device!\n");
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

/* Function to call the matrix vector multiply kernel */
extern "C" int matrix_vector_multiplies(int *n_mat,
					int *n_row,
					int *n_col,
					void **d_lhs_vector,
					void **d_matrix,
					void **d_rhs_vector)
{
  /* Syntax <<<NumBlocks, ThreadsPerBlock>>> */
  /* Max no threads per dimension per block 1024 */
  /* Max size of grid in each dimension 65535 */
  /* Warp size 32 */
  int NumBlocks;
  int ThreadsPerBlock_1D;
  dim3 ThreadsPerBlock_2D;
  cudaError_t cuda_status;

  int method = 2;

  size_t shared_mem_size;

  if (method == 3)
    {
      ThreadsPerBlock_1D = 1024;
      NumBlocks = (*n_mat * *n_row)/ThreadsPerBlock_1D;
      
      if ( (*n_mat * *n_row)%ThreadsPerBlock_1D != 0)
	NumBlocks += 1;
      
      /* Launch kernel */
      MultiMatVecMultiply1<<<NumBlocks, ThreadsPerBlock_1D>>>
	(*n_mat,
	 *n_row,
	 *n_col,
	 (double *)(*d_lhs_vector),
	 (double *)(*d_matrix),
	 (double *)(*d_rhs_vector));
    }
  else if (method == 2)
    {
      ThreadsPerBlock_1D = *n_row;
      NumBlocks = *n_mat; /* (*n_mat * *n_row)/ThreadsPerBlock_1D;
      
      if ( (*n_mat * *n_row)%ThreadsPerBlock_1D != 0)
      NumBlocks += 1; */
      
      /* Launch kernel */
      MultiMatVecMultiply2<<<NumBlocks, ThreadsPerBlock_1D>>>
	(*n_mat,
	 *n_row,
	 *n_col,
	 (double *)(*d_lhs_vector),
	 (double *)(*d_matrix),
	 (double *)(*d_rhs_vector));
      
    }
  else
    {
      NumBlocks = *n_mat;
      
      ThreadsPerBlock_2D.x = *n_row;
      ThreadsPerBlock_2D.y = 2;
      ThreadsPerBlock_2D.y = 1;
      
      shared_mem_size = sizeof(*d_matrix) * 
	ThreadsPerBlock_2D.x * ThreadsPerBlock_2D.y;

      MultiMatVecMultiply3
	<<<NumBlocks, ThreadsPerBlock_2D, shared_mem_size>>>
	(*n_mat,
	 *n_row,
	 *n_col,
	 (double *)(*d_lhs_vector),
	 (double *)(*d_matrix),
	 (double *)(*d_rhs_vector));
    }


  /* Test for error at kernel launch */
  cuda_status = cudaGetLastError();
  if (cuda_status != cudaSuccess)
    {
      printf("Kernel launch failure!\n");
      return EXIT_FAILURE;
    }
  
  /* Ensure synchronisation */
  /* Note: this should pick up kernel errors */
  cuda_status = cudaDeviceSynchronize();
  if (cuda_status != cudaSuccess)
    {
      printf("Streams failed to synchronise!\n");
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
 
 
