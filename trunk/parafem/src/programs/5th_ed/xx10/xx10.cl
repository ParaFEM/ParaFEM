#ifdef cl_khr_fp64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#else
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

// Naive Matrix vector multiplication kernel using 2D NDRange
// Global work size should be set to [N,M] where result matrix
// C is an MxN matrix.
__kernel void  MultiMatMatMultiply_col( const int n_colsArowsB,
					__global const double *matrixA, 
					__global const double *matrixB, 
					__global double *matrixC )
{
  double accum;
  int col = get_global_id(0);
  int row = get_global_id(1);
  int n_colsB = get_global_size(0);
  int n_rowsA = get_global_size(1);

  // Each work-item (thread) calculates a single element of the
  // C result matrix by doing the dot product of an entire row
  // from the A matrix with an entire column from the B matrix.
  // Hence we have n_rowsA * n_colsB work items. Use a 2D NDRange.

  accum = 0.0;
  // i gives a column in A and a row in B. Stride in to the matrices
  // assuming column major order (from fortran)
  for ( int i=0; i<n_colsArowsB; i++ ) {
    accum += matrixA[row + n_rowsA*i] * matrixB[i + n_colsArowsB*col];
  }
  matrixC[col*n_rowsA + row] = accum;

};

