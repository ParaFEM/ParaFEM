#ifdef cl_khr_fp64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#else
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

// Naive Matrix vector multiplication kernel using 2D NDRange
// Global work size should be set to [N,M] where result matrix
// C is an MxN matrix. Assumes fortran column-major array layout.
__kernel void MatMatMultiply_col_2d( const int n_rowsA,
				     const int n_colsArowsB,
				     const int n_colsB,
				     __global const double *matrixA, 
				     __global const double *matrixB, 
				     __global double *matrixC )
{
  double accum;
  int col = get_global_id(0);
  int row = get_global_id(1);

  if ( col < n_colsB && row < n_rowsA ) {

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
  }

};


// Assumes fortran column-major array layout. Will exhibit poor coalesced memory reads.
__kernel void  __attribute__((reqd_work_group_size(60,1,1))) MatMatMultiply_col_1d( const int n_rowsA, 
										    const int n_colsArowsB,
										    const int n_colsB,
										    __global const double *matrixA, 
										    __global const double *matrixB, 
										    __global double *matrixC )
{
  __local double a_row_local[60];	    // Hard-coded for paraFEM's 60x60 matrix
  double accum;

  int global_thread_id = get_global_id(0);  // work-item's global thread ID
  int local_thread_id  = get_local_id(0);   // work-item's thread ID in this work-group

  // Offset in to matrix A to get to start of a column
  int a_col_offset = local_thread_id*n_rowsA;

  // Offset in to matrix B to get to start of a column
  int b_col_offset = global_thread_id*n_colsArowsB;

  // Offset in to matrix C to get to start of a column
  int c_col_offset = global_thread_id*n_rowsA;

  int local_group_size = get_local_size(0); // size of this work-group (should be 60)

  // There is a global work-item (thread) for each column in a
  // particular row of the C result (and hence also B) matrix.
  // The number of rows is current fixed (in the xx10 code) to
  // be 60. The number of columns increases with dataset size.
  // Each work-item (thread) calculates the elements in an entire
  // column of the C result matrix. Each element is calculated by
  // doing a dot product of an A matrix row with a B matrix column.
  // The work-items in a work group (of size 60) operate on the same
  // row of matrix A. Hence we copy that row in to local memory.

  if (global_thread_id < n_colsB) {

    // For each A matrix row
    for ( int ac_row=0; ac_row<n_rowsA; ac_row++ ) {

      // Load matrix A's row values in to local memory
      if (local_thread_id < 60)
	a_row_local[local_thread_id] = matrixA[a_col_offset+ac_row];
      
      // Sync threads
      barrier(CLK_LOCAL_MEM_FENCE);

      // Do the Arow.Bcol for the current A matrix's row and B matrix column
      accum = 0.0;

      for ( int b_row=0; b_row<n_colsArowsB; b_row++ ) {
	accum += a_row_local[b_row] * matrixB[b_col_offset+b_row];
      }

      // Result in matrix C
      matrixC[c_col_offset+ac_row] = accum;

      // Sync threads
      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }
}
