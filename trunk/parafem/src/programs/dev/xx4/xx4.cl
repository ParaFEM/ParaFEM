#ifdef cl_khr_fp64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#else
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

// Naive Matrix vector multiplication kernel
__kernel void  MultiMatVecMultiply_col_naive( const int n_mat, const int n_row, const int n_col,
					      __global const double *matrix, 
					      __global const double *x_vector, 
					      __global double *y_vector )
{
  double accum;
  int matrix_id;
  int row_id;
  int global_thread_id = get_global_id(0);
  const int mat_size = n_col*n_row;

  // Each work-item (thread) calculates a single element of the
  // Y result vector by doing the dot product of an entire row
  // from a matrix with an entire X vector. Hence we have
  // n_mat*n_row work-items (threads). 

  if (global_thread_id < n_mat*n_row) {

    // Which matrix are we in
    matrix_id = global_thread_id / n_row;

    // Which row are we in within that matrix
    row_id = global_thread_id % n_row;

    // Do the A.x for the current matrix's row and X vector
    accum = 0.0;
    for ( int col=0; col<n_col; col++ ) {

      // Stride to current matrix: matrix_id*mat_size
      // Stride to column 'col' in current row of matrix (using column-major order): col*n_row+row_id
      // Stride to current X vector's row for the current col: matrix_id*n_col + col
      accum += matrix[(matrix_id*mat_size) + col*n_row+row_id  ] * x_vector[matrix_id*n_col + col];
    }

    // There's a work-item (thread) per result vector (Y) row
    y_vector[global_thread_id] = accum;
  }
};

// Matrix vector multiplication kernel using local memory for vector re-use
__kernel void  __attribute__((reqd_work_group_size(60,1,1))) MultiMatVecMultiply_col_local( const int n_mat, const int n_row, const int n_col,
											    __global const double *matrix, 
											    __global const double *x_vector, 
											    __global double *y_vector )
{
  __local double x_vector_local[60];
  double accum;
  int matrix_id;
  int matrix_offset;
  int row_id;
  int global_thread_id = get_global_id(0);
  int local_thread_id  = get_local_id(0);   // work-item's ID in this work-group
  int local_group_size = get_local_size(0); // size of this work-group (should be a 60)

  // Each work-item (thread) calculates a single element of the
  // Y result vector by doing the dot product of an entire row
  // from a matrix with an entire X vector. Hence we have
  // n_mat*n_row work-items (threads). 
  //
  // We now have a 1D work-group of size n_rows. Each work-item
  // (thread) in the work-group will operate on a row from the
  // same matrix. They will all use the same X vector. Hence the
  // particular X vector can be copied in to local memory to be
  // shared amongst the threads in the work-group.

  // Which matrix are we in (should now be same at get_group_id(0))
  matrix_id = global_thread_id / n_row;

  // Load X vector values shared vector
  if (local_thread_id < 60)
    x_vector_local[local_thread_id] = x_vector[matrix_id*n_col+local_thread_id];

  // Sunc threads
  barrier(CLK_LOCAL_MEM_FENCE);

  if (global_thread_id < n_mat*n_row) {

    // Which row are we in within our matrix
    row_id = global_thread_id % n_row;
    matrix_offset = matrix_id * n_col * n_row;

    // Do the A.x for the current matrix's row and X vector
    accum = 0.0;
    for ( int col=0; col<n_col; col++ ) {

      // Stride to current matrix: matrix_id*matrix_size
      // Stride to column 'col' in current row of matrix (using column-major order): col*n_row+row_id
      // X vector is now local so just use the current column to access it
      accum += matrix[matrix_offset + col*n_row+row_id] * x_vector_local[col];
    }

    // There's a work-item (thread) per result vector (Y) row
    y_vector[global_thread_id] = accum;
  }
};
