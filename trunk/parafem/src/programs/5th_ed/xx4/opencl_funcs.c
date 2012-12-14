/* ----------------------------------------------------------------------
 * opencl_funcs.c
 * ----------------------------------------------------------------------
 * C wrappers for various OpenCL functions that can be called from
 * fortran.
 *
 * Most of the args are pointers because fortran use call-by-reference
 * for function args.
 *
 * Note, in fortran we store a cl_mem object as a c_ptr from the ISO C
 * Binding. Hence when you pass in a fortran var of type c_ptr to a C
 * function the C function should see this as a void** type.
 *
 * To access the cl_mem object in the C code:
 *
 * func_(void **device_ptr )
 * {
 *   // Either of the following work (which is 'correct'?)
 *   cl_mem memobj = *(cl_mem *)device_ptr;
 *   cl_mem memobj = (cl_mem)((cl_mem*)*device_ptr);
 *
 *  // When assigning a cl_mem obj to an 'out arg' pointer we do
 *  *device_ptr = (void *)memobj;
 * }
 *   
 * ----------------------------------------------------------------------
 */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sys/types.h>
#include<math.h>

#include <CL/cl.h>

// See our makefile. Enable AMD clBLAS functionality.
#ifdef AMDBLAS
#include<clAmdBlas.h>
#endif

// ---------------------------------------------------------------------
// Report function names
//#define FNAMES

#ifdef FNAMES
#define SNAME( _str ) { printf( "Start %s()...\n", _str ); fflush(stdout); }
#define ENAME( _str ) { printf( "End %s()\n", _str ); fflush(stdout); }
#else
#define SNAME( _str )
#define ENAME( _str )
#endif

#ifdef _MIN
#undef _MIN
#endif
#define _MIN(_a,_b) ((_a)<(_b)?(_a):(_b))

// ---------------------------------------------------------------------

/* ----------------------------------------------------------------------
 * OpenCL context info required by various cl functions
 * ----------------------------------------------------------------------
 */
typedef struct {
  cl_uint		platform_count;
  cl_platform_id	*platforms;	// All platforms
  cl_platform_id	platform;	// The platform we are using
  cl_uint		device_count;
  cl_device_id		*devices;	// All the devices
  cl_device_id		device;		// The device we are using
  cl_context_properties props[3];	// {CL_CONTEXT_PLATFORM, 0, 0};
  cl_context		ctx;
  cl_command_queue	queue;
  cl_program		program;
  cl_kernel		kernel;
  int			clblas;
  cl_char		device_name[512];
  cl_uint		device_vendor_id;
  cl_long		device_max_mem_alloc;
  cl_long		device_global_mem;
} OCLinfo;

/* ----------------------------------------------------------------------
 * Private globals
 * ----------------------------------------------------------------------
 */

// Access context etc between function calls (not thread safe)
static OCLinfo *g_oclinfo = NULL;

// Names of OCL mem buffer types for error messages
static char *memtype_names[] = {"RO", "WO", "RW"};

/* ----------------------------------------------------------------------
 * Functions to be called from fortran (hence all args are pointers)
 * ----------------------------------------------------------------------
 */

/* ----------------------------------------------------------------------
 * init_opencl_()
 * ----------------------------------------------------------------------
 * Initialize a device for running OpenCL kernels. This does the usual
 * things - selects a device, creates a context, etc. Only one device
 * is used and currently only the type (GPU or CPU can be specified).
 * A specific device number (0,...,num_devices-1) can be requested or
 * we can loop over all devices and use the first one that allows a
 * context to be successfully created.
 *
 * Note that if the GPU is set to be in exclusive mode (by sys admin)
 * then you will only find this out when you try to create a context.
 * It will fail if the device is already in use. Hence the code should
 * get all devices and loop over them if you can't guarantee you have
 * all of the GPUs to yourself.
 *
 * If required the AMD clBlas library is initialized.
 *
 * The context, command queue etc is stored locally for other functions
 * to use. This isn't really a general toolkit - only one context and
 * command queue can currently be used. You might need more.
 *
 * This should be the first 'opencl' function called by the host code.
 *
 * Args:
 * In
 *	int *use_gpu:	1=GPU device, 0=CPU device
 *	int *use_device_id: 0,1,2... Which device to use if more than one
 *			exists. If -1 then use first free device.
 *	int *use_blas:	1=Init AMD clBlas library, 0=don't!
 *
 * Returns 1 on success, 0 on error.
 * ----------------------------------------------------------------------
 */
int init_opencl_(const int *use_gpu, const int *use_device_id, const int *use_blas)
{
  OCLinfo	*oclinfo;
  cl_int	err;
  cl_int	i;

  SNAME( "init_opencl" );

  // Private struct to hold OCL items between functions
  oclinfo = (OCLinfo *)calloc(1, sizeof(OCLinfo));
  if ( !oclinfo )
    return 0;

  // Get number of platforms. Usually only one but could have multiple vendors.
  err = clGetPlatformIDs(0, NULL, &oclinfo->platform_count);
  if ( oclinfo->platform_count == 0 ) {
    printf( "Error: No OpenCL platforms found. Is OpenCL installed?\n" );
    free(oclinfo);
    return 0;
  }
  // Get all the platform IDs.
  oclinfo->platforms = (cl_platform_id *)malloc(oclinfo->platform_count*sizeof(cl_platform_id));
  err = clGetPlatformIDs(oclinfo->platform_count, oclinfo->platforms, NULL);
  if ( err != CL_SUCCESS ) {
    printf( "Failed to get any OpenCL platforms\n" );
    free(oclinfo);
    return 0;
  }

  // Here we should use clGetPlatformInfo() to check things like the
  // CL_VENDOR_STRING to look for a specific platform. For now we
  // assume only one platform and use the first given to us.
  oclinfo->platform = oclinfo->platforms[0];
  printf( "Using platform id 0 of %d\n", oclinfo->platform_count );

  // There may be multiple device. Find out how many.
  err = clGetDeviceIDs(oclinfo->platform, (*use_gpu?CL_DEVICE_TYPE_GPU:CL_DEVICE_TYPE_CPU), 0, 
		       NULL, &oclinfo->device_count);
  
  if ( err != CL_SUCCESS || oclinfo->device_count == 0 ) {
    printf( "Error: No OpenCL devices reported\n" );
    return 0;
  }
    
  // Now get all the device_ids
  oclinfo->devices = (cl_device_id *)malloc( oclinfo->device_count * sizeof(cl_device_id) );
  err = clGetDeviceIDs(oclinfo->platform, (*use_gpu?CL_DEVICE_TYPE_GPU:CL_DEVICE_TYPE_CPU), 
		       oclinfo->device_count, oclinfo->devices, NULL);
  if ( err != CL_SUCCESS ) {
    printf( "Failed to get any OpenCL %s devices\n", (*use_gpu?"GPU":"CPU") );
    free(oclinfo->platforms);
    free(oclinfo->devices);
    free(oclinfo);
    return 0;
  }

  // We'll deal with the devices below. First set up the context properties
  // since we'll try to create a context as part of the device availability
  // check.
  oclinfo->props[0] = CL_CONTEXT_PLATFORM;
  oclinfo->props[1] = (cl_context_properties)oclinfo->platform;
  oclinfo->props[2] = 0;

  // Now either use a requested device id (0,1,...) or search for the first free device.
  if ( *use_device_id > -1 ) {

    // Use the requested device
    if ( *use_device_id >= oclinfo->device_count ) {
      printf( "Error: Only %d devices present. Requested device id %d out of range.\n",
	      oclinfo->device_count, *use_device_id );
      free(oclinfo->platforms);
      free(oclinfo->devices);
      free(oclinfo);
      return 0;
    }

    // Note, there is no guarantee the user has requested a device that is not
    // already in use. Only when we try to create the context will we detect this.
    // Getting the CL_DEVICE_AVAILABLE flag doesn't seem to work.
    oclinfo->device = oclinfo->devices[*use_device_id];
    printf( "Using requested device id %d\n", *use_device_id  );

    // Create the context. Error check done below.
    oclinfo->ctx = clCreateContext(oclinfo->props, 1, &oclinfo->device, NULL, NULL, &err);
  }
  else {

    // Loop over all devices. We'll use the first one that allows a context to be created.
    for ( i=0; i<oclinfo->device_count; i++ ) {
      oclinfo->device = oclinfo->devices[i];
      oclinfo->ctx = clCreateContext(oclinfo->props, 1, &oclinfo->device, NULL, NULL, &err);
      if ( oclinfo->ctx != 0 && err == CL_SUCCESS ) {
	printf( "Using device %d\n", i );
	break;
      }
    }
  }
  if ( err != CL_SUCCESS ) {
    printf( "Failed to get an OpenCL context on any device\n" );
    free(oclinfo->platforms);
    if (oclinfo->devices)
      free(oclinfo->devices);
    free(oclinfo);
    return 0;
  }

  // Create command queue in the context
  oclinfo->queue = clCreateCommandQueue(oclinfo->ctx, oclinfo->device, 0, &err);
  if ( err != CL_SUCCESS ) {
    printf( "Failed to create an OpenCL command queue\n" );
    clReleaseContext(oclinfo->ctx);
    if (oclinfo->devices)
      free(oclinfo->devices);
    free(oclinfo);
    return 0;
  }

  // Get device name string
  err = clGetDeviceInfo(oclinfo->device, CL_DEVICE_NAME, 
			sizeof(oclinfo->device_name), &oclinfo->device_name, NULL );
  if ( err != CL_SUCCESS || oclinfo->device_name[0] == '\0' ) {
    printf( "Failed to get CL_DEVICE_NAME\n" );
  }
  // Query device for mem alloc sizes (add other properties as required).
  // We might want the max memory alloc and global size for error checking.
  err = clGetDeviceInfo(oclinfo->device, CL_DEVICE_MAX_MEM_ALLOC_SIZE, 
			sizeof(oclinfo->device_max_mem_alloc), 
			&oclinfo->device_max_mem_alloc, NULL );
  if ( err != CL_SUCCESS || oclinfo->device_max_mem_alloc == 0 ) {
    printf( "Failed to get CL_DEVICE_MAX_MEM_ALLOC_SIZE\n" );
  }
  err = clGetDeviceInfo(oclinfo->device, CL_DEVICE_GLOBAL_MEM_SIZE, 
			sizeof(oclinfo->device_global_mem), 
			&oclinfo->device_global_mem, NULL );
  if ( err != CL_SUCCESS || oclinfo->device_global_mem == 0 ) {
    printf( "Failed to get CL_DEVICE_GLOBAL_MEM_SIZE\n" );
  }
  err = clGetDeviceInfo(oclinfo->device, CL_DEVICE_VENDOR_ID, 
			sizeof(oclinfo->device_vendor_id), 
			&oclinfo->device_vendor_id, NULL );
  if ( err != CL_SUCCESS ) {
    printf( "Failed to get CL_DEVICE_GLOBAL_MEM_SIZE\n" );
  }


  printf( "Device name:          %s\nDevice vendor id:     %d\nDevice global mem:    %ld\nDevice max mem alloc: %ld\n", 
	  oclinfo->device_name,
	  oclinfo->device_vendor_id,
	  oclinfo->device_global_mem, 
	  oclinfo->device_max_mem_alloc ); 
  fflush(stdout);

  // Set up AMD Blas library
#ifdef AMDBLAS
  if ( *use_blas ) {
    err = clAmdBlasSetup();
    if (err != CL_SUCCESS) {
      printf("Failed to set up in clAmdBlasSetup(): %d\n", err);
      clReleaseCommandQueue(oclinfo->queue);
      clReleaseContext(oclinfo->ctx);
      free(oclinfo->platforms);
      if (oclinfo->devices)
	free(oclinfo->devices);
      free(oclinfo);
      return 0;
    }
    oclinfo->clblas = 1;
  }
  else {
    oclinfo->clblas = 0;
  }
#else
  oclinfo->clblas = 0;
#endif

  // Save in global ptr
  g_oclinfo = oclinfo;

  ENAME( "init_opencl" );
  return 1;
}

/* ----------------------------------------------------------------------
 * stop_opencl_()
 * ----------------------------------------------------------------------
 * Clean up the OpenCL objects and shutdown the AMD clBlas library if
 * it has been initialized.
 *
 * Someone else should release any device mem objects before calling this.
 * This should really be the last 'opencl' function you call.
 *
 * Returns 1 on success, 0 on error.
 * ----------------------------------------------------------------------
 */
int stop_opencl_()
{
  SNAME( "stop_opencl" );

#ifdef AMDBLAS
  if ( g_oclinfo->clblas )
    clAmdBlasTeardown();
#endif
  
  clReleaseProgram(g_oclinfo->program);
  clReleaseKernel(g_oclinfo->kernel);
  clReleaseCommandQueue(g_oclinfo->queue);
  clReleaseContext(g_oclinfo->ctx);
  free(g_oclinfo->platforms);
  if (g_oclinfo->devices)
    free(g_oclinfo->devices);
  free(g_oclinfo);
  g_oclinfo = NULL;

  ENAME( "stop_opencl" );
  return 1;
}

/* ----------------------------------------------------------------------
 * compile_kernel_from_file_()
 * ----------------------------------------------------------------------
 * Read OpenCL kernel source from a text file and compile the code.
 * The source code can contain more than one kernel but only the named
 * kernel will be compiled in to the current program/context.
 *
 * The oclinfo struct would need to be modified to support multiple
 * kernels. Currently we only use one kernel at a time.
 *
 * The build options string is passed through to clBuildProgram().
 *
 * If an error occurs during compilation the build log from the compiler
 * is printed to stdout.
 *
 * Returns 1 on success, 0 on error.
 * ----------------------------------------------------------------------
 */
int compile_kernel_from_file_( const char *filename, const char *kernelname, const char *build_options )
{
  cl_int err;
  FILE	 *kfp;
  size_t kfilesize;
  char	 *ksource;
  size_t numread;
  size_t log_size;
  char   *build_log;
  
  SNAME( "compile_kernel_from_file" );

  // Load kernel source from file
  kfp = fopen(filename, "r");
  if( !kfp ) {
    printf("Could not open kernel source file %s\n", filename );
    return 0;
  }

  // Find out how big the file is then read it
  fseek(kfp , 0 , SEEK_END);
  kfilesize = (long)ftell(kfp);
  rewind(kfp);
  ksource = (char *)malloc((kfilesize+1)*sizeof(char));
  numread = fread(ksource, 1, kfilesize, kfp);
  ksource[kfilesize] = '\0';
  if ( numread < kfilesize ) {
    printf( "Error reading kernel source file %s. Only read %ld of %ld bytes\n", 
	    filename, numread, kfilesize );
    free(ksource);
    return 0;
  }
  fclose(kfp);
  // printf( "%s\n", ksource ); fflush(stdout);

  // Create the program object for the context from the source text
  g_oclinfo->program = clCreateProgramWithSource(g_oclinfo->ctx, 1, (const char **) &ksource, NULL, &err);
  if ( err != CL_SUCCESS ) {
    printf( "Failed to create program (err %d)\n", err );
    free(ksource);
    return 0;
  }
  
  // Have finished with the source text
  free(ksource);
  
  // Build (compile and link) the program executable
  err = clBuildProgram(g_oclinfo->program, 0, NULL, build_options, NULL, NULL);
  
  // Get the build log (warnings, errors, ... )
  clGetProgramBuildInfo(g_oclinfo->program, g_oclinfo->device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
  build_log = (char *)malloc(log_size);
  clGetProgramBuildInfo(g_oclinfo->program, g_oclinfo->device, CL_PROGRAM_BUILD_LOG, log_size, build_log, NULL);
  
  if ( err != CL_SUCCESS ) {    
    printf("Failed to build program executable (err %d)\n", err);
    printf("%s\n", build_log);
    free(build_log);
    return 0;
  }
  printf( "Build Log:\n%s\n", build_log );
  free(build_log);

  // Create the kernel. The caller must supply the kernel name defined in the source file.
  g_oclinfo->kernel = clCreateKernel(g_oclinfo->program, kernelname, &err);
  if (!g_oclinfo->kernel || err != CL_SUCCESS)
  {
    printf("Failed to create compute kernel %s (err %d)\n", kernelname, err);
    return 0;
  }

  ENAME( "compile_kernel_from_file" );
  return 1;
}


/* ----------------------------------------------------------------------
 * allocate_memory_on_gpu_()
 * ----------------------------------------------------------------------
 * Allocate a device memory object of the requested size and type.
 * OpenCL supports several read/write modes for the mem obj but currently
 * we are only interested on CL_MEM_READ_ONLY, CL_MEM_WRITE_ONLY and
 * CL_MEM_READ_WRITE. You might want to experiment with host mapped
 * arrays. We pass in a 0,1,2 flag to represent the above types.
 *
 * OpenCL allows a host array to be copied in to the mem obj at creation
 * time. Again, we don't currently do that but instead split this in to
 * creation and population functions. See below for reading from and
 * writing to a device mem object.
 *
 * Args:
 * In
 *	int *n_elems:	   Number of elements required in mem obj
 *	int *elem_size:    Size in bytes of one element
 *	int memflag:	   Our flag: 0=CL_MEM_READ_ONLY
 *				     1=CL_MEM_WRITE_ONLY
 *				     2=CL_MEM_READ_WRITE
 * Out
 *	void **device_ptr: Ptr to device mem obj
 *
 * Returns 1 on success, 0 on error.
 * ----------------------------------------------------------------------
 */
int allocate_memory_on_gpu_( const int *n_elems, const int *elem_size, const int *memflag, void **device_ptr )
{
  // Assume OCLinfo exists
  cl_int err;
  cl_mem memobj;
  cl_mem_flags memtype;


  SNAME( "allocate_memory_on_gpu" );

  switch(*memflag) {
  case 0:
    memtype = CL_MEM_READ_ONLY;
    break;
  case 1:
    memtype = CL_MEM_WRITE_ONLY;
    break;
  case 2:
    memtype = CL_MEM_READ_WRITE;
    break;
  default:
    printf( "Invalid memflag arg %d\n", *memflag );
    return 0;
  }
  memobj = clCreateBuffer( g_oclinfo->ctx, memtype, (*n_elems * *elem_size), NULL, &err );

  if ( err != CL_SUCCESS ) {
    printf( "Failed to allocate %d x %d bytes %s device memory (err %d)\n", 
	    *n_elems, *elem_size, memtype_names[*memflag], err );
    return 0;
  }

  // printf( "Allocated %d x %d bytes on GPU\n", *n_elems, *elem_size ); fflush(stdout);

  // We store a cl_mem as a fortran c_ptr which corresponds to a void* in C.
  *device_ptr = (void *)memobj;

  ENAME( "allocate_memory_on_gpu" );
  return 1;
}

/* ----------------------------------------------------------------------
 * free_memory_on_gpu_()
 * ----------------------------------------------------------------------
 * Release a device mem obj
 *
 * Args:
 * In
 *	void **device_ptr: Ptr to device mem obj
 *
 * Returns 1 on success, 0 on error
 * ----------------------------------------------------------------------
 */
int free_memory_on_gpu_( void **device_ptr )
{
  cl_int err;
  cl_mem memobj;

  SNAME( "free_memory_on_gpu" );

  memobj = *(cl_mem *)device_ptr;

  err = clReleaseMemObject(memobj);
  if( err != CL_SUCCESS ) {
    printf( "Failed to release device memory\n" );
    return 0;
  }

  ENAME( "free_memory_on_gpu" );
  return 1;
}

/* ----------------------------------------------------------------------
 * copy_data_to_gpu_()
 * ----------------------------------------------------------------------
 * Do a blocking write of memory from the host array to a device mem obj.
 * The host array and device mem obj must already have been allocated.
 *
 * Args:
 * In
 *	int *n_elems:	   Number of elements to copy
 *	int *elem_size:    Size in bytes of one element
 *	void *host_data:   Source array on the host
 * Out
 *	void **device_ptr: Ptr to device mem obj
 *
 * Returns 1 on success, 0 on error.
 * ----------------------------------------------------------------------
 */
int copy_data_to_gpu_( const int *n_elems, const int *elem_size, const void *host_data, void **device_ptr )
{
  cl_int err;
  cl_mem memobj;

  SNAME( "copy_data_to_gpu" );

  memobj = *(cl_mem *)device_ptr;

  err = clEnqueueWriteBuffer(g_oclinfo->queue, memobj, CL_TRUE, 0, (*n_elems * *elem_size), host_data,
			     0, NULL, NULL );

  if ( err != CL_SUCCESS ) {
    printf( "Failed to copy %d x %d bytes to GPU (err %d)\n", *n_elems, *elem_size, err );
    return 0;
  }

  // printf( "Copied %d x %d bytes to GPU\n", *n_elems, *elem_size ); fflush(stdout);

  ENAME( "copy_data_to_gpu" );
  return 1;
}

/* ----------------------------------------------------------------------
 * copy_data_from_gpu_()
 * ----------------------------------------------------------------------
 * Do a blocking read of memory from a device mem obj to a host array.
 * The host array and device mem obj must already have been allocated.
 *
 * Args:
 * In
 *	int *n_elems:	   Number of elements to copy
 *	int *elem_size:    Size in bytes of one element
 *	void **device_ptr: Ptr to device mem obj
 * Out
 *	void *host_data:   Resulting array on the host
 *
 * Returns 1 on success, 0 on error.
 * ----------------------------------------------------------------------
 */
int copy_data_from_gpu_( const int *n_elems, const int *elem_size, const void **device_ptr, void *host_data )
{
  cl_int err;
  cl_mem memobj;

  SNAME( "copy_data_from_gpu" );

  memobj = *(cl_mem *)device_ptr;

  err = clEnqueueReadBuffer(g_oclinfo->queue, memobj, CL_TRUE, 0, (*n_elems * *elem_size), host_data,
			    0, NULL, NULL );

  if ( err != CL_SUCCESS ) {
    printf( "Failed to copy %d x %d bytes from GPU (err %d)\n", *n_elems, *elem_size, err );
    return 0;
  }

  ENAME( "copy_data_from_gpu" );
  return 1;
}

/* ----------------------------------------------------------------------
 * calc_array_blocks_for_gpu_() - currently unused
 * ----------------------------------------------------------------------
 * Specific to ParaFEM p121. Given the number of elements in the arrays
 * containing all the matrices and vectors to be multiplied, see if the
 * entire arrays can be loaded on to the GPU. If not, work out how many
 * blocks the arrays must be split in to. We must also allow room for
 * the result vector (the result of multiplying a matrix by a vector).
 * All three items need space on the GPU at once.
 *
 * In the code the X, Y vectors refer to the RHS and LHS vectors:
 *
 * vec_y = Mat . vec_x
 *
 * Some of our test cases (those where we write our own kernels) really
 * want the entire arrays to be on the GPU, rather than transferring
 * individual matrices and vectors to the GPU for multiplication (we do
 * that for one of the AMD clBlas tests). But for larger test cases it
 * isn't possible to copy the entire arrays. So we'll do the best we
 * can and split the arrays in to large blocks. Here we calculate how
 * many blocks and the offset to each block. Note this has nothing to
 * do with work group sizes in OpenCL (NDRange etc). 
 *
 * We split the arrays in to blocks of complete matrices and vectors
 * (in p121 each matrix is 60x60 and a vector is 60x1). They are stored
 * in fortran column-major order so for an OpenCL kernel to have a 
 * complete matrix row to work with we must pass a complete matrix to
 * the kernel.
 *
 * Args: See below
 *
 * Returns 1 on success, 0 on error.
 * ----------------------------------------------------------------------
 */
int calc_array_blocks_for_gpu_( const int *n_matsA,	// Total number of matrices (and X,Y vectors)
				const int *n_rowsA,	// Num rows in one matrix (and rows of one Y vector)
				const int *n_colsA,	// Num cols in one matrix (and rows of one X vector)
				const int *elem_size,	// Size of an element (in bytes) in an array and vector
				int *num_blocks,	// Number of blocks in to which arrays are divided
				int *num_matvecxy_per_block,// Number of [matrix+x_vec+y_vec] per block
				int *mat_block_elems,	// Size (elems) of a matrix block (last one may be smaller)
				int *vec_y_block_elems, // Size (elems) of a  Y vector block (last one may be smaller)
				int *vec_x_block_elems,  // Size (elems) of an X vector block (last one may be smaller)
				int *num_matvecxy_in_last_block  // Last/remainder mvblock may be smaller
  )
{
  // Use long for byte arithmetic
  cl_int err;
  long device_max_mem_alloc_orig = 0;
  long device_max_mem_alloc = 0;
  long device_global_mem = 0;
  long mat_bytes;
  long vec_x_bytes;
  long vec_y_bytes;
  long matvec_xy_bytes;
  long matvec_xy_bytes_all;
  int  max_num_matvecxy_per_block;
  int  max_mem_objs;

  SNAME( "calc_array_blocks_for_gpu" );

  // *** FUDGE ***
  // NVidia hardware lets you allocate more than CL_DEVICE_MAX_MEM_ALLOC_SIZE.
  // Hence we just allow the code to use the entire global memory.
  device_max_mem_alloc_orig = g_oclinfo->device_max_mem_alloc;
  device_max_mem_alloc      = g_oclinfo->device_global_mem;
  device_global_mem         = g_oclinfo->device_global_mem;

  // Number of memory buffer objects we could allocate on the device.
  // For now we just use one.
  max_mem_objs = device_global_mem / device_max_mem_alloc;

  // Size in bytes of a matrix, an X vector and a Y vector
  mat_bytes       = (*n_rowsA) * (*n_colsA) * (*elem_size);
  vec_x_bytes     = (*n_colsA) * (*elem_size);
  vec_y_bytes     = (*n_rowsA) * (*elem_size);

  // Size in bytes needed to store a matrix + a Y vector + an X vector
  matvec_xy_bytes     = mat_bytes + vec_x_bytes + vec_y_bytes;
  matvec_xy_bytes_all = (*n_matsA) * matvec_xy_bytes;

  // Return values
  // -------------

  // Max number of [mat,Xvec,Yvec] triples that can fit in to the device mem
  max_num_matvecxy_per_block = (int)(device_max_mem_alloc / matvec_xy_bytes);
  if ( max_num_matvecxy_per_block == 0 ) {
    printf( "Error: A single matrix+Xvec+Yvec cannot fit on the GPU!\n" );
    return 0;
  }

  // Number of blocks we need to store all mat,Xvec,Yvec data.
  // The caller will need to do this many iterations to
  // process all mat,vecs. 
  *num_blocks = ((*n_matsA) / max_num_matvecxy_per_block) + (((*n_matsA) % max_num_matvecxy_per_block) == 0?0:1);
  *num_matvecxy_per_block = _MIN( *n_matsA, max_num_matvecxy_per_block );
  *num_matvecxy_in_last_block = (*n_matsA) % max_num_matvecxy_per_block;

  // Size of blocks in elements (not bytes)
  *mat_block_elems   = *num_matvecxy_per_block * (*n_rowsA) * (*n_colsA);
  *vec_x_block_elems = *num_matvecxy_per_block * (*n_colsA);
  *vec_y_block_elems = *num_matvecxy_per_block * (*n_rowsA);

#if 1
  printf( \
"Total num mat,vecx,y entries:             %d\n" \
"Num bytes per single mat,vecx,y entry:    %ld\n" \
"Num bytes for all matrices and vectors:   %ld\n" \
"Device MAX_MEM_ALLOC bytes:               %ld\n" \
"Device GLOBAL_MEM bytes:                  %ld\n" \
"Num max-sized mem objects:                %d\n" \
"Max Num mat,vecx,y entries per mem-block: %d\n" \
"Num mem-blocks needed:                    %d\n" \
"Num mat,vecx,y entries per mem-block:     %d\n" \
"Num mat,vecx,y entries in last mem-block: %d\n" \
"Num elems per block for matrices:         %d\n" \
"Num elems per block for x vectors:        %d\n" \
"Num elems per block for y vectors:        %d\n",
*n_matsA, matvec_xy_bytes, matvec_xy_bytes_all, device_max_mem_alloc_orig, device_global_mem, max_mem_objs,
max_num_matvecxy_per_block, *num_blocks, *num_matvecxy_per_block, *num_matvecxy_in_last_block, 
*mat_block_elems, *vec_x_block_elems, *vec_y_block_elems ); fflush(stdout);
#endif

  return 1;
  ENAME( "calc_array_blocks_for_gpu" );  
}

/* ----------------------------------------------------------------------
 * blas_dgemv_()
 * ----------------------------------------------------------------------
 * Use the AMD clBlas library to do matrix vector multiplies. We require
 * that the matrix and vector are already loaded in to mem objects on the
 * device.
 *
 * AMDBLAS must be defined in the makefile
 *
 * The resulting vector will be written to a device mem obj. Someone else
 * is responsible for transferring that back to the host.
 *
 * Dgemv will do Y=aA.X+bY but we only want Y = A.X so we set a=1.0, b=0.0.
 *
 * Return 1 on success, 0 on error (or if AMDBLAS undefined)
 * ----------------------------------------------------------------------
 */
int blas_dgemv_( const int *n_rowsA,		// Rows in the matrix (and rows in the Y vector)
		 const int *n_colsA,		// Cols in the matrix (and rows in the X vector
		 const void **A_device_ptr,	// Device mem obj containing the matrix
		 const void **X_device_ptr,	// Device mem obj containing the X vector
		 void **Y_device_ptr		// Device mem obj to hold resulting Y vector
  )
{
#ifdef AMDBLAS
  cl_int err;
  cl_event event = NULL;
  cl_mem memobj_A;
  cl_mem memobj_X;
  cl_mem memobj_Y;
  
  // Blas args. Assumes Fortran array ordering.
  const clAmdBlasOrder order = clAmdBlasColumnMajor;
  const clAmdBlasTranspose transA = clAmdBlasNoTrans;
  const size_t lda = *n_rowsA;	// Leading dimension of matrix A
				// Stride between consecutive elems in a row
  const int incx = 1;		// Stride through X vector
  const int incy = 1;		// Stride through Y vector
  const double alpha = 1.0;	// Scalar coeffs
  const double beta = 0.0;

  SNAME( "blas_dgemv" );

  memobj_A = *(cl_mem *)A_device_ptr;
  memobj_X = *(cl_mem *)X_device_ptr;
  memobj_Y = *(cl_mem *)Y_device_ptr;

  err = clAmdBlasDgemv(order, transA, *n_rowsA, *n_colsA, alpha, memobj_A, lda, memobj_X, 0, incx,
		       beta, memobj_Y, 0, incy, 1, &g_oclinfo->queue, 0, NULL, &event);

  if (err != CL_SUCCESS) {
    printf("clAmdBlasDgemv() failed with %d\n", err);
    return 0;
  }
  else {
    // Wait for calculations to finish
    err = clWaitForEvents(1, &event);
    if (err != CL_SUCCESS) {
      printf("clWaitForEvents() failed with %d\n", err);
      return 0;
    }
  }

  ENAME( "blas_dgemv" );
  return 1;
#else
  return 0;
#endif
}

/* ----------------------------------------------------------------------
 * mem_copy_blas_dgemv_async_() - experimental
 * ----------------------------------------------------------------------
 * Use the AMD clBlas library to do matrix vector multiplies but try to
 * overlap calls to the library with async memory copies. We are supplied
 * the array of all matrices and vectors, none of which have been copied
 * to the device. We transfer and single matrix and vector asynchronously
 * then multiply them once the copy has finished. The read back from the
 * device occurs once the clBlas library has returned. This is also done
 * asynchronously. The next call to the clBlas library must wait until
 * this readback has finished.
 *
 * AMDBLAS must be defined in the makefile
 *
 * Dgemv will do Y=aA.X+bY but we only want Y = A.X so we set a=1.0, b=0.0.
 *
 * Return 1 on success, 0 on error (or if AMDBLAS undefined)
 * ----------------------------------------------------------------------
 */
int mem_copy_blas_dgemv_async_( const int *n_matsA,	// Total number of matrices (and X,Y vectors)
				const int *n_rowsA,	// Rows in one matrix (and rows of one Y vector)
				const int *n_colsA,	// Cols in one matrix (and rows of one X vector)
				const void *all_A_host_ptr,	// Entire array of matrices
				const void *all_X_host_ptr,	// Entire array of X vectors
				void *all_Y_host_ptr,		// Entire array of resulting Y vectors
				const void **A_device_ptr,	// Device memobj to use for one matrix
				const void **X_device_ptr,	// Device memobj to use for one X vector
				void **Y_device_ptr		// Device memobj to use for one Y vector
  )
{
#ifdef AMDBLAS
  cl_int err;
  cl_mem memobj_A;
  cl_mem memobj_X;
  cl_mem memobj_Y;
  cl_bool blockflag = CL_FALSE;

  // Create an event obj for
  // 0: write matrix to device
  // 1: write vector to device
  // 2: Multiply using clBlas
  // 3: Read result vector from device
  cl_event events[4] = { NULL, NULL, NULL, NULL };

  // Blas args. Assumes Fortran array ordering.
  const clAmdBlasOrder order = clAmdBlasColumnMajor;
  const clAmdBlasTranspose transA = clAmdBlasNoTrans;
  const size_t lda = *n_rowsA;	// Leading dimension of matrix A
  const int incx = 1;		// Stride through X vector
  const int incy = 1;		// Stride through Y vector
  const double alpha = 1.0;	// Scalar coeffs
  const double beta = 0.0;

  // Host array starting points
  const double *hA = (double *)all_A_host_ptr;
  const double *hX = (double *)all_X_host_ptr;
        double *hY = (double *)all_Y_host_ptr;

  // Offsets in to the host arrays. We assume arrays of doubles.
  const int matrix_offset = *n_rowsA * *n_colsA;
  const int x_offset = *n_colsA;
  const int y_offset = *n_rowsA;

  SNAME( "mem_copy_blas_dgemv_async" );

  // These should not have been loaded with data. We're going to do that.
  // We do want them to have been allocated on the device though by now.
  memobj_A = *(cl_mem *)A_device_ptr;
  memobj_X = *(cl_mem *)X_device_ptr;
  memobj_Y = *(cl_mem *)Y_device_ptr;


  // Loop over all of the matrices and copy each one and the associated X vector
  // to the device asynchronously. Then multiply them.
  for ( int matvec_id = 0; matvec_id < *n_matsA; matvec_id++ ) {

    // Copy current A matrix to device (asynchronously if blockflag is false)
    err  = clEnqueueWriteBuffer(g_oclinfo->queue, memobj_A, blockflag, 0, matrix_offset*sizeof(double), 
				&hA[matvec_id*matrix_offset], 0, NULL, &events[0] );
    err |= clEnqueueWriteBuffer(g_oclinfo->queue, memobj_X, blockflag, 0, x_offset*sizeof(double), 
				&hX[matvec_id*x_offset], 0, NULL, &events[1] );
    if ( err != CL_SUCCESS ) {
      printf( "\nFailed to enqueue non-blocking write of matrix A or X vector %d to GPU\n", matvec_id);
      return 0;
    }

    // Wait for read-op from previous iteration
    if ( events[3] )
      err = clWaitForEvents(1, &events[3]);

    // Assume that the kernel in here won't be called until the previous events have completed.
    err = clAmdBlasDgemv(order, transA, *n_rowsA, *n_colsA, alpha, memobj_A, lda, memobj_X, 0, incx,
			 beta, memobj_Y, 0, incy, 1, &g_oclinfo->queue, 2, events, &events[2]);
    if (err != CL_SUCCESS) {
      printf("\nclAmdBlasDgemv() failed with %d\n", err);
      return 0;
    }
    
    // Read the results once the clBlas event ([2]) has completed
    err = clEnqueueReadBuffer(g_oclinfo->queue, memobj_Y, blockflag, 0, y_offset*sizeof(double),
			      &hY[matvec_id*y_offset], 1, &events[2], &events[3] );
    if ( err != CL_SUCCESS ) {
      printf( "\nFailed to enqueue non-blocking read of Y vector %d from GPU\n", matvec_id);
      return 0;
    }
  }

  // Wait for the final read event to complete
  err = clWaitForEvents(1, &events[3]);

  // Wait for the commands to complete
  clFinish(g_oclinfo->queue);

  ENAME( "mem_copy_blas_dgemv_async_" );
  return 1;
#else
  return 0;
#endif
}

/* ----------------------------------------------------------------------
 * blas_dgemm_()
 * ----------------------------------------------------------------------
 * Use the AMD clBlas library to do matrix matrix multiplies. We require
 * that the matrices are already loaded in to mem objects on the device.
 *
 * AMDBLAS must be defined in the makefile
 *
 * The resulting matrix will be written to a device mem obj. Someone else
 * is responsible for transferring that back to the host.
 *
 * Dgemm will do C=aAB+bC but we only want C = AB so we set a=1.0, b=0.0.
 *
 * Return 1 on success, 0 on error (or if AMDBLAS undefined)
 * ----------------------------------------------------------------------
 */
int blas_dgemm_( const int *n_rowsA,		// Rows in matrix A
		 const int *n_colsB,		// Cols in matrix B
		 const int *n_colsArowsB,	// Cols in matrix A == Rows in matrix B
		 const void **A_device_ptr,	// Device mem obj containing matrix A
		 const void **B_device_ptr,	// Device mem obj containing matrix B
		 void **C_device_ptr		// Device mem obj to hold resulting C matrix
  )
{
#ifdef AMDBLAS
  cl_int err;
  cl_event event = NULL;
  cl_mem memobj_A;
  cl_mem memobj_B;
  cl_mem memobj_C;
  
  // Blas args. Assumes Fortran array ordering.
  const clAmdBlasOrder order = clAmdBlasColumnMajor;
  const clAmdBlasTranspose transAB = clAmdBlasNoTrans;
  const size_t lda = *n_rowsA;		// Leading dimension of matrix A
  const size_t ldb = *n_rowsA;		// Leading dimension of matrix B
  const size_t ldc = *n_rowsA;		// Leading dimension of matrix C
  const double alpha = 1.0;		// Scalar coeffs
  const double beta = 0.0;

  SNAME( "blas_dgemm" );

  memobj_A = *(cl_mem *)A_device_ptr;
  memobj_B = *(cl_mem *)B_device_ptr;
  memobj_C = *(cl_mem *)C_device_ptr;

  err = clAmdBlasDgemm(order, transAB, transAB, *n_rowsA, *n_colsB, *n_colsArowsB,
		       alpha, memobj_A, lda, memobj_B, ldb, beta, memobj_C, ldc,
		       1, &g_oclinfo->queue, 0, NULL, &event);

  if (err != CL_SUCCESS) {
    printf("clAmdBlasDgemm() failed with %d\n", err);
    return 0;
  }
  else {
    // Wait for calculations to finish
    err = clWaitForEvents(1, &event);
    if (err != CL_SUCCESS) {
      printf("clWaitForEvents() failed with %d\n", err);
      return 0;
    }
  }

  ENAME( "blas_dgemm" );
  return 1;
#else
  return 0;
#endif
}

/* ----------------------------------------------------------------------
 * matrix_vector_multiplies_()
 * ----------------------------------------------------------------------
 * We can't use the BLAS library to do all of the matrix/vector mults in
 * one go - i.e., when we transfer the entire array of matrices and the
 * array of vectors to the GPU in one go. The BLAS library does not allow
 * you to offset in to an array of matrices to pick out a particular matrix
 * to use during the multiplication (although you can do this for the vector
 * arg). Hence we use our own y=A.x kernel (x,y vectors, A matrix).
 *
 * The input args are an array of matrices and an array of vectors. We
 * write out another array of vectors.
 *
 * This is not the kernel code. Pass kernel args to the device and then 
 * enqueue the kernel. We'll wait for it to finish.
 *
 * The caller can then retrieve the results vector (we don't do that here).
 *
 * Args:
 * In
 *	int *num_mats: Number of matrices stored in the array of matrices
 *		       (also gives the number of vectors in the array of vecs)
 *	int *num_rows: Number of rows in ONE matrix in the array
 *		       (also gives the number of rows in ONE of the Y vectors)
 *	int *num_cols: Number of cols in ONE matrix in the array
 *		       (also gives the number of rows in ONE of the X vectors)
 *	void **A_device_ptr: Pointer to (void *) representing the A matrix's
 *			     OpenCL mem object.
 *	void **X_device_ptr: Pointer to (void *) representing the X vector's
 *			     OpenCL mem object.
 * Out
 *	void **Y_device_ptr: Pointer to (void *) representing the Y vector's
 *			     OpenCL mem object.
 *
 * Returns 1 on succes, 0 on error
 * ----------------------------------------------------------------------
 */
int matrix_vector_multiplies_(const int *num_mats, const int *num_rows, const int *num_cols,
			      const void **A_device_ptr, const void **X_device_ptr, void **Y_device_ptr )
{
  cl_int err;
  size_t global[1];                      // global domain size  
  size_t local[1];                       // local  domain size  

  SNAME( "matrix_vector_multiplies" );

  // Set the arguments to our compute kernel
  err  = clSetKernelArg(g_oclinfo->kernel, 0, sizeof(int), (const void *)num_mats);
  err |= clSetKernelArg(g_oclinfo->kernel, 1, sizeof(int), (const void *)num_rows);
  err |= clSetKernelArg(g_oclinfo->kernel, 2, sizeof(int), (const void *)num_cols);
  err |= clSetKernelArg(g_oclinfo->kernel, 3, sizeof(cl_mem),(cl_mem *)A_device_ptr );
  err |= clSetKernelArg(g_oclinfo->kernel, 4, sizeof(cl_mem),(cl_mem *)X_device_ptr );
  err |= clSetKernelArg(g_oclinfo->kernel, 5, sizeof(cl_mem),(cl_mem *)Y_device_ptr );
  if (err != CL_SUCCESS) {
    printf("Error: Failed to set kernel arguments (err %d)\n", err);
    return 0;
  }
  
  // Submit the kernel to the queue

  // Use a 1D work-group of size 60 work-items (threads). Each thread in the
  // work group will do the dot product of a row in the matrix with the vector.
  // In total we need to process every row in every matrix.
  // Note, if '&local' is replaced by NULL, OpenCL will calculate a local work group size.
  global[0] = *num_mats * *num_rows;
  local[0]  = 60;			// Must divide the global size evenly (each matrix has 60 rows)
  err = clEnqueueNDRangeKernel(g_oclinfo->queue, g_oclinfo->kernel, 1, NULL, global, local, 0, NULL, NULL);

  if (err != CL_SUCCESS) {
    printf("Error: Failed to execute kernel (err %d)\n", err);
    return 0;
  }
  
  // Wait for the commands to complete before reading back results
  clFinish(g_oclinfo->queue);
    
  ENAME( "matrix_vector_multiplies" );

  return 1;
}

/* ----------------------------------------------------------------------
 * matrix_matrix_multiplies_2d_()
 * ----------------------------------------------------------------------
 * Perform matrix-matrix multiplication: C=AB using a 2D NDRange
 *
 * Not all of the input args (rows, cols) are passed through to the kernel.
 * Instead they set the global work size, which the kernel can then retrieve
 * using the standard OpenCL functions.
 *
 * This is not the kernel code. Pass kernel args to the device and then 
 * enqueue the kernel. We'll wait for it to finish.
 *
 * The caller can then retrieve the results vector (we don't do that here).
 *
 * Args:
 * In
 *	int *num_rowsA: Number of rows in matrix A
 *	int *num_colsB: Number of cols in matrix B
 *	int *num_colsArowsB: Common number of cols in A, rows in B
 *	void **A_device_ptr: Pointer to (void *) representing the A matrix's
 *			     OpenCL mem object.
 *	void **B_device_ptr: Pointer to (void *) representing the B matrix's
 *			     OpenCL mem object.
 * Out
 *	void **C_device_ptr: Pointer to (void *) representing the C matrix's
 *			     OpenCL mem object.
 *
 * Returns 1 on succes, 0 on error
 * ----------------------------------------------------------------------
 */
int matrix_matrix_multiplies_2d_( const int *n_rowsA,		// Rows in matrix A
				  const int *n_colsB,		// Cols in matrix B
				  const int *n_colsArowsB,	// Cols in matrix A == Rows in matrix B
				  const void **A_device_ptr, 
				  const void **B_device_ptr, 
				  void **C_device_ptr )
{
  cl_int err;
  size_t global[2];	 // global domain size  
  size_t local[2];	 // local work-group size  
  size_t num_workgroups[2]; // Temp var to calculate global work size

  SNAME( "matrix_matrix_multiplies_2d" );

  // Set the arguments to our compute kernel
  err  = clSetKernelArg(g_oclinfo->kernel, 0, sizeof(int), (const void *)n_rowsA);
  err |= clSetKernelArg(g_oclinfo->kernel, 1, sizeof(int), (const void *)n_colsArowsB);
  err |= clSetKernelArg(g_oclinfo->kernel, 2, sizeof(int), (const void *)n_colsB);
  err |= clSetKernelArg(g_oclinfo->kernel, 3, sizeof(cl_mem),(cl_mem *)A_device_ptr );
  err |= clSetKernelArg(g_oclinfo->kernel, 4, sizeof(cl_mem),(cl_mem *)B_device_ptr );
  err |= clSetKernelArg(g_oclinfo->kernel, 5, sizeof(cl_mem),(cl_mem *)C_device_ptr );
  if (err != CL_SUCCESS) {
    printf("Error: Failed to set kernel arguments (err %d)\n", err);
    fflush(stdout);
    return 0;
  }

  // Work-groups of threads
  local[0]  = 16;
  local[1]  = 16;

  // Total amount of work. Must be divisible by work-group sizes.
  num_workgroups[0] = ((*n_colsB)/local[0]) + ((((*n_colsB) % local[0])==0)?0:1);
  global[0] = num_workgroups[0] * local[0];
  num_workgroups[1] = ((*n_rowsA)/local[1]) + ((((*n_rowsA) % local[1])==0)?0:1);
  global[1] = num_workgroups[1] * local[1];

  // Submit the kernel to the queue. Let OpenCL set the local work-group size.
  err = clEnqueueNDRangeKernel(g_oclinfo->queue, g_oclinfo->kernel, 2, NULL, global, local, 0, NULL, NULL);

  if (err != CL_SUCCESS) {
    printf("Error: Failed to execute kernel (err %d)\n", err);
    fflush(stdout);
    return 0;
  }

  // Wait for the commands to complete before reading back results
  clFinish(g_oclinfo->queue);

  ENAME( "matrix_matrix_multiplies_2d" );

  return 1;
}

/* ----------------------------------------------------------------------
 * matrix_matrix_multiplies_1d_()
 * ----------------------------------------------------------------------
 * Perform matrix-matrix multiplication: C=AB using a 1D NDRange
 *
 * Not all of the input args (rows, cols) are passed through to the kernel.
 * Instead they set the global work size, which the kernel can then retrieve
 * using the standard OpenCL functions.
 *
 * This is not the kernel code. Pass kernel args to the device and then 
 * enqueue the kernel. We'll wait for it to finish.
 *
 * The caller can then retrieve the results vector (we don't do that here).
 *
 * Args:
 * In
 *	int *num_rowsA: Number of rows in matrix A
 *	int *num_colsB: Number of cols in matrix B
 *	int *num_colsArowsB: Common number of cols in A, rows in B
 *	void **A_device_ptr: Pointer to (void *) representing the A matrix's
 *			     OpenCL mem object.
 *	void **B_device_ptr: Pointer to (void *) representing the B matrix's
 *			     OpenCL mem object.
 * Out
 *	void **C_device_ptr: Pointer to (void *) representing the C matrix's
 *			     OpenCL mem object.
 *
 * Returns 1 on succes, 0 on error
 * ----------------------------------------------------------------------
 */
int matrix_matrix_multiplies_1d_( const int *n_rowsA,		// Rows in matrix A
				  const int *n_colsB,		// Cols in matrix B
				  const int *n_colsArowsB,	// Cols in matrix A == Rows in matrix B
				  const void **A_device_ptr, 
				  const void **B_device_ptr, 
				  void **C_device_ptr )
{
  cl_int err;
  size_t global[1];	 // global domain size  
  size_t local[1];	 // local work-group size  
  size_t num_workgroups; // Temp var to calculate global work size

  SNAME( "matrix_matrix_multiplies_1d" );

  // Set the arguments to our compute kernel
  err  = clSetKernelArg(g_oclinfo->kernel, 0, sizeof(int), (const void *)n_rowsA);
  err |= clSetKernelArg(g_oclinfo->kernel, 1, sizeof(int), (const void *)n_colsArowsB);
  err |= clSetKernelArg(g_oclinfo->kernel, 2, sizeof(int), (const void *)n_colsB);
  err |= clSetKernelArg(g_oclinfo->kernel, 3, sizeof(cl_mem),(cl_mem *)A_device_ptr );
  err |= clSetKernelArg(g_oclinfo->kernel, 4, sizeof(cl_mem),(cl_mem *)B_device_ptr );
  err |= clSetKernelArg(g_oclinfo->kernel, 5, sizeof(cl_mem),(cl_mem *)C_device_ptr );
  if (err != CL_SUCCESS) {
    printf("Error: Failed to set kernel arguments (err %d)\n", err);
    return 0;
  }

  // Number of work-items (threads) in a work-group required by this kernel
  local[0]  = 60;

  // Total amount of work. One work-item per matrix B,C column but rounded
  // up to nearest multiple of local[0], as required by OpenCL. The kernel
  // must be able to handle (ignore) the extra threads.
  num_workgroups = ((*n_colsB)/local[0]) + ((((*n_colsB) % local[0])==0)?0:1);
  global[0] = num_workgroups * local[0];

  // Submit the kernel to the queue. Let OpenCL set the local work-group size.
  err = clEnqueueNDRangeKernel(g_oclinfo->queue, g_oclinfo->kernel, 1, NULL, global, local, 0, NULL, NULL);

  if (err != CL_SUCCESS) {
    printf("Error: Failed to execute kernel (err %d)\n", err);
    return 0;
  }

  // Wait for the commands to complete before reading back results
  clFinish(g_oclinfo->queue);

  ENAME( "matrix_matrix_multiplies_1d" );

  return 1;
}

/* ----------------------------------------------------------------------
 * c_compare_vecs()
 * ----------------------------------------------------------------------
 * Compare two 'double' vectors, to within a built-in tolerance.
 *
 * Args:
 * In
 *	int *dim:	Number of elements in both arrays
 *	double vec_a[]: Vector to compare
 *	double vec_b[]: Vector to compare
 *
 * Returns 1 if considered equal, 0 otherwise
 * ----------------------------------------------------------------------
 */
int c_compare_vecs_( const int *dim, const double vec_a[], const double vec_b[] )
{
  int cdim = *dim;
  int r;

  for ( r=0; r<cdim; r++ ) {
    if ( fabs(vec_a[r] - vec_b[r]) > (double)0.00001 ) {
      printf( "[%d]: %lf != %lf\n", r, vec_a[r], vec_b[r] );
      return 0;
    }
  }

  return 1;
}
