MPI_STUBS

MPI_STUBS are FORTRAN90 routines that can look enough like the real MPI routines to satisfy the compiler and the linker/loader. For very simple MPI programs that can run correctly on one process, the resulting executable program might even run, and might even be correct. However, the main purpose of this library is simple to make it possible to compile and debug MPI programs on systems that do not actually support MPI. 

If your compiler supports MPI, then you may not need to copy the include file. Otherwise, copy the include file and the source code of the library, rename the include file to match the name of the include file you are invoking, and compile the library along with your program. 

These MPI_STUBS were downloaded from the following website on 13th June 2012:

http://people.sc.fsu.edu/~jburkardt/classes/acs2_2008/mpi/mpi_stubs/mpi_stubs.html


USAGE IN PARAFEM

Initially, we are trying to build parts of ParaFEM under Windows for a particular project. At some point, this capability will be better integrated into ParaFEM.