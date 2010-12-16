MODULES

The modules held in this directory and subdirectories will be used to create
several libraries that will run on different architectures, as follows:

  lparafem-mpi.a       Parallel execution using MPI
  lparafem-serial.a    Serial execution
  lparafem-omp.a       Parallel execution using OpenMP
  lparafem-gpu.a       Parallel execution using GPUs
  

DIRECTORY STRUCTURE

Modules that are common to the MPI, SERIAL, OMP and GPU libraries are contained 
in the folder:

  /parafem/src/modules/shared

All other modules are contained in:

  /parafem/src/modules/serial
  /parafem/src/modules/mpi
  /parafem/src/modules/omp
  /parafem/src/modules/gpu
  
Each library will contain the shared modules and the architecture specific 
modules. The MPI, SERIAL, OMP and GPU folders will contain modules with the 
same names as each other. The function and subroutines held within the modules 
will also use the same names. Internally, the executable statements will be 
those required by the different platforms.

The rationale behind this is that the user will be able to compile essentially
the same driver program (p121, p122 etc) for different platforms without having
to modify USE statements, subroutine calls or subroutine arguments within the
driver program.

When the driver program is compiled, the user must be careful to ensure that
the compiler links to the correct library for the platform.


AUTHOR

lee.margetts@manchester.ac.uk