# MAKEFILE src/modules/makefile
# Author: Louise M. Lever, Mark Filipiak
# -----------------------------------------------------------------------
# Compiles common "shared" files and specific "mpi" files
# Builds as a STATIC library libParaFEM_mpi.a
# -----------------------------------------------------------------------
# Options:
# (default) : Compile and build library
# install: Copy modules into global include directory
#          Copy library into global lib directory
# -----------------------------------------------------------------------
# Notes:
# (1) VPATH required to find "shared" modules when building lib rule
# -----------------------------------------------------------------------

VPATH = shared
.SUFFIXES: .f90 .o 
%.o : %.mod

include ../../build/mk_defs.inc
include ../../build/$(MACHINE).inc

all:	release-mpi release-choose_solvers release-parafem_petsc release-serial release-bem_lib release-bem_lib_p

release-mpi:
	@echo
	@echo "Building MODULE RELEASE (MPI)"
	@cd mpi && $(MAKE) release && cd ..
	@echo "Done MODULE RELEASE (MPI)"

debug-mpi:
	@echo
	@echo "Building MODULE DEBUG (MPI)"
	@cd mpi && $(MAKE) debug && cd ..
	@echo "Done MODULE DEBUG (MPI)"

release-choose_solvers:
	@echo
	@echo "Building MODULE RELEASE (CHOOSE_SOLVERS)"
	@cd choose_solvers && $(MAKE) release && cd ..
	@echo "Done MODULE RELEASE (CHOOSE_SOLVERS)"

debug-choose_solvers:
	@echo
	@echo "Building MODULE DEBUG (CHOOSE_SOLVERS)"
	@cd choose_solvers && $(MAKE) debug && cd ..
	@echo "Done MODULE DEBUG (CHOOSE_SOLVERS)"

release-parafem_petsc: release-mpi
	@echo
	@echo "Building MODULE RELEASE (PARAFEM_PETSC)"
	@cd parafem_petsc && $(MAKE) release && cd ..
	@echo "Done MODULE RELEASE (PARAFEM_PETSC)"

debug-parafem_petsc: debug-mpi
	@echo
	@echo "Building MODULE DEBUG (PARAFEM_PETSC)"
	@cd parafem_petsc && $(MAKE) debug && cd ..
	@echo "Done MODULE DEBUG (PARAFEM_PETSC)"

release-serial:
	@echo
	@echo "Building MODULE RELEASE (SERIAL)"
	@cd serial && $(MAKE) release && cd ..
	@echo "Done MODULE RELEASE (SERIAL)"

debug-serial:
	@echo
	@echo "Building MODULE DEBUG (SERIAL)"
	@cd serial && $(MAKE) debug && cd ..
	@echo "Done MODULE DEBUG (SERIAL)"

release-bem_lib:
	@echo
	@echo "Building MODULE RELEASE (BEM_LIB)"
	@cd bem && $(MAKE) release && cd ..
	@echo "Done MODULE RELEASE (BEM_LIB)"

debug-bem_lib:
	@echo
	@echo "Building MODULE DEBUG (BEM_LIB)"
	@cd bem && $(MAKE) debug && cd ..
	@echo "Done MODULE DEBUG (BEM_LIB)"

release-bem_lib_p:
	@echo
	@echo "Building MODULE RELEASE (BEM_LIB_P)"
	@cd bem_p && $(MAKE) release && cd ..
	@echo "Done MODULE RELEASE (BEM_LIB_P)"

debug-bem_lib_p:
	@echo
	@echo "Building MODULE DEBUG (BEM_LIB_P)"
	@cd bem_p && $(MAKE) debug && cd ..
	@echo "Done MODULE DEBUG (BEM_LIB_P)"

clean:
	@cd mpi && $(MAKE) clean && cd ..
	@cd choose_solvers && $(MAKE) clean && cd ..
	@cd parafem_petsc && $(MAKE) clean && cd ..
	@cd serial && $(MAKE) clean && cd ..
	@cd bem && $(MAKE) clean && cd ..
	@cd bem_p && $(MAKE) clean && cd ..

execlean:
	@cd mpi && $(MAKE) execlean && cd ..
	@cd choose_solvers && $(MAKE) execlean && cd ..
	@cd parafem_petsc && $(MAKE) execlean && cd ..
	@cd serial && $(MAKE) execlean && cd ..
	@cd bem && $(MAKE) execlean && cd ..
	@cd bem_p && $(MAKE) execlean && cd ..

install:
	@echo
	@echo "Installing: MPI CHOOSE_SOLVERS PARAFEM_PETSC SERIAL BEM BEM_P libraries"
	$(MAKE) install-mpi install-choose_solvers install-parafem_petsc install-serial install-bem_lib install-bem_lib_p

install-mpi:
	@cd mpi && $(MAKE) install & cd ..

install-choose_solvers:
	@cd choose_solvers && $(MAKE) install & cd ..

install-parafem_petsc:
	@cd parafem_petsc && $(MAKE) install & cd ..

install-serial:
	@cd serial && $(MAKE) install & cd ..

install-bem_lib:
	@cd bem && $(MAKE) install && cd ..

install-bem_lib_p:
	@cd bem_p && $(MAKE) install && cd ..
