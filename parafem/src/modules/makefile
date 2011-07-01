# MAKEFILE src/modules/makefile
# Author: Louise M. Lever
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

include ../../include/mk_defs.inc
include ../../include/$(MACHINE).inc

all:	release-mpi release-bem_lib release-bem_lib_p

release-mpi:
	@echo
	@echo "Building MODULE RELEASE (MPI)"
#	@cd mpi && $(MAKE) release && cd ..
	@echo "Done MODULE RELEASE (MPI)"

debug-mpi:
	@echo
	@echo "Building MODULE DEBUG (MPI)"
#	@cd mpi && $(MAKE) debug && cd ..
	@echo "Done MODULE DEBUG (MPI)"

release-bem_lib:
	@echo
	@echo "Building MODULE RELEASE (BEM_LIB)"
#	@cd bem && $(MAKE) release && cd ..
	@echo "Done MODULE RELEASE (BEM_LIB)"

debug-bem_lib:
	@echo
	@echo "Building MODULE DEBUG (BEM_LIB)"
#	@cd bem & $(MAKE) debug && cd ..
	@echo "Done MODULE DEBUG (BEM_LIB)"

release-bem_lib_p:
	@echo
	@echo "Building MODULE RELEASE (BEM_LIB_P)"
#	@cd bem_p && $(MAKE) release && cd ..
	@echo "Done MODULE RELEASE (BEM_LIB_P)"

debug-bem_lib_p:
	@echo
	@echo "Building MODULE DEBUG (BEM_LIB_P)"
#	@cd bem_p & $(MAKE) debug && cd ..
	@echo "Done MODULE DEBUG (BEM_LIB_P)"

clean:
	@cd mpi && $(MAKE) clean && cd ..
	@cd bem && $(MAKE) clean && cd ..
	@cd bem_p && $(MAKE) clean && cd ..

install-mpi:
	@cd mpi && $(MAKE) install & cd ..

install-bem_lib:
	@cd bem && $(MAKE) install && cd ..

install-bem_lib_p:
	@cd bem_p && $(MAKE) install && cd ..
