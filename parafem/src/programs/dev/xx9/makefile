# MAKEFILE src/programs/5th_ed/xx9/makefile
# Author: Louise M. Lever; Lee Margetts
# -----------------------------------------------------------------------
# Compiles xx9 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx9 executable as RELEASE edition
#   debug: Build xx9 executable as DEBUG edition
#   install: Copy xx9 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx9 executable
#   relink: Rebuild xx9 from object files and libraries
# -----------------------------------------------------------------------
# Notes:
# -----------------------------------------------------------------------

include ../../../../build/mk_defs.inc
include ../../../../build/$(MACHINE).inc 

# Add CUBLAS library
XX9_GPU_CU_LIB = $(GPU_CU_LIB) -lcublas

# CUDA CUBLAS wrappers for Fortran/C interop.
XX9_GPU_CU_FTN_C = $(GPU_CU_HEADER)/fortran.c
XX9_GPU_CU_FTN_H = $(GPU_CU_HEADER)/fortran.h
XX9_GPU_CU_INC   = $(GPU_CU_INC) -I$(GPU_CU_ROOT)

MODULE	= ../../../../include/mpi
LIB	= ../../../../lib

.SUFFIXES: .o .f90

%.o: %.mod

.PHONY: fortranfiles

all:	release

release:
	@echo
	@echo "Building xx9 RELEASE"
	$(MAKE) xx9 \
	FFLAGS="-c -I$(MODULE) $(XX9_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX9_MPI_REL_LDFLAGS) $(XX9_GPU_CU_LIB)"
	@echo "Done xx9 RELEASE"

debug:
	@echo
	@echo "Building xx9 DEBUG"
	$(MAKE) xx9 \
	FFLAGS="-c -I$(MODULE) $(XX9_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX9_MPI_DBG_LDFLAGS) $(XX9_GPU_CU_LIB)"
	@echo "Done xx9 DEBUG"

mathsgpu.o: mathsgpu.f90
	$(FC) -c mathsgpu.f90 $(FFLAGS)

cuda_helpers.o: cuda_helpers.cu
	$(GPU_CU) -c cuda_helpers.cu $(GPU_CU_FLAGS) $(XX9_GPU_CU_INC) 

fortran.o: $(XX9_GPU_CU_FTN_C) $(XX9_GPU_CU_FTN_H)
	$(CC) -DCUBLAS_GFORTRAN -c $(XX9_GPU_CU_FTN_C) $(CFLAGS) $(XX9_GPU_CU_INC)

xx9: mathsgpu.o xx9.o cuda_helpers.o fortran.o
	$(FC) xx9.o mathsgpu.o cuda_helpers.o fortran.o -o xx9 $(LDFLAGS) $(XX9_GPU_CU_INC)

xx9.o: xx9.f90
	$(FC) xx9.f90 $(FFLAGS)

#mathsgpu.o: mathsgpu.f90
#$(FC) mathsgpu.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx9 ../../../../bin/xx9

relink: xx9

install:
	cp xx9 ../../../../bin
