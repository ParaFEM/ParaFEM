# MAKEFILE src/programs/5th_ed/xx4/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles xx4 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx4 executable as RELEASE edition
#   debug: Build xx4 executable as DEBUG edition
#   install: Copy xx4 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx4 executable
#   relink: Rebuild xx4 from object files and libraries
# -----------------------------------------------------------------------
# Notes:
# -----------------------------------------------------------------------

include ../../../../build/mk_defs.inc
include ../../../../build/$(MACHINE).inc

MODULE	= ../../../../include/mpi
LIB	= ../../../../lib
.SUFFIXES: .o .f90 .c

# C compiler flags not set in parafem included makefiles
XX4_MPI_REL_CFLAGS=-O3
XX4_MPI_DBG_CFLAGS=-g

COBJS = opencl_funcs.o
EXES = xx4

all:	release

release:
	@echo
	@echo "Building $(EXES) RELEASE"
	$(MAKE) $(EXES) \
	FFLAGS="-c -I$(MODULE) $(XX4_MPI_REL_FFLAGS)" \
	CFLAGS="-c $(XX4_MPI_REL_CFLAGS) $(GPU_OCL_INC) $(GPU_OCL_FLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX4_MPI_REL_LDFLAGS) $(GPU_OCL_LIB)"
	@echo "Done $(EXES) RELEASE"

debug:
	@echo
	@echo "Building $(EXES) DEBUG"
	$(MAKE) $(EXES) \
	FFLAGS="-c -I$(MODULE) $(XX4_MPI_DBG_FFLAGS)" \
	CFLAGS="-c $(XX4_MPI_DBG_CFLAGS) $(GPU_OCL_INC) $(GPU_OCL_FLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX4_MPI_DBG_LDFLAGS) $(GPU_OCL_LIB)"
	@echo "Done $(EXES) DEBUG"

xx4: $(COBJS) xx4.o
	$(FC) xx4.o $(COBJS) -o $@ $(LDFLAGS)


# Fortran files
xx4.o: xx4.f90
	$(FC) $(FFLAGS) $<


# OpenCL C files
opencl_funcs.o: opencl_funcs.c
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f *.o

execlean:
	rm -f $(EXES) ../../../../bin/xx4

relink: $(EXES)

install:
	cp $(EXES) ../../../../bin
