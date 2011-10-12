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

# No env var is available even if ATI/CUDA modules are loaded.
# Ensure that the correct module *is* loaded in the job script
# when running the xx4 code so that LD_LIBRARY_PATH points to
# the required libOpenCL.so which we dynamically link in here.

# AMD OpenCL
#OCLROOT     = /export/ati-app-sdk/lnx64/2.5
#OCLBLASROOT = $(HOME)/clAmdBlas-1.4.182
#OCLBLASLIB  = -L$(OCLBLASROOT)/lib64 -lclAmdBlas
#OCLINC      = -I$(OCLROOT)/include -I$(OCLBLASROOT)/include
#OCLLIB      = -L$(OCLROOT)/lib/x86_64/ -lOpenCL $(OCLBLASLIB)
#OCLFLAGS    = -DAMDBLAS

# Nvidia OpenCL
OCLROOT  = /export/cuda/cuda3.2.16
OCLINC   = -I$(OCLROOT)/include
OCLLIB   = -lOpenCL

# C compiler flags not set in parafem included makefiles
CC=mpicc
XX4_MPI_REL_CFLAGS=-O3
XX4_MPI_DBG_CFLAGS=-g

COBJS = opencl_funcs.o
FOBJS = xx4.o
OBJS = $(COBJS) $(FOBJS)
EXES = xx4

all:	release

release:
	@echo
	@echo "Building $(EXES) RELEASE"
	$(MAKE) $(EXES) \
	FFLAGS="-c -I$(MODULE) $(XX4_MPI_REL_FFLAGS)" \
	CFLAGS="-c $(XX4_MPI_REL_CFLAGS) $(OCLINC) $(OCLFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX4_MPI_REL_LDFLAGS) $(OCLLIB) -lm"
	@echo "Done $(EXES) RELEASE"

debug:
	@echo
	@echo "Building $(EXES) DEBUG"
	$(MAKE) $(EXES) \
	FFLAGS="-c -I$(MODULE) $(XX4_MPI_DBG_FFLAGS)" \
	CFLAGS="-c $(XX4_MPI_DBG_CFLAGS) $(OCLINC) $(OCLFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX4_MPI_DBG_LDFLAGS) $(OCLLIB) -lm"
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
