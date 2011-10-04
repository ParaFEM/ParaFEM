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

MODULE	= ../../../../include/mpi
LIB	= ../../../../lib -L/export/cuda/cuda4.0.17/lib64 -lcublas
.SUFFIXES: .o .f90 .c

CUBLAS_INCLUDE=-I/export/cuda/cuda4.0.17/include
CFLAGS=-O3
CC=nvcc

all:	release

release:
	@echo
	@echo "Building xx9 RELEASE"
	$(MAKE) xx9 \
	FFLAGS="-c -I$(MODULE) $(XX9_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX9_MPI_REL_LDFLAGS)"
	@echo "Done xx9 RELEASE"

debug:
	@echo
	@echo "Building xx9 DEBUG"
	$(MAKE) xx9 \
	FFLAGS="-c -I$(MODULE) $(XX9_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX9_MPI_DBG_LDFLAGS)"
	@echo "Done xx9 DEBUG"

fortran.o: fortran.c fortran.h
	$(CC) -c fortran.c $(CFLAGS) $(CUBLAS_INCLUDE)

xx9: xx9.o fortran.o
	$(FC) xx9.o fortran.o -o xx9 $(LDFLAGS)

xx9.o: xx9.f90
	$(FC) xx9.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx9 ../../../../bin/xx9

relink: xx9

install:
	cp xx9 ../../../../bin
