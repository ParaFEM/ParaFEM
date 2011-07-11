# MAKEFILE src/tools/postprocessing/mat2numpe/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles mat2numpe program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build mat2numpe executable as RELEASE edition
#   debug: Build mat2numpe executable as DEBUG edition
#   install: Copy mat2numpe executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local mat2numpe executable
#   relink: Rebuild mat2numpe from object files and libraries
# -----------------------------------------------------------------------
# Notes:
# -----------------------------------------------------------------------

include ../../../../build/mk_defs.inc
include ../../../../build/$(MACHINE).inc

MODULE	= ../../../../include/mpi
LIB	= ../../../../lib

.SUFFIXES: .o .f90

all:	release

release:
	@echo
	@echo "Building mat2numpe RELEASE"
	$(MAKE) mat2numpe \
	FFLAGS="-c -I$(MODULE) $(MAT2NUMPE_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(MAT2NUMPE_MPI_REL_LDFLAGS)"
	@echo "Done mat2numpe RELEASE"

debug:
	@echo
	@echo "Building mat2numpe DEBUG"
	$(MAKE) mat2numpe \
	FFLAGS="-c -I$(MODULE) $(MAT2NUMPE_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(MAT2NUMPE_MPI_DBG_LDFLAGS)"
	@echo "Done mat2numpe DEBUG"

mat2numpe: mat2numpe.o
	$(FC) mat2numpe.o -o mat2numpe $(LDFLAGS)

mat2numpe.o: mat2numpe.f90
	$(FC) mat2numpe.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f mat2numpe ../../../../bin/mat2numpe

relink: mat2numpe

install:
	cp mat2numpe ../../../../bin
