# MAKEFILE src/programs/rfem/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles rfem-solve program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build rfem-solve executable as RELEASE edition
#   debug: Build rfem-solve executable as DEBUG edition
#   install: Copy rfem-solve executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local rfem-solve executable
#   relink: Rebuild rfem-solve from object files and libraries
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
	@echo "Building rfem-solve RELEASE"
	$(MAKE) rfem-solve \
	FFLAGS="-c -I$(MODULE) $(RFEM-SOLVE_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) \
	$(RFEM-SOLVE_MPI_REL_LDFLAGS)"
	@echo "Done rfem-solve RELEASE"

debug:
	@echo
	@echo "Building rfem-solve DEBUG"
	$(MAKE) rfem-solve \
	FFLAGS="-c -I$(MODULE) $(RFEM-SOLVE_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) \ 
	$(RFEM-SOLVE_MPI_DBG_LDFLAGS)"
	@echo "Done rfem-solve DEBUG"

rfem-solve: rfem-solve.o
	$(FC) rfem-solve.o -o rfem-solve $(LDFLAGS)

rfem-solve.o: rfem-solve.f90
	$(FC) rfem-solve.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f rfem-solve ../../../../bin/rfem-solve

relink: rfem-solve

install:
	cp rfem-solve ../../../../bin
