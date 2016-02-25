# MAKEFILE src/programs/rfem/rfemsolve/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles rfemsolve program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build rfemsolve executable as RELEASE edition
#   debug: Build rfemsolve executable as DEBUG edition
#   install: Copy rfemsolve executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local rfemsolve executable
#   relink: Rebuild rfemsolve from object files and libraries
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
	@echo "Building rfemsolve RELEASE"
	$(MAKE) rfemsolve \
	FFLAGS="-c -I$(MODULE) $(RFEMSOLVE_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(RFEMSOLVE_MPI_REL_LDFLAGS)"
	@echo "Done rfemsolve RELEASE"

debug:
	@echo
	@echo "Building rfemsolve DEBUG"
	$(MAKE) rfemsolve \
	FFLAGS="-c -I$(MODULE) $(RFEMSOLVE_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(RFEMSOLVE_MPI_DBG_LDFLAGS)"
	@echo "Done rfemsolve DEBUG"

rfemsolve: rfemsolve.o
	$(FC) rfemsolve.o -o rfemsolve $(LDFLAGS)

rfemsolve.o: rfemsolve.f90
	$(FC) rfemsolve.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f rfemsolve ../../../../bin/rfemsolve

relink: rfemsolve

install:
	cp rfemsolve ../../../../bin
