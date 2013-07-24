# MAKEFILE src/programs/rfem/rfemsolve_te/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles rfemsolve_te program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build rfemsolve_te executable as RELEASE edition
#   debug: Build rfemsolve_te executable as DEBUG edition
#   install: Copy rfemsolve_te executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local rfemsolve_te executable
#   relink: Rebuild rfemsolve_te from object files and libraries
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
	@echo "Building rfemsolve_te RELEASE"
	$(MAKE) rfemsolve_te \
	FFLAGS="-c -I$(MODULE) $(RFEMSOLVE_TE_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(RFEMSOLVE_TE_MPI_REL_LDFLAGS)"
	@echo "Done rfemsolve_te RELEASE"

debug:
	@echo
	@echo "Building rfemsolve_te DEBUG"
	$(MAKE) rfemsolve_te \
	FFLAGS="-c -I$(MODULE) $(RFEMSOLVE_TE_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(RFEMSOLVE_TE_MPI_DBG_LDFLAGS)"
	@echo "Done rfemsolve_te DEBUG"

rfemsolve_te: rfemsolve_te.o
	$(FC) rfemsolve_te.o -o rfemsolve_te $(LDFLAGS)

rfemsolve_te.o: rfemsolve_te.f90
	$(FC) rfemsolve_te.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f rfemsolve_te ../../../../bin/rfemsolve_te

relink: rfemsolve_te

install:
	cp rfemsolve_te ../../../../bin
