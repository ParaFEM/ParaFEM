# MAKEFILE src/tools/preprocessing/rfemcube/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles rfemcube program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build rfemcube executable as RELEASE edition
#   debug: Build rfemcube executable as DEBUG edition
#   install: Copy rfemcube executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local rfemcube executable
#   relink: Rebuild rfemcube from object files and libraries
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
	@echo "Building rfemcube RELEASE"
	$(MAKE) rfemcube \
	FFLAGS="-c -I$(MODULE) $(RFEMCUBE_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(RFEMCUBE_MPI_REL_LDFLAGS) \
	-l gaf77"
	@echo "Done rfemcube RELEASE"

debug:
	@echo
	@echo "Building rfemcube DEBUG"
	$(MAKE) rfemcube \
	FFLAGS="-c -I$(MODULE) $(RFEMCUBE_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(RFEMCUBE_MPI_DBG_LDFLAGS) \
	         -l gaf77"
	@echo "Done rfemcube DEBUG"

rfemcube: rfemcube.o
	$(FC) rfemcube.o -o rfemcube $(LDFLAGS) 

rfemcube.o: rfemcube.f90
	$(FC) rfemcube.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f rfemcube ../../../../bin/rfemcube

relink: rfemcube

install:
	$(INSTALL_CMD) rfemcube ../../../../bin
