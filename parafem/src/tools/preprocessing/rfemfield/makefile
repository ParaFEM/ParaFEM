# MAKEFILE src/tools/preprocessing/rfemfield/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles rfemfield program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build rfemfield executable as RELEASE edition
#   debug: Build rfemfield executable as DEBUG edition
#   install: Copy rfemfield executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local rfemfield executable
#   relink: Rebuild rfemfield from object files and libraries
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
	@echo "Building rfemfield RELEASE"
	$(MAKE) rfemfield \
	FFLAGS="-c -I$(MODULE) $(RFEMFIELD_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(RFEMFIELD_MPI_REL_LDFLAGS) \
	-l gaf77"
	@echo "Done rfemfield RELEASE"

debug:
	@echo
	@echo "Building rfemfield DEBUG"
	$(MAKE) rfemfield \
	FFLAGS="-c -I$(MODULE) $(RFEMFIELD_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(RFEMFIELD_MPI_DBG_LDFLAGS) \
	         -l gaf77"
	@echo "Done rfemfield DEBUG"

rfemfield: rfemfield.o
	$(FC) rfemfield.o -o rfemfield $(LDFLAGS) 

rfemfield.o: rfemfield.f90
	$(FC) rfemfield.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f rfemfield ../../../../bin/rfemfield

relink: rfemfield

install:
	$(INSTALL_CMD) rfemfield ../../../../bin
