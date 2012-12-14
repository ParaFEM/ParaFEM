# MAKEFILE src/tools/preprocessing/rfembc/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles rfembc program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build rfembc executable as RELEASE edition
#   debug: Build rfembc executable as DEBUG edition
#   install: Copy rfembc executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local rfembc executable
#   relink: Rebuild rfembc from object files and libraries
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
	@echo "Building rfembc RELEASE"
	$(MAKE) rfembc \
	FFLAGS="-c -I$(MODULE) $(RFEMBC_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(RFEMBC_MPI_REL_LDFLAGS) \
	-lgaf77"
	@echo "Done rfembc RELEASE"

debug:
	@echo
	@echo "Building rfembc DEBUG"
	$(MAKE) rfembc \
	FFLAGS="-c -I$(MODULE) $(RFEMBC_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(RFEMBC_MPI_DBG_LDFLAGS) \
	         -lgaf77"
	@echo "Done rfembc DEBUG"

rfembc: rfembc.o
	$(FC) rfembc.o -o rfembc $(LDFLAGS) 

rfembc.o: rfembc.f90
	$(FC) rfembc.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f rfembc ../../../../bin/rfembc

relink: rfembc

install:
	$(INSTALL_CMD) rfembc ../../../../bin
