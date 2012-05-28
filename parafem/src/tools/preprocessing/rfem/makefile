# MAKEFILE src/tools/preprocessing/rfem/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles rfem program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build rfem executable as RELEASE edition
#   debug: Build rfem executable as DEBUG edition
#   install: Copy rfem executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local rfem executable
#   relink: Rebuild rfem from object files and libraries
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
	@echo "Building rfem RELEASE"
	$(MAKE) rfem \
	FFLAGS="-c -I$(MODULE) $(RFEM_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(RFEM_MPI_REL_LDFLAGS) \
	-l gaf77"
	@echo "Done rfem RELEASE"

debug:
	@echo
	@echo "Building rfem DEBUG"
	$(MAKE) rfem \
	FFLAGS="-c -I$(MODULE) $(RFEM_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(RFEM_MPI_DBG_LDFLAGS) \
	         -l gaf77"
	@echo "Done rfem DEBUG"

rfem: rfem.o
	$(FC) rfem.o -o rfem $(LDFLAGS) 

rfem.o: rfem.f90
	$(FC) rfem.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f rfem ../../../../bin/rfem

relink: rfem

install:
	cp rfem ../../../../bin
