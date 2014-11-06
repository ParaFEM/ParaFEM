# MAKEFILE src/programs/5th_ed/xx15/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles xx15 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx15 executable as RELEASE edition
#   debug: Build xx15 executable as DEBUG edition
#   install: Copy xx15 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx15 executable
#   relink: Rebuild xx15 from object files and libraries
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
	@echo "Building xx15 RELEASE"
	$(MAKE) xx15 \
	FFLAGS="-c -I$(MODULE) $(XX15_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX15_MPI_REL_LDFLAGS)"
	@echo "Done xx15 RELEASE"

debug:
	@echo
	@echo "Building xx15 DEBUG"
	$(MAKE) xx15 \
	FFLAGS="-c -I$(MODULE) $(XX15_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX15_MPI_DBG_LDFLAGS)"
	@echo "Done xx15 DEBUG"

xx15: xx15.o
	$(FC) xx15.o -o xx15 $(LDFLAGS)

xx15.o: xx15.f90
	$(FC) xx15.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx15 ../../../../bin/xx15

relink: xx15

install:
	cp xx15 ../../../../bin
