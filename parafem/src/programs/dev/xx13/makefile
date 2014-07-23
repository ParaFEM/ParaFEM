# MAKEFILE src/programs/5th_ed/xx13/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles xx13 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx13 executable as RELEASE edition
#   debug: Build xx13 executable as DEBUG edition
#   install: Copy xx13 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx13 executable
#   relink: Rebuild xx13 from object files and libraries
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
	@echo "Building XX13 RELEASE"
	$(MAKE) xx13 \
	FFLAGS="-c -I$(MODULE) $(XX13_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX13_MPI_REL_LDFLAGS)"
	@echo "Done xx13 RELEASE"

debug:
	@echo
	@echo "Building xx13 DEBUG"
	$(MAKE) xx13 \
	FFLAGS="-c -I$(MODULE) $(XX13_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX13_MPI_DBG_LDFLAGS)"
	@echo "Done xx13 DEBUG"

xx13: xx13.o
	$(FC) xx13.o -o xx13 $(LDFLAGS)

xx13.o: xx13.f90
	$(FC) xx13.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx13 ../../../../bin/xx13

relink: xx13

install:
	cp xx13 ../../../../bin

test:
	./test.sh xx13
	
