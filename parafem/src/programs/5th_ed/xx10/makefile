# MAKEFILE src/programs/5th_ed/xx10/makefile
# Author: Louise M. Lever; Lee Margetts
# -----------------------------------------------------------------------
# Compiles xx10 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx10 executable as RELEASE edition
#   debug: Build xx10 executable as DEBUG edition
#   install: Copy xx10 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx10 executable
#   relink: Rebuild xx10 from object files and libraries
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
	@echo "Building xx10 RELEASE"
	$(MAKE) xx10 \
	FFLAGS="-c -I$(MODULE) $(XX10_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX10_MPI_REL_LDFLAGS)"
	@echo "Done xx10 RELEASE"

debug:
	@echo
	@echo "Building xx10 DEBUG"
	$(MAKE) xx10 \
	FFLAGS="-c -I$(MODULE) $(XX10_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX10_MPI_DBG_LDFLAGS)"
	@echo "Done xx10 DEBUG"

xx10: xx10.o
	$(FC) xx10.o -o xx10 $(LDFLAGS)

xx10.o: xx10.f90
	$(FC) xx10.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx10 ../../../../bin/xx10

relink: xx10

install:
	cp xx10 ../../../../bin
