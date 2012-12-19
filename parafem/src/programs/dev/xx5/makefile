# MAKEFILE src/programs/5th_ed/xx5/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles xx5 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx5 executable as RELEASE edition
#   debug: Build xx5 executable as DEBUG edition
#   install: Copy xx5 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx5 executable
#   relink: Rebuild xx5 from object files and libraries
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
	@echo "Building xx5 RELEASE"
	$(MAKE) xx5 \
	FFLAGS="-c -I$(MODULE) $(XX5_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX5_MPI_REL_LDFLAGS)"
	@echo "Done xx5 RELEASE"

debug:
	@echo
	@echo "Building xx5 DEBUG"
	$(MAKE) xx5 \
	FFLAGS="-c -I$(MODULE) $(XX5_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX5_MPI_DBG_LDFLAGS)"
	@echo "Done xx5 DEBUG"

xx5: xx5.o
	$(FC) xx5.o -o xx5 $(LDFLAGS)

xx5.o: xx5.f90
	$(FC) xx5.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx5 ../../../../bin/xx5

relink: xx5

install:
	cp xx5 ../../../../bin
