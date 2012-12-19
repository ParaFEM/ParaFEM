# MAKEFILE src/programs/5th_ed/xx2/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles xx2 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx2 executable as RELEASE edition
#   debug: Build xx2 executable as DEBUG edition
#   install: Copy xx2 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx2 executable
#   relink: Rebuild xx2 from object files and libraries
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
	@echo "Building xx2 RELEASE"
	$(MAKE) xx2 \
	FFLAGS="-c -I$(MODULE) $(XX2_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX2_MPI_REL_LDFLAGS)"
	@echo "Done xx2 RELEASE"

debug:
	@echo
	@echo "Building xx2 DEBUG"
	$(MAKE) xx2 \
	FFLAGS="-c -I$(MODULE) $(XX2_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX2_MPI_DBG_LDFLAGS)"
	@echo "Done xx2 DEBUG"

xx2: xx2.o
	$(FC) xx2.o -o xx2 $(LDFLAGS)

xx2.o: xx2.f90
	$(FC) xx2.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx2 ../../../../bin/xx2

relink: xx2

install:
	cp xx2 ../../../../bin
