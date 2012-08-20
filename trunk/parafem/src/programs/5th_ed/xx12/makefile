# MAKEFILE src/programs/5th_ed/xx12/makefile
# Author: Louise M. Lever
#         Llion M. Evans
# -----------------------------------------------------------------------
# Compiles xx12 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx12 executable as RELEASE edition
#   debug: Build xx12 executable as DEBUG edition
#   install: Copy xx12 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx12 executable
#   relink: Rebuild xx12 from object files and libraries
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
	@echo "Building xx12 RELEASE"
	$(MAKE) xx12 \
	FFLAGS="-c -I$(MODULE) $(XX12_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX12_MPI_REL_LDFLAGS)"
	@echo "Done xx12 RELEASE"

debug:
	@echo
	@echo "Building xx12 DEBUG"
	$(MAKE) xx12 \
	FFLAGS="-c -I$(MODULE) $(XX12_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX12_MPI_DBG_LDFLAGS)"
	@echo "Done xx12 DEBUG"

xx12: xx12.o
	$(FC) xx12.o -o xx12 $(LDFLAGS)

xx12.o: xx12.f90
	$(FC) xx12.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx12 ../../../../bin/xx12

relink: xx12

install:
	cp xx12 ../../../../bin
