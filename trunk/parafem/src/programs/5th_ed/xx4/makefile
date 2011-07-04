# MAKEFILE src/programs/5th_ed/xx4/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles xx4 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx4 executable as RELEASE edition
#   debug: Build xx4 executable as DEBUG edition
#   install: Copy xx4 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx4 executable
#   relink: Rebuild xx4 from object files and libraries
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
	@echo "Building xx4 RELEASE"
	$(MAKE) xx4 \
	FFLAGS="-c -I$(MODULE) $(XX4_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX4_MPI_REL_LDFLAGS)"
	@echo "Done xx4 RELEASE"

debug:
	@echo
	@echo "Building xx4 DEBUG"
	$(MAKE) xx4 \
	FFLAGS="-c -I$(MODULE) $(XX4_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX4_MPI_DBG_LDFLAGS)"
	@echo "Done xx4 DEBUG"

xx4: xx4.o
	$(FC) xx4.o -o xx4 $(LDFLAGS)

xx4.o: xx4.f90
	$(FC) xx4.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx4 ../../../../bin/xx4

relink: xx4

install:
	cp xx4 ../../../../bin
