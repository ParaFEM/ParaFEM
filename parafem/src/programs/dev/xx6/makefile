# MAKEFILE src/programs/5th_ed/xx6/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles xx6 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx6 executable as RELEASE edition
#   debug: Build xx6 executable as DEBUG edition
#   install: Copy xx6 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx6 executable
#   relink: Rebuild xx6 from object files and libraries
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
	@echo "Building xx6 RELEASE"
	$(MAKE) xx6 \
	FFLAGS="-c -I$(MODULE) $(XX6_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) -larpack_$(MACHINE) $(XX6_MPI_REL_LDFLAGS)"
	@echo "Done xx6 RELEASE"

debug:
	@echo
	@echo "Building xx6 DEBUG"
	$(MAKE) xx6 \
	FFLAGS="-c -I$(MODULE) $(XX6_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) -larpack_$(MACHINE) $(XX6_MPI_DBG_LDFLAGS)"
	@echo "Done xx6 DEBUG"

xx6: xx6.o
	$(FC) xx6.o -o xx6 $(LDFLAGS)

xx6.o: xx6.f90
	$(FC) xx6.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx6 ../../../../bin/xx6

relink: xx6

install:
	cp xx6 ../../../../bin
