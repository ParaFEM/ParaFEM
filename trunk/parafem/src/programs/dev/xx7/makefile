# MAKEFILE src/programs/5th_ed/xx7/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles xx7 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx7 executable as RELEASE edition
#   debug: Build xx7 executable as DEBUG edition
#   install: Copy xx7 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx7 executable
#   relink: Rebuild xx7 from object files and libraries
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
	@echo "Building xx7 RELEASE"
	$(MAKE) xx7 \
	FFLAGS="-c -I$(MODULE) $(XX7_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX7_MPI_REL_LDFLAGS)"
	@echo "Done xx7 RELEASE"

debug:
	@echo
	@echo "Building xx7 DEBUG"
	$(MAKE) xx7 \
	FFLAGS="-c -I$(MODULE) $(XX7_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX7_MPI_DBG_LDFLAGS)"
	@echo "Done xx7 DEBUG"

xx7: xx7.o
	$(FC) xx7.o -o xx7 $(LDFLAGS)

xx7.o: xx7.f90
	$(FC) xx7.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx7 ../../../../bin/xx7

relink: xx7

install:
	cp xx7 ../../../../bin
