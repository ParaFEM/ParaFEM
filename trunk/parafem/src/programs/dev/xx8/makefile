# MAKEFILE src/programs/5th_ed/xx8/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles xx8 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx8 executable as RELEASE edition
#   debug: Build xx8 executable as DEBUG edition
#   install: Copy xx8 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx8 executable
#   relink: Rebuild xx8 from object files and libraries
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
	@echo "Building XX8 RELEASE"
	$(MAKE) xx8 \
	FFLAGS="-c -I$(MODULE) $(XX8_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX8_MPI_REL_LDFLAGS)"
	@echo "Done xx8 RELEASE"

debug:
	@echo
	@echo "Building xx8 DEBUG"
	$(MAKE) xx8 \
	FFLAGS="-c -I$(MODULE) $(XX8_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX8_MPI_DBG_LDFLAGS)"
	@echo "Done xx8 DEBUG"

xx8: xx8.o
	$(FC) xx8.o -o xx8 $(LDFLAGS)

xx8.o: xx8.f90
	$(FC) xx8.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx8 ../../../../bin/xx8

relink: xx8

install:
	cp xx8 ../../../../bin
