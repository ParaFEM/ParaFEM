# MAKEFILE src/programs/dev/xx19/makefile
# Author: Lee Margetts
# -----------------------------------------------------------------------
# Compiles xx19 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx19 executable as RELEASE edition
#   debug: Build xx19 executable as DEBUG edition
#   install: Copy xx19 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx19 executable
#   relink: Rebuild xx19 from object files and libraries
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
	@echo "Building xx19 RELEASE"
	$(MAKE) xx19 \
	FFLAGS="-c -I$(MODULE) $(XX19_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX19_MPI_REL_LDFLAGS)"
	@echo "Done xx19 RELEASE"

debug:
	@echo
	@echo "Building xx19 DEBUG"
	$(MAKE) xx19 \
	FFLAGS="-c -I$(MODULE) $(XX19_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX19_MPI_DBG_LDFLAGS)"
	@echo "Done xx19 DEBUG"

xx19: xx19.o
	$(FC) xx19.o -o xx19 $(LDFLAGS)

xx19.o: xx19.f90
	$(FC) xx19.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx19 ../../../../bin/xx19

relink: xx19

install:
	cp xx19 ../../../../bin
