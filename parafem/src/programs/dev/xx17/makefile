# MAKEFILE src/programs/dev/xx17/makefile
# Author: Louise M. Lever, Lee Margetts
# -----------------------------------------------------------------------
# Compiles xx17 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx17 executable as RELEASE edition
#   debug: Build xx17 executable as DEBUG edition
#   install: Copy xx17 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx17 executable
#   relink: Rebuild xx17 from object files and libraries
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
	@echo "Building XX17 RELEASE"
	$(MAKE) xx17 \
	FFLAGS="-c -I$(MODULE) $(XX17_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX17_MPI_REL_LDFLAGS)"
	@echo "Done xx17 RELEASE"

debug:
	@echo
	@echo "Building xx17 DEBUG"
	$(MAKE) xx17 \
	FFLAGS="-c -I$(MODULE) $(XX17_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX17_MPI_DBG_LDFLAGS)"
	@echo "Done xx17 DEBUG"

xx17: xx17.o
	$(FC) xx17.o -o xx17 $(LDFLAGS)

xx17.o: xx17.f90
	$(FC) xx17.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx17 ../../../../bin/xx17

relink: xx17

install:
	cp xx17 ../../../../bin

test:
	./test.sh xx17
	
