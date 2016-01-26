# MAKEFILE src/programs/dev/xx18/makefile
# Author: Louise M. Lever, Lee Margetts
# -----------------------------------------------------------------------
# Compiles xx18 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx18 executable as RELEASE edition
#   debug: Build xx18 executable as DEBUG edition
#   install: Copy xx18 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx18 executable
#   relink: Rebuild xx18 from object files and libraries
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
	@echo "Building XX18 RELEASE"
	$(MAKE) xx18 \
	FFLAGS="-c -I$(MODULE) $(XX18_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX18_MPI_REL_LDFLAGS)"
	@echo "Done xx18 RELEASE"

debug:
	@echo
	@echo "Building xx18 DEBUG"
	$(MAKE) xx18 \
	FFLAGS="-c -I$(MODULE) $(XX18_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX18_MPI_DBG_LDFLAGS)"
	@echo "Done xx18 DEBUG"

xx18: xx18.o
	$(FC) xx18.o -o xx18 $(LDFLAGS)

xx18.o: xx18.f90
	$(FC) xx18.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx18 ../../../../bin/xx18

relink: xx18

install:
	cp xx18 ../../../../bin

test:
	./test.sh xx18
	
