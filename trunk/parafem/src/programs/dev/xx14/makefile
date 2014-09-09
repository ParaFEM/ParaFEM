# MAKEFILE src/programs/5th_ed/xx14/makefile
# Author: Anton Shterenlikht, 9-SEP-2014.
# Based on p121 by Louise M. Lever
# -----------------------------------------------------------------------
# Compiles xx14 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx14 executable as RELEASE edition
#   debug: Build xx14 executable as DEBUG edition
#   install: Copy xx14 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx14 executable
#   relink: Rebuild xx14 from object files and libraries
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
	@echo "Building XX14 RELEASE"
	$(MAKE) xx14 \
	FFLAGS="-c -I$(MODULE) $(XX14_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX14_MPI_REL_LDFLAGS)"
	@echo "Done xx14 RELEASE"

debug:
	@echo
	@echo "Building xx14 DEBUG"
	$(MAKE) xx14 \
	FFLAGS="-c -I$(MODULE) $(XX14_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX14_MPI_DBG_LDFLAGS)"
	@echo "Done xx14 DEBUG"

xx14: xx14.o
	$(FC) xx14.o -o xx14 $(LDFLAGS)

xx14.o: xx14.f90
	$(FC) xx14.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx14 ../../../../bin/xx14

relink: xx14

install:
	cp xx14 ../../../../bin

test:
	./test.sh xx14
	
