# MAKEFILE src/programs/dev/xx16/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles xx16 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx16 executable as RELEASE edition
#   debug: Build xx16 executable as DEBUG edition
#   install: Copy xx16 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx16 executable
#   relink: Rebuild xx16 from object files and libraries
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
	@echo "Building XX16 RELEASE"
	$(MAKE) xx16 \
	FFLAGS="-c -I$(MODULE) $(XX16_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX16_MPI_REL_LDFLAGS)"
	@echo "Done xx16 RELEASE"

debug:
	@echo
	@echo "Building xx16 DEBUG"
	$(MAKE) xx16 \
	FFLAGS="-c -I$(MODULE) $(XX16_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX16_MPI_DBG_LDFLAGS)"
	@echo "Done xx16 DEBUG"

xx16: xx16.o
	$(FC) xx16.o -o xx16 $(LDFLAGS)

xx16.o: xx16.f90
	$(FC) xx16.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx16 ../../../../bin/xx16

relink: xx16

install:
	cp xx16 ../../../../bin

test:
	./test.sh xx16
	
