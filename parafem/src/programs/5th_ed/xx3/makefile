# MAKEFILE src/programs/5th_ed/xx3/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles xx3 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx3 executable as RELEASE edition
#   debug: Build xx3 executable as DEBUG edition
#   install: Copy xx3 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx3 executable
#   relink: Rebuild xx3 from object files and libraries
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
	@echo "Building xx3 RELEASE"
	$(MAKE) xx3 \
	FFLAGS="-c -I$(MODULE) $(XX3_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX3_MPI_REL_LDFLAGS)"
	@echo "Done xx3 RELEASE"

debug:
	@echo
	@echo "Building xx3 DEBUG"
	$(MAKE) xx3 \
	FFLAGS="-c -I$(MODULE) $(XX3_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX3_MPI_DBG_LDFLAGS)"
	@echo "Done xx3 DEBUG"

xx3: xx3.o
	$(FC) xx3.o -o xx3 $(LDFLAGS)

xx3.o: xx3.f90
	$(FC) xx3.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx3 ../../../../bin/xx3

relink: xx3

install:
	cp xx3 ../../../../bin
