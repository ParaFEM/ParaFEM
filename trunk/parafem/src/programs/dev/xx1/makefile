# MAKEFILE src/programs/5th_ed/xx1/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles xx1 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx1 executable as RELEASE edition
#   debug: Build xx1 executable as DEBUG edition
#   install: Copy xx1 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx1 executable
#   relink: Rebuild xx1 from object files and libraries
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
	@echo "Building xx1 RELEASE"
	$(MAKE) xx1 \
	FFLAGS="-c -I$(MODULE) $(XX1_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX1_MPI_REL_LDFLAGS)"
	@echo "Done xx1 RELEASE"

debug:
	@echo
	@echo "Building xx1 DEBUG"
	$(MAKE) xx1 \
	FFLAGS="-c -I$(MODULE) $(XX1_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX1_MPI_DBG_LDFLAGS)"
	@echo "Done xx1 DEBUG"

xx1: pf-elastic.o xx1.o
	$(FC) xx1.o pf-elastic.o -o xx1 $(LDFLAGS)

pf-elastic.o:
#	$(FC) ../../../umats/parafem/pf-elastic.f $(FFLAGS)
	$ bgf77 ../../../umats/parafem/pf-elastic.f $(FFLAGS)

xx1.o: xx1.f90
	$(FC) xx1.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx1 ../../../../bin/xx1

relink: xx1

install:
	cp xx1 ../../../../bin
