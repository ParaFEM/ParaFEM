# MAKEFILE src/programs/4th_ed/p1210_4/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p1210_4 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p1210_4 executable as RELEASE edition
#   debug: Build p1210_4 executable as DEBUG edition
#   install: Copy p1210_4 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p1210_4 executable
#   relink: Rebuild p1210_4 from object files and libraries
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
	@echo "Building p1210_4 RELEASE"
	$(MAKE) p1210_4 \
	FFLAGS="-c -I$(MODULE) $(P1210_4_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P1210_4_MPI_REL_LDFLAGS)"
	@echo "Done p1210_4 RELEASE"

debug:
	@echo
	@echo "Building p1210_4 DEBUG"
	$(MAKE) p1210_4 \
	FFLAGS="-c -I$(MODULE) $(P1210_4_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P1210_4_MPI_DBG_LDFLAGS)"
	@echo "Done p1210_4 DEBUG"

p1210_4: p1210_4.o
	$(FC) p1210_4.o -o p1210_4 $(LDFLAGS)

p1210_4.o: p1210_4.f90
	$(FC) p1210_4.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p1210_4 ../../../../bin/p1210_4

relink: p1210_4

install:
	cp p1210_4 ../../../../bin
