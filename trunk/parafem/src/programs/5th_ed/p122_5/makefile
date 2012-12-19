# MAKEFILE src/programs/5th_ed/p122_5/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p122_5 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p122_5 executable as RELEASE edition
#   debug: Build p122_5 executable as DEBUG edition
#   install: Copy p122_5 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p122_5 executable
#   relink: Rebuild p122_5 from object files and libraries
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
	@echo "Building P122_5 RELEASE"
	$(MAKE) p122_5 \
	FFLAGS="-c -I$(MODULE) $(P122_5_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P122_5_MPI_REL_LDFLAGS)"
	@echo "Done p122_5 RELEASE"

debug:
	@echo
	@echo "Building p122_5 DEBUG"
	$(MAKE) p122_5 \
	FFLAGS="-c -I$(MODULE) $(P122_5_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P122_5_MPI_DBG_LDFLAGS)"
	@echo "Done p122_5 DEBUG"

p122_5: p122_5.o
	$(FC) p122_5.o -o p122_5 $(LDFLAGS)

p122_5.o: p122_5.f90
	$(FC) p122_5.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p122_5 ../../../../bin/p122_5

relink: p122_5

install:
	cp p122_5 ../../../../bin
