# MAKEFILE src/programs/5th_ed/p122/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p122 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p122 executable as RELEASE edition
#   debug: Build p122 executable as DEBUG edition
#   install: Copy p122 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p122 executable
#   relink: Rebuild p122 from object files and libraries
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
	@echo "Building P122 RELEASE"
	$(MAKE) p122 \
	FFLAGS="-c -I$(MODULE) $(P122_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P122_MPI_REL_LDFLAGS)"
	@echo "Done p122 RELEASE"

debug:
	@echo
	@echo "Building p122 DEBUG"
	$(MAKE) p122 \
	FFLAGS="-c -I$(MODULE) $(P122_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P122_MPI_DBG_LDFLAGS)"
	@echo "Done p122 DEBUG"

p122: p122.o
	$(FC) p122.o -o p122 $(LDFLAGS)

p122.o: p122.f90
	$(FC) p122.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p122 ../../../../bin/p122

relink: p122

install:
	cp p122 ../../../../bin
