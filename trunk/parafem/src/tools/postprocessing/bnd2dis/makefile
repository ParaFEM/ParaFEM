# MAKEFILE src/tools/postprocessing/bnd2dis/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles bnd2dis program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build bnd2dis executable as RELEASE edition
#   debug: Build bnd2dis executable as DEBUG edition
#   install: Copy bnd2dis executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local bnd2dis executable
#   relink: Rebuild bnd2dis from object files and libraries
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
	@echo "Building bnd2dis RELEASE"
	$(MAKE) bnd2dis \
	FFLAGS="-c -I$(MODULE) $(BND2DIS_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(BND2DIS_MPI_REL_LDFLAGS)"
	@echo "Done bnd2dis RELEASE"

debug:
	@echo
	@echo "Building bnd2dis DEBUG"
	$(MAKE) bnd2dis \
	FFLAGS="-c -I$(MODULE) $(BND2DIS_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(BND2DIS_MPI_DBG_LDFLAGS)"
	@echo "Done bnd2dis DEBUG"

bnd2dis: bnd2dis.o
	$(FC) bnd2dis.o -o bnd2dis $(LDFLAGS)

bnd2dis.o: bnd2dis.f90
	$(FC) bnd2dis.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f bnd2dis ../../../../bin/bnd2dis

relink: bnd2dis

install:
	cp bnd2dis ../../../../bin
