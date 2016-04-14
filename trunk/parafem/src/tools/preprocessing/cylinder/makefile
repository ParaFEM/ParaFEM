# MAKEFILE src/tools/preprocessing/cylinder/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles cylinder program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build cylinder executable as RELEASE edition
#   debug: Build cylinder executable as DEBUG edition
#   install: Copy cylinder executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local cylinder executable
#   relink: Rebuild cylinder from object files and libraries
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
	@echo "Building cylinder RELEASE"
	$(MAKE) cylinder \
	FFLAGS="-c -I$(MODULE) $(P12MESHGEN_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P12MESHGEN_MPI_REL_LDFLAGS)"
	@echo "Done cylinder RELEASE"

debug:
	@echo
	@echo "Building cylinder DEBUG"
	$(MAKE) cylinder \
	FFLAGS="-c -I$(MODULE) $(P12MESHGEN_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P12MESHGEN_MPI_DBG_LDFLAGS)"
	@echo "Done cylinder DEBUG"

cylinder: jacobif.o mesh.o cylinder.o
	$(FC) cylinder.o jacobif.o mesh.o -o cylinder $(LDFLAGS)

jacobif.o: jacobif.f90
	$(FC) jacobif.f90 $(FFLAGS)

mesh.o: mesh.f90
	$(FC) mesh.f90 $(FFLAGS)

cylinder.o: cylinder.f90
	$(FC) cylinder.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f cylinder ../../../../bin/cylinder

relink: cylinder

install:
	cp cylinder ../../../../bin
