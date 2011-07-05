# MAKEFILE src/tools/preprocessing/mg2d/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles mg2d program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build mg2d executable as RELEASE edition
#   debug: Build mg2d executable as DEBUG edition
#   install: Copy mg2d executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local mg2d executable
#   relink: Rebuild mg2d from object files and libraries
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
	@echo "Building mg2d RELEASE"
	$(MAKE) mg2d \
	FFLAGS="-c -I$(MODULE) $(MG2D_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(MG2D_MPI_REL_LDFLAGS)"
	@echo "Done mg2d RELEASE"

debug:
	@echo
	@echo "Building mg2d DEBUG"
	$(MAKE) mg2d \
	FFLAGS="-c -I$(MODULE) $(MG2D_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(MG2D_MPI_DBG_LDFLAGS)"
	@echo "Done mg2d DEBUG"

mg2d: mg2d.o
	$(FC) mg2d.o -o mg2d $(LDFLAGS)

mg2d.o: mg2d.f90
	$(FC) mg2d.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f mg2d ../../../../bin/mg2d

relink: mg2d

install:
	cp mg2d ../../../../bin
