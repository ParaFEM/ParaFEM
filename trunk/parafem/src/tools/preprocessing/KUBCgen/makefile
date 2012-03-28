# MAKEFILE src/tools/preprocessing/KUBCgen/makefile
# Author: Louise M. Lever
#         Lee Margetts (modified for KUBCgen)
# -----------------------------------------------------------------------
# Compiles KUBCgen program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build KUBCgen executable as RELEASE edition
#   debug: Build KUBCgen executable as DEBUG edition
#   install: Copy KUBCgen executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local KUBCgen executable
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
	@echo "Building KUBCgen RELEASE"
	$(MAKE) KUBCgen \
	FFLAGS="-c -I$(MODULE) $(KUBCgen_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(KUBCgen_MPI_REL_LDFLAGS)"
	@echo "Done KUBCgen RELEASE"

debug:
	@echo
	@echo "Building KUBCgen DEBUG"
	$(MAKE) KUBCgen \
	FFLAGS="-c -I$(MODULE) $(KUBCgen_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(KUBCgen_MPI_DBG_LDFLAGS)"
	@echo "Done KUBCgen DEBUG"

KUBCgen: KUBCgen.o
	$(FC) KUBCgen.o -o KUBCgen $(LDFLAGS)

KUBCgen.o: KUBCgen.f90
	$(FC) KUBCgen.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f KUBCgen ../../../../bin/KUBCgen

relink: KUBCgen

install:
	cp KUBCgen ../../../../bin
