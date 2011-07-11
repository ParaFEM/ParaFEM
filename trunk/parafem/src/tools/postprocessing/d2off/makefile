# MAKEFILE src/tools/postprocessing/d2off/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles d2off program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build d2off executable as RELEASE edition
#   debug: Build d2off executable as DEBUG edition
#   install: Copy d2off executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local d2off executable
#   relink: Rebuild d2off from object files and libraries
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
	@echo "Building d2off RELEASE"
	$(MAKE) d2off \
	FFLAGS="-c -I$(MODULE) $(D2OFF_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(D2OFF_MPI_REL_LDFLAGS)"
	@echo "Done d2off RELEASE"

debug:
	@echo
	@echo "Building d2off DEBUG"
	$(MAKE) d2off \
	FFLAGS="-c -I$(MODULE) $(D2OFF_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(D2OFF_MPI_DBG_LDFLAGS)"
	@echo "Done d2off DEBUG"

d2off: d2off.o
	$(FC) d2off.o -o d2off $(LDFLAGS)

d2off.o: d2off.f90
	$(FC) d2off.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f d2off ../../../../bin/d2off

relink: d2off

install:
	cp d2off ../../../../bin
