# MAKEFILE src/programs/4th_ed/p124/makefile
# Author: Louise M. Lever, Lee Margetts
# -----------------------------------------------------------------------
# Compiles p124 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p124 executable as RELEASE edition
#   debug: Build p124 executable as DEBUG edition
#   install: Copy p124 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p124 executable
#   relink: Rebuild p124 from object files and libraries
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
	@echo "Building p124 RELEASE"
	$(MAKE) p124 \
	FFLAGS="-c -I$(MODULE) $(P124_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P124_MPI_REL_LDFLAGS)"
	@echo "Done p124 RELEASE"

debug:
	@echo
	@echo "Building p124 DEBUG"
	$(MAKE) p124 \
	FFLAGS="-c -I$(MODULE) $(P124_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P124_MPI_DBG_LDFLAGS)"
	@echo "Done p124 DEBUG"

p124: p124.o
	$(FC) p124.o -o p124 $(LDFLAGS)

p124.o: p124.f90
	$(FC) p124.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p124 ../../../../bin/p124

relink: p124

install:
	cp p124 ../../../../bin
