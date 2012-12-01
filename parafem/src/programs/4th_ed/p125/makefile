# MAKEFILE src/programs/4th_ed/p125/makefile
# Author: Louise M. Lever, Lee Margetts
# -----------------------------------------------------------------------
# Compiles p125 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p125 executable as RELEASE edition
#   debug: Build p125 executable as DEBUG edition
#   install: Copy p125 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p125 executable
#   relink: Rebuild p125 from object files and libraries
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
	@echo "Building p125 RELEASE"
	$(MAKE) p125 \
	FFLAGS="-c -I$(MODULE) $(P125_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P125_MPI_REL_LDFLAGS)"
	@echo "Done p125 RELEASE"

debug:
	@echo
	@echo "Building p125 DEBUG"
	$(MAKE) p125 \
	FFLAGS="-c -I$(MODULE) $(P125_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P125_MPI_DBG_LDFLAGS)"
	@echo "Done p125 DEBUG"

p125: p125.o
	$(FC) p125.o -o p125 $(LDFLAGS)

p125.o: p125.f90
	$(FC) p125.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p125 ../../../../bin/p125

relink: p125

install:
	cp p125 ../../../../bin
