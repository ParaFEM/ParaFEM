# MAKEFILE src/programs/4th_ed/p125_4/makefile
# Author: Louise M. Lever, Lee Margetts
# -----------------------------------------------------------------------
# Compiles p125_4 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p125_4 executable as RELEASE edition
#   debug: Build p125_4 executable as DEBUG edition
#   install: Copy p125_4 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p125_4 executable
#   relink: Rebuild p125_4 from object files and libraries
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
	@echo "Building p125_4 RELEASE"
	$(MAKE) p125_4 \
	FFLAGS="-c -I$(MODULE) $(P125_4_4_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P125_4_4_MPI_REL_LDFLAGS)"
	@echo "Done p125_4 RELEASE"

debug:
	@echo
	@echo "Building p125_4 DEBUG"
	$(MAKE) p125_4 \
	FFLAGS="-c -I$(MODULE) $(P125_4_4_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P125_4_4_MPI_DBG_LDFLAGS)"
	@echo "Done p125_4 DEBUG"

p125_4: p125_4.o
	$(FC) p125_4.o -o p125_4 $(LDFLAGS)

p125_4.o: p125_4.f90
	$(FC) p125_4.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p125_4 ../../../../bin/p125_4

relink: p125_4

install:
	cp p125_4 ../../../../bin
