# MAKEFILE src/programs/5th_ed/p125_5/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p125_5 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p125_5 executable as RELEASE edition
#   debug: Build p125_5 executable as DEBUG edition
#   install: Copy p125_5 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p125_5 executable
#   relink: Rebuild p125_5 from object files and libraries
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
	@echo "Building P125_5 RELEASE"
	$(MAKE) p125_5 \
	FFLAGS="-c -I$(MODULE) $(P125_5_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P125_5_MPI_REL_LDFLAGS)"
	@echo "Done p125_5 RELEASE"

debug:
	@echo
	@echo "Building p125_5 DEBUG"
	$(MAKE) p125_5 \
	FFLAGS="-c -I$(MODULE) $(P125_5_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P125_5_MPI_DBG_LDFLAGS)"
	@echo "Done p125_5 DEBUG"

p125_5: p125_5.o
	$(FC) p125_5.o -o p125_5 $(LDFLAGS)

p125_5.o: p125_5.f90
	$(FC) p125_5.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p125_5 ../../../../bin/p125_5

relink: p125_5

install:
	cp p125_5 ../../../../bin
