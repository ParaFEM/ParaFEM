# MAKEFILE src/programs/5th_ed/p128/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p128 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p128 executable as RELEASE edition
#   debug: Build p128 executable as DEBUG edition
#   install: Copy p128 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p128 executable
#   relink: Rebuild p128 from object files and libraries
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
	@echo "Building P128 RELEASE"
	$(MAKE) p128 \
	FFLAGS="-c -I$(MODULE) $(P128_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P128_MPI_REL_LDFLAGS)"
	@echo "Done p128 RELEASE"

debug:
	@echo
	@echo "Building p128 DEBUG"
	$(MAKE) p128 \
	FFLAGS="-c -I$(MODULE) $(P128_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P128_MPI_DBG_LDFLAGS)"
	@echo "Done p128 DEBUG"

p128: p128.o
	$(FC) p128.o -o p128 $(LDFLAGS)

p128.o: p128.f90
	$(FC) p128.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p128 ../../../../bin/p128

relink: p128

install:
	cp p128 ../../../../bin
