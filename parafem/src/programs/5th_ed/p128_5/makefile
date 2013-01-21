# MAKEFILE src/programs/5th_ed/p128_5/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p128_5 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p128_5 executable as RELEASE edition
#   debug: Build p128_5 executable as DEBUG edition
#   install: Copy p128_5 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p128_5 executable
#   relink: Rebuild p128_5 from object files and libraries
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
	@echo "Building P128_5 RELEASE"
	$(MAKE) p128_5 \
	FFLAGS="-c -I$(MODULE) $(P128_5_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) -larpack_$(MACHINE) $(P128_5_MPI_REL_LDFLAGS)"
	@echo "Done p128_5 RELEASE"

debug:
	@echo
	@echo "Building p128_5 DEBUG"
	$(MAKE) p128_5 \
	FFLAGS="-c -I$(MODULE) $(P128_5_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) -larpact_$(MACHINE) $(P128_5_MPI_DBG_LDFLAGS)"
	@echo "Done p128_5 DEBUG"

p128_5: p128_5.o
	$(FC) p128_5.o -o p128_5 $(LDFLAGS)

p128_5.o: p128_5.f90
	$(FC) p128_5.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p128_5 ../../../../bin/p128_5

relink: p128_5

install:
	cp p128_5 ../../../../bin
