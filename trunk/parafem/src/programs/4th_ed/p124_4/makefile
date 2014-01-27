# MAKEFILE src/programs/4th_ed/p124_4/makefile
# Author: Louise M. Lever, Lee Margetts
# -----------------------------------------------------------------------
# Compiles p124_4 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p124_4 executable as RELEASE edition
#   debug: Build p124_4 executable as DEBUG edition
#   install: Copy p124_4 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p124_4 executable
#   relink: Rebuild p124_4 from object files and libraries
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
	@echo "Building p124_4 RELEASE"
	$(MAKE) p124_4 \
	FFLAGS="-c -I$(MODULE) $(P124_4_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P124_4_MPI_REL_LDFLAGS)"
	@echo "Done p124_4 RELEASE"

debug:
	@echo
	@echo "Building p124_4 DEBUG"
	$(MAKE) p124_4 \
	FFLAGS="-c -I$(MODULE) $(P124_4_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P124_4_MPI_DBG_LDFLAGS)"
	@echo "Done p124_4 DEBUG"

p124_4: p124_4.o
	$(FC) p124_4.o -o p124_4 $(LDFLAGS)

p124_4.o: p124_4.f90
	$(FC) p124_4.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p124_4 ../../../../bin/p124_4

relink: p124_4

install:
	cp p124_4 ../../../../bin
