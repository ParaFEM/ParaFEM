# MAKEFILE src/programs/5th_ed/p124_5/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p124_5 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p124_5 executable as RELEASE edition
#   debug: Build p124_5 executable as DEBUG edition
#   install: Copy p124_5 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p124_5 executable
#   relink: Rebuild p124_5 from object files and libraries
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
	@echo "Building P124_5 RELEASE"
	$(MAKE) p124_5 \
	FFLAGS="-c -I$(MODULE) $(P124_5_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P124_5_MPI_REL_LDFLAGS)"
	@echo "Done p124_5 RELEASE"

debug:
	@echo
	@echo "Building p124_5 DEBUG"
	$(MAKE) p124_5 \
	FFLAGS="-c -I$(MODULE) $(P124_5_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P124_5_MPI_DBG_LDFLAGS)"
	@echo "Done p124_5 DEBUG"

p124_5: p124_5.o
	$(FC) p124_5.o -o p124_5 $(LDFLAGS)

p124_5.o: p124_5.f90
	$(FC) p124_5.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p124_5 ../../../../bin/p124_5

relink: p124_5

install:
	cp p124_5 ../../../../bin
