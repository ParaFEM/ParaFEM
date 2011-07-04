# MAKEFILE src/programs/4th_ed/p1210/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p1210 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p1210 executable as RELEASE edition
#   debug: Build p1210 executable as DEBUG edition
#   install: Copy p1210 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p1210 executable
#   relink: Rebuild p1210 from object files and libraries
# -----------------------------------------------------------------------
# Notes:
# -----------------------------------------------------------------------

include ../../../../include/mk_defs.inc
include ../../../../include/$(MACHINE).inc

MODULE	= ../../../../include
LIB	= ../../../../lib
.SUFFIXES: .o .f90

all:	release

release:
	@echo
	@echo "Building p1210 RELEASE"
	$(MAKE) p1210 \
	FFLAGS="-c -I$(MODULE) $(P1210_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P1210_MPI_REL_LDFLAGS)"
	@echo "Done p1210 RELEASE"

debug:
	@echo
	@echo "Building p1210 DEBUG"
	$(MAKE) p1210 \
	FFLAGS="-c -I$(MODULE) $(P1210_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P1210_MPI_DBG_LDFLAGS)"
	@echo "Done p1210 DEBUG"

p1210: p1210.o
	$(FC) p1210.o -o p1210 $(LDFLAGS)

p1210.o: p1210.f90
	$(FC) p1210.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p1210 ../../../../bin/p1210

relink: p1210

install:
	cp p1210 ../../../../bin
