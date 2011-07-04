# MAKEFILE src/programs/4th_ed/p126/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p126 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p126 executable as RELEASE edition
#   debug: Build p126 executable as DEBUG edition
#   install: Copy p126 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p126 executable
#   relink: Rebuild p126 from object files and libraries
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
	@echo "Building p126 RELEASE"
	$(MAKE) p126 \
	FFLAGS="-c -I$(MODULE) $(P126_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P126_MPI_REL_LDFLAGS)"
	@echo "Done p126 RELEASE"

debug:
	@echo
	@echo "Building p126 DEBUG"
	$(MAKE) p126 \
	FFLAGS="-c -I$(MODULE) $(P126_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P126_MPI_DBG_LDFLAGS)"
	@echo "Done p126 DEBUG"

p126: p126.o
	$(FC) p126.o -o p126 $(LDFLAGS)

p126.o: p126.f90
	$(FC) p126.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p126 ../../../../bin/p126

relink: p126

install:
	cp p126 ../../../../bin
