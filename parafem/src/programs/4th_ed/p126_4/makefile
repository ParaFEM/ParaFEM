# MAKEFILE src/programs/4th_ed/p126_4/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p126_4 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p126_4 executable as RELEASE edition
#   debug: Build p126_4 executable as DEBUG edition
#   install: Copy p126_4 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p126_4 executable
#   relink: Rebuild p126_4 from object files and libraries
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
	@echo "Building p126_4 RELEASE"
	$(MAKE) p126_4 \
	FFLAGS="-c -I$(MODULE) $(P126_4_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P126_4_MPI_REL_LDFLAGS)"
	@echo "Done p126_4 RELEASE"

debug:
	@echo
	@echo "Building p126_4 DEBUG"
	$(MAKE) p126_4 \
	FFLAGS="-c -I$(MODULE) $(P126_4_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P126_4_MPI_DBG_LDFLAGS)"
	@echo "Done p126_4 DEBUG"

p126_4: p126_4.o
	$(FC) p126_4.o -o p126_4 $(LDFLAGS)

p126_4.o: p126_4.f90
	$(FC) p126_4.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p126_4 ../../../../bin/p126_4

relink: p126_4

install:
	cp p126_4 ../../../../bin
