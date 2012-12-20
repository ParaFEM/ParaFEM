# MAKEFILE src/programs/5th_ed/p126_5/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p126_5 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p126_5 executable as RELEASE edition
#   debug: Build p126_5 executable as DEBUG edition
#   install: Copy p126_5 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p126_5 executable
#   relink: Rebuild p126_5 from object files and libraries
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
	@echo "Building P126_5 RELEASE"
	$(MAKE) p126_5 \
	FFLAGS="-c -I$(MODULE) $(P126_5_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P126_5_MPI_REL_LDFLAGS)"
	@echo "Done p126_5 RELEASE"

debug:
	@echo
	@echo "Building p126_5 DEBUG"
	$(MAKE) p126_5 \
	FFLAGS="-c -I$(MODULE) $(P126_5_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P126_5_MPI_DBG_LDFLAGS)"
	@echo "Done p126_5 DEBUG"

p126_5: p126_5.o
	$(FC) p126_5.o -o p126_5 $(LDFLAGS)

p126_5.o: p126_5.f90
	$(FC) p126_5.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p126_5 ../../../../bin/p126_5

relink: p126_5

install:
	cp p126_5 ../../../../bin
