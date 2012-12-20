# MAKEFILE src/programs/5th_ed/p129_5/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p129_5 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p129_5 executable as RELEASE edition
#   debug: Build p129_5 executable as DEBUG edition
#   install: Copy p129_5 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p129_5 executable
#   relink: Rebuild p129_5 from object files and libraries
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
	@echo "Building P129_5 RELEASE"
	$(MAKE) p129_5 \
	FFLAGS="-c -I$(MODULE) $(P129_5_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P129_5_MPI_REL_LDFLAGS)"
	@echo "Done p129_5 RELEASE"

debug:
	@echo
	@echo "Building p129_5 DEBUG"
	$(MAKE) p129_5 \
	FFLAGS="-c -I$(MODULE) $(P129_5_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P129_5_MPI_DBG_LDFLAGS)"
	@echo "Done p129_5 DEBUG"

p129_5: p129_5.o
	$(FC) p129_5.o -o p129_5 $(LDFLAGS)

p129_5.o: p129_5.f90
	$(FC) p129_5.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p129_5 ../../../../bin/p129_5

relink: p129_5

install:
	cp p129_5 ../../../../bin
