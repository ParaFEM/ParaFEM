# MAKEFILE src/programs/4th_ed/p129_4/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p129_4 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p129_4 executable as RELEASE edition
#   debug: Build p129_4 executable as DEBUG edition
#   install: Copy p129_4 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p129_4 executable
#   relink: Rebuild p129_4 from object files and libraries
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
	@echo "Building p129_4 RELEASE"
	$(MAKE) p129_4 \
	FFLAGS="-c -I$(MODULE) $(P129_4_4_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P129_4_4_MPI_REL_LDFLAGS)"
	@echo "Done p129_4 RELEASE"

debug:
	@echo
	@echo "Building p129_4 DEBUG"
	$(MAKE) p129_4 \
	FFLAGS="-c -I$(MODULE) $(P129_4_4_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P129_4_4_MPI_DBG_LDFLAGS)"
	@echo "Done p129_4 DEBUG"

p129_4: p129_4.o
	$(FC) p129_4.o -o p129_4 $(LDFLAGS)

p129_4.o: p129_4.f90
	$(FC) p129_4.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p129_4 ../../../../bin/p129_4

relink: p129_4

install:
	cp p129_4 ../../../../bin
