# MAKEFILE src/programs/5th_ed/p129/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p129 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p129 executable as RELEASE edition
#   debug: Build p129 executable as DEBUG edition
#   install: Copy p129 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p129 executable
#   relink: Rebuild p129 from object files and libraries
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
	@echo "Building P129 RELEASE"
	$(MAKE) p129 \
	FFLAGS="-c -I$(MODULE) $(P129_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P129_MPI_REL_LDFLAGS)"
	@echo "Done p129 RELEASE"

debug:
	@echo
	@echo "Building p129 DEBUG"
	$(MAKE) p129 \
	FFLAGS="-c -I$(MODULE) $(P129_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P129_MPI_DBG_LDFLAGS)"
	@echo "Done p129 DEBUG"

p129: p129.o
	$(FC) p129.o -o p129 $(LDFLAGS)

p129.o: p129.f90
	$(FC) p129.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p129 ../../../../bin/p129

relink: p129

install:
	cp p129 ../../../../bin
