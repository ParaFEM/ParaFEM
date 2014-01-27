# MAKEFILE src/programs/5th_ed/p123/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p123 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p123 executable as RELEASE edition
#   debug: Build p123 executable as DEBUG edition
#   install: Copy p123 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p123 executable
#   relink: Rebuild p123 from object files and libraries
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
	@echo "Building P123 RELEASE"
	$(MAKE) p123 \
	FFLAGS="-c -I$(MODULE) $(P123_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P123_MPI_REL_LDFLAGS)"
	@echo "Done p123 RELEASE"

debug:
	@echo
	@echo "Building p123 DEBUG"
	$(MAKE) p123 \
	FFLAGS="-c -I$(MODULE) $(P123_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P123_MPI_DBG_LDFLAGS)"
	@echo "Done p123 DEBUG"

p123: p123.o
	$(FC) p123.o -o p123 $(LDFLAGS)

p123.o: p123.f90
	$(FC) p123.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p123 ../../../../bin/p123

relink: p123

install:
	cp p123 ../../../../bin
