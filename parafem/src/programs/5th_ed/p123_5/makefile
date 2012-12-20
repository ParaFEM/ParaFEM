# MAKEFILE src/programs/5th_ed/p123_5/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p123_5 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p123_5 executable as RELEASE edition
#   debug: Build p123_5 executable as DEBUG edition
#   install: Copy p123_5 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p123_5 executable
#   relink: Rebuild p123_5 from object files and libraries
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
	@echo "Building P123_5 RELEASE"
	$(MAKE) p123_5 \
	FFLAGS="-c -I$(MODULE) $(P123_5_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P123_5_MPI_REL_LDFLAGS)"
	@echo "Done p123_5 RELEASE"

debug:
	@echo
	@echo "Building p123_5 DEBUG"
	$(MAKE) p123_5 \
	FFLAGS="-c -I$(MODULE) $(P123_5_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P123_5_MPI_DBG_LDFLAGS)"
	@echo "Done p123_5 DEBUG"

p123_5: p123_5.o
	$(FC) p123_5.o -o p123_5 $(LDFLAGS)

p123_5.o: p123_5.f90
	$(FC) p123_5.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p123_5 ../../../../bin/p123_5

relink: p123_5

install:
	cp p123_5 ../../../../bin
