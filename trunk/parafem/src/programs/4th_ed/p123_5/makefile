# MAKEFILE src/programs/4th_ed/p123_4/makefile
# Author: Louise M. Lever, Lee Margetts
# -----------------------------------------------------------------------
# Compiles p123_4 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p123_4 executable as RELEASE edition
#   debug: Build p123_4 executable as DEBUG edition
#   install: Copy p123_4 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p123_4 executable
#   relink: Rebuild p123_4 from object files and libraries
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
	@echo "Building p123_4 RELEASE"
	$(MAKE) p123_4 \
	FFLAGS="-c -I$(MODULE) $(P123_4_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P123_4_MPI_REL_LDFLAGS)"
	@echo "Done p123_4 RELEASE"

debug:
	@echo
	@echo "Building p123_4 DEBUG"
	$(MAKE) p123_4 \
	FFLAGS="-c -I$(MODULE) $(P123_4_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P123_4_MPI_DBG_LDFLAGS)"
	@echo "Done p123_4 DEBUG"

p123_4: p123_4.o
	$(FC) p123_4.o -o p123_4 $(LDFLAGS)

p123_4.o: p123_4.f90
	$(FC) p123_4.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p123_4 ../../../../bin/p123_4

relink: p123_4

install:
	cp p123_4 ../../../../bin
