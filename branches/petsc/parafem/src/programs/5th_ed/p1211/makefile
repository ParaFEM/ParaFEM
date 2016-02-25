# MAKEFILE src/programs/5th_ed/p1211/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p1211 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p1211 executable as RELEASE edition
#   debug: Build p1211 executable as DEBUG edition
#   install: Copy p1211 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p1211 executable
#   relink: Rebuild p1211 from object files and libraries
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
	@echo "Building P1211 RELEASE"
	$(MAKE) p1211 \
	FFLAGS="-c -I$(MODULE) $(P1211_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P1211_MPI_REL_LDFLAGS)"
	@echo "Done p1211 RELEASE"

debug:
	@echo
	@echo "Building p1211 DEBUG"
	$(MAKE) p1211 \
	FFLAGS="-c -I$(MODULE) $(P1211_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P1211_MPI_DBG_LDFLAGS)"
	@echo "Done p1211 DEBUG"

p1211: p1211.o
	$(FC) p1211.o -o p1211 $(LDFLAGS)

p1211.o: p1211.f90
	$(FC) p1211.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p1211 ../../../../bin/p1211

relink: p1211

install:
	cp p1211 ../../../../bin
