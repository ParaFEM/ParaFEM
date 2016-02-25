# MAKEFILE src/tools/postprocessing/ndttrget/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles ttbttr program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build upvw2dis executable as RELEASE edition
#   debug: Build upvw2dis executable as DEBUG edition
#   install: Copy upvw2dis executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local upvw2dis executable
#   relink: Rebuild upvw2dis from object files and libraries
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
	@echo "Building ndttrget RELEASE"
	$(MAKE) ndttrget \
	FFLAGS="-c -I$(MODULE) $(NDTTRGET_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(NDTTRGET_MPI_REL_LDFLAGS)"
	@echo "Done ndttrget RELEASE"

debug:
	@echo
	@echo "Building ndttrget DEBUG"
	$(MAKE) ndttrget \
	FFLAGS="-c -I$(MODULE) $(NDTTRGET_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(NDTTRGET_MPI_DBG_LDFLAGS)"
	@echo "Done ndttrget DEBUG"

ndttrget: ndttrget.o
	$(FC) ndttrget.o -o ndttrget $(LDFLAGS)

ndttrget.o: ndttrget.f90
	$(FC) ndttrget.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f ndttrget ../../../../bin/ndttrget

relink: ndttrget

install:
	cp ndttrget ../../../../bin
