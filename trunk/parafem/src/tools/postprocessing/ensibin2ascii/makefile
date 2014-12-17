# MAKEFILE src/tools/postprocessing/ensibin2ascii/makefile
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
	@echo "Building ensibin2ascii RELEASE"
	$(MAKE) ensibin2ascii \
	FFLAGS="-c -I$(MODULE) $(ENSIBIN2ASCII_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(ENSIBIN2ASCII_MPI_REL_LDFLAGS)"
	@echo "Done ensibin2ascii RELEASE"

debug:
	@echo
	@echo "Building ensibin2ascii DEBUG"
	$(MAKE) ensibin2ascii \
	FFLAGS="-c -I$(MODULE) $(ENSIBIN2ASCII_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(ENSIBIN2ASCII_MPI_DBG_LDFLAGS)"
	@echo "Done ensibin2ascii DEBUG"

ensibin2ascii: ensibin2ascii.o
	$(FC) ensibin2ascii.o -o ensibin2ascii $(LDFLAGS)

ensibin2ascii.o: ensibin2ascii.f90
	$(FC) ensibin2ascii.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f ensibin2ascii ../../../../bin/ensibin2ascii

relink: ensibin2ascii

install:
	cp ensibin2ascii ../../../../bin
