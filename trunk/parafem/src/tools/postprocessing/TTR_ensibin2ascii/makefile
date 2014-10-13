# MAKEFILE src/tools/postprocessing/TTR_ensibin2ascii/makefile
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
	@echo "Building TTR_ensibin2ascii RELEASE"
	$(MAKE) TTR_ensibin2ascii \
	FFLAGS="-c -I$(MODULE) $(TTR_ENSIBIN2ASCII_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(TTR_ENSIBIN2ASCII_MPI_REL_LDFLAGS)"
	@echo "Done TTR_ensibin2ascii RELEASE"

debug:
	@echo
	@echo "Building TTR_ensibin2ascii DEBUG"
	$(MAKE) TTR_ensibin2ascii \
	FFLAGS="-c -I$(MODULE) $(TTR_ENSIBIN2ASCII_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(TTR_ENSIBIN2ASCII_MPI_DBG_LDFLAGS)"
	@echo "Done TTR_ensibin2ascii DEBUG"

TTR_ensibin2ascii: TTR_ensibin2ascii.o
	$(FC) TTR_ensibin2ascii.o -o TTR_ensibin2ascii $(LDFLAGS)

TTR_ensibin2ascii.o: TTR_ensibin2ascii.f90
	$(FC) TTR_ensibin2ascii.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f TTR_ensibin2ascii ../../../../bin/TTR_ensibin2ascii

relink: TTR_ensibin2ascii

install:
	cp TTR_ensibin2ascii ../../../../bin
