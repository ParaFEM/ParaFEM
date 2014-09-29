# MAKEFILE src/tools/postprocessing/pf2ensibin/makefile
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
	@echo "Building pf2ensibin RELEASE"
	$(MAKE) pf2ensibin \
	FFLAGS="-c -I$(MODULE) $(PF2ENSIBIN_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(PF2ENSIBIN_MPI_REL_LDFLAGS)"
	@echo "Done pf2ensibin RELEASE"

debug:
	@echo
	@echo "Building pf2ensibin DEBUG"
	$(MAKE) pf2ensibin \
	FFLAGS="-c -I$(MODULE) $(PF2ENSIBIN_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(PF2ENSIBIN_MPI_DBG_LDFLAGS)"
	@echo "Done pf2ensibin DEBUG"

pf2ensibin: pf2ensibin.o
	$(FC) pf2ensibin.o -o pf2ensibin $(LDFLAGS)

pf2ensibin.o: pf2ensibin.f90
	$(FC) pf2ensibin.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f pf2ensibin ../../../../bin/pf2ensibin

relink: pf2ensibin

install:
	cp pf2ensibin ../../../../bin
