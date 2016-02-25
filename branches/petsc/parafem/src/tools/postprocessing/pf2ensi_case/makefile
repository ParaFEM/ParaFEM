# MAKEFILE src/tools/postprocessing/pf2ensibin/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles pf2ensi_case program executable as RELEASE or DEBUG
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
	@echo "Building pf2ensi_case RELEASE"
	$(MAKE) pf2ensi_case \
	FFLAGS="-c -I$(MODULE) $(PF2ENSI_CASE_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(PF2ENSI_CASE_MPI_REL_LDFLAGS)"
	@echo "Done pf2ensi_case RELEASE"

debug:
	@echo
	@echo "Building pf2ensibin DEBUG"
	$(MAKE) pf2ensi_case \
	FFLAGS="-c -I$(MODULE) $(PF2ENSI_CASE_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(PF2ENSI_CASE_MPI_DBG_LDFLAGS)"
	@echo "Done pf2ensi_case DEBUG"

pf2ensi_case: pf2ensi_case.o
	$(FC) pf2ensi_case.o -o pf2ensi_case $(LDFLAGS)

pf2ensi_case.o: pf2ensi_case.f90
	$(FC) pf2ensi_case.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f pf2ensi_case ../../../../bin/pf2ensi_case

relink: pf2ensi_case

install:
	cp pf2ensi_case ../../../../bin
