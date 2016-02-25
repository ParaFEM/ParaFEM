# MAKEFILE src/tools/postprocessing/ttrb2ttr/makefile
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
	@echo "Building ttrb2ttr RELEASE"
	$(MAKE) ttrb2ttr \
	FFLAGS="-c -I$(MODULE) $(TTRB2TTR_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(TTRB2TTR_MPI_REL_LDFLAGS)"
	@echo "Done ttrb2ttr RELEASE"

debug:
	@echo
	@echo "Building ttrb2ttr DEBUG"
	$(MAKE) ttrb2ttr \
	FFLAGS="-c -I$(MODULE) $(TTRB2TTR_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(TTRB2TTR_MPI_DBG_LDFLAGS)"
	@echo "Done ttrb2ttr DEBUG"

ttrb2ttr: ttrb2ttr.o
	$(FC) ttrb2ttr.o -o ttrb2ttr $(LDFLAGS)

ttrb2ttr.o: ttrb2ttr.f90
	$(FC) ttrb2ttr.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f ttrb2ttr ../../../../bin/ttrb2ttr

relink: ttrb2ttr

install:
	cp ttrb2ttr ../../../../bin
