# MAKEFILE src/tools/postprocessing/upvw2dis/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles upvw2dis program executable as RELEASE or DEBUG
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
	@echo "Building upvw2dis RELEASE"
	$(MAKE) upvw2dis \
	FFLAGS="-c -I$(MODULE) $(UPVW2DIS_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(UPVW2DIS_MPI_REL_LDFLAGS)"
	@echo "Done upvw2dis RELEASE"

debug:
	@echo
	@echo "Building upvw2dis DEBUG"
	$(MAKE) upvw2dis \
	FFLAGS="-c -I$(MODULE) $(UPVW2DIS_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(UPVW2DIS_MPI_DBG_LDFLAGS)"
	@echo "Done upvw2dis DEBUG"

upvw2dis: upvw2dis.o
	$(FC) upvw2dis.o -o upvw2dis $(LDFLAGS)

upvw2dis.o: upvw2dis.f90
	$(FC) upvw2dis.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f upvw2dis ../../../../bin/upvw2dis

relink: upvw2dis

install:
	cp upvw2dis ../../../../bin
