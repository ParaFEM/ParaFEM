# MAKEFILE src/programs/5th_ed/p127_5/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p127_5 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p127_5 executable as RELEASE edition
#   debug: Build p127_5 executable as DEBUG edition
#   install: Copy p127_5 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p127_5 executable
#   relink: Rebuild p127_5 from object files and libraries
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
	@echo "Building P127_5 RELEASE"
	$(MAKE) p127_5 \
	FFLAGS="-c -I$(MODULE) $(P127_5_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P127_5_MPI_REL_LDFLAGS)"
	@echo "Done p127_5 RELEASE"

debug:
	@echo
	@echo "Building p127_5 DEBUG"
	$(MAKE) p127_5 \
	FFLAGS="-c -I$(MODULE) $(P127_5_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P127_5_MPI_DBG_LDFLAGS)"
	@echo "Done p127_5 DEBUG"

p127_5: p127_5.o
	$(FC) p127_5.o -o p127_5 $(LDFLAGS)

p127_5.o: p127_5.f90
	$(FC) p127_5.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p127_5 ../../../../bin/p127_5

relink: p127_5

install:
	cp p127_5 ../../../../bin
