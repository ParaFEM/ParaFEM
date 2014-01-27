# MAKEFILE src/programs/5th_ed/p127/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p127 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p127 executable as RELEASE edition
#   debug: Build p127 executable as DEBUG edition
#   install: Copy p127 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p127 executable
#   relink: Rebuild p127 from object files and libraries
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
	@echo "Building P127 RELEASE"
	$(MAKE) p127 \
	FFLAGS="-c -I$(MODULE) $(P127_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P127_MPI_REL_LDFLAGS)"
	@echo "Done p127 RELEASE"

debug:
	@echo
	@echo "Building p127 DEBUG"
	$(MAKE) p127 \
	FFLAGS="-c -I$(MODULE) $(P127_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P127_MPI_DBG_LDFLAGS)"
	@echo "Done p127 DEBUG"

p127: p127.o
	$(FC) p127.o -o p127 $(LDFLAGS)

p127.o: p127.f90
	$(FC) p127.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p127 ../../../../bin/p127

relink: p127

install:
	cp p127 ../../../../bin
