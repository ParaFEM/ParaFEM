# MAKEFILE src/programs/bsd/prog83/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles prog83 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build prog83 executable as RELEASE edition
#   debug: Build prog83 executable as DEBUG edition
#   install: Copy prog83 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local prog83 executable
#   relink: Rebuild prog83 from object files and libraries
# -----------------------------------------------------------------------
# Notes:
# -----------------------------------------------------------------------

include ../../../../include/mk_defs.inc
include ../../../../include/$(MACHINE).inc

MODULE	= ../../../../include
LIB	= ../../../../lib
.SUFFIXES: .o .f90

all:	release

release:
	@echo
	@echo "Building prog83 RELEASE"
	$(MAKE) prog83 \
	FFLAGS="-c -I$(MODULE) $(PROG83_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) -lParaFEM_bem_p.$(VERSION) $(PROG83_MPI_REL_LDFLAGS)"
	@echo "Done prog83 RELEASE"

debug:
	@echo
	@echo "Building prog83 DEBUG"
	$(MAKE) prog83 \
	FFLAGS="-c -I$(MODULE) $(PROG83_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) -lParaFEM_bem_p_D.$(VERSION) $(PROG83_MPI_DBG_LDFLAGS)"
	@echo "Done prog83 DEBUG"

prog83: prog83.o
	$(FC) prog83.o -o prog83 $(LDFLAGS)

prog83.o: prog83.f90
	$(FC) prog83.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f prog83 ../../../../bin/prog83

relink: prog83

install:
	cp prog83 ../../../../bin
