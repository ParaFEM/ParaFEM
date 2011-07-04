# MAKEFILE src/programs/bsd/prog82/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles prog82 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build prog82 executable as RELEASE edition
#   debug: Build prog82 executable as DEBUG edition
#   install: Copy prog82 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local prog82 executable
#   relink: Rebuild prog82 from object files and libraries
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
	@echo "Building prog82 RELEASE"
	$(MAKE) prog82 \
	FFLAGS="-c -I$(MODULE) $(PROG82_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_bem.$(VERSION) $(PROG82_MPI_REL_LDFLAGS)"
	@echo "Done prog82 RELEASE"

debug:
	@echo
	@echo "Building prog82 DEBUG"
	$(MAKE) prog82 \
	FFLAGS="-c -I$(MODULE) $(PROG82_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_bem_D.$(VERSION) $(PROG82_MPI_DBG_LDFLAGS)"
	@echo "Done prog82 DEBUG"

prog82: prog82.o
	$(FC) prog82.o -o prog82 $(LDFLAGS)

prog82.o: prog82.f90
	$(FC) prog82.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f prog82 ../../../../bin/prog82

relink: prog82

install:
	cp prog82 ../../../../bin
