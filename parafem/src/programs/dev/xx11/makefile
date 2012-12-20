# MAKEFILE src/programs/5th_ed/xx11/makefile
# Author: Louise M. Lever
#         Llion M. Evans
# -----------------------------------------------------------------------
# Compiles xx11 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx11 executable as RELEASE edition
#   debug: Build xx11 executable as DEBUG edition
#   install: Copy xx11 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx11 executable
#   relink: Rebuild xx11 from object files and libraries
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
	@echo "Building xx11 RELEASE"
	$(MAKE) xx11 \
	FFLAGS="-c -I$(MODULE) $(XX11_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX11_MPI_REL_LDFLAGS)"
	@echo "Done xx11 RELEASE"

debug:
	@echo
	@echo "Building xx11 DEBUG"
	$(MAKE) xx11 \
	FFLAGS="-c -I$(MODULE) $(XX11_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX11_MPI_DBG_LDFLAGS)"
	@echo "Done xx11 DEBUG"

xx11: xx11.o
	$(FC) xx11.o -o xx11 $(LDFLAGS)

xx11.o: xx11.f90
	$(FC) xx11.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx11 ../../../../bin/xx11

relink: xx11

install:
	cp xx11 ../../../../bin
