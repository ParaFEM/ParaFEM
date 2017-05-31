# MAKEFILE src/programs/5th_ed/xx12_b/makefile
# Author: Louise M. Lever
#         Llion M. Evans
# -----------------------------------------------------------------------
# Compiles xx12_b program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx12_b executable as RELEASE edition
#   debug: Build xx12_b executable as DEBUG edition
#   install: Copy xx12_b executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx12_b executable
#   relink: Rebuild xx12_b from object files and libraries
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
	@echo "Building xx12_b RELEASE"
	$(MAKE) xx12_b \
	FFLAGS="-c -I$(MODULE) $(XX12_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX12_MPI_REL_LDFLAGS)"
	@echo "Done xx12_b RELEASE"

debug:
	@echo
	@echo "Building xx12_b DEBUG"
	$(MAKE) xx12_b \
	FFLAGS="-c -I$(MODULE) $(XX12_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX12_MPI_DBG_LDFLAGS)"
	@echo "Done xx12_b DEBUG"

xx12_b: xx12_b.o
	$(FC) xx12_b.o -o xx12_b $(LDFLAGS)

xx12_b.o: xx12_b.f90
	$(FC) xx12_b.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx12_b ../../../../bin/xx12_b

relink: xx12_b

install:
	cp xx12_b ../../../../bin
