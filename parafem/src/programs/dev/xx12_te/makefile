# MAKEFILE src/programs/dev/xx12_te/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles rfemsolve_te program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx12_te executable as RELEASE edition
#   debug: Build xx12_te executable as DEBUG edition
#   install: Copy xx12_te executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx12_te executable
#   relink: Rebuild xx12_te from object files and libraries
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
	@echo "Building xx12_te RELEASE"
	$(MAKE) xx12_te \
	FFLAGS="-c -I$(MODULE) $(XX12_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX12_MPI_REL_LDFLAGS)"
	@echo "Done xx12_te RELEASE"

debug:
	@echo
	@echo "Building xx12_te DEBUG"
	$(MAKE) xx12_te \
	FFLAGS="-c -I$(MODULE) $(XX12_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX12_MPI_DBG_LDFLAGS)"
	@echo "Done xx12_te DEBUG"

xx12_te: xx12_te.o
	$(FC) xx12_te.o -o xx12_te $(LDFLAGS)

xx12_te.o: xx12_te.f90
	$(FC) xx12_te.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx12_te ../../../../bin/xx12_te

relink: xx12_te

install:
	cp xx12_te ../../../../bin
