# MAKEFILE src/programs/dev/xx15/makefile
# Author: Louise M. Lever, Mark Filipiak
# -----------------------------------------------------------------------
# Compiles xx15 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with DEBUG STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx15 executable as RELEASE edition
#   debug: Build xx15 executable as DEBUG edition
#   install: Copy xx15 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx15 executable
#   relink: Rebuild xx15 from object files and libraries
# -----------------------------------------------------------------------
# Notes:
# -----------------------------------------------------------------------

include ../../../../build/mk_defs.inc
include ../../../../build/$(MACHINE).inc

MODULE_MPI            = ../../../../include/mpi
MODULE_CHOOSE_SOLVERS = ../../../../include/choose_solvers
MODULE_PARAFEM_PETSC  = ../../../../include/parafem_petsc
LIB	= ../../../../lib
.SUFFIXES: .o .F90 .f90

all:	release

release:
	@echo
	@echo "Building xx15 RELEASE"
	$(MAKE) xx15 \
	FFLAGS="-c -I$(MODULE_MPI) -I$(MODULE_CHOOSE_SOLVERS) -I$(MODULE_PARAFEM_PETSC) $(XX15_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) -lParaFEM_choose_solvers.$(VERSION) -lParaFEM_parafem_petsc.$(VERSION) $(XX15_MPI_REL_LDFLAGS)"
#	$(MAKE) xx15_original \
#	FFLAGS="-c -I$(MODULE_MPI) $(XX15_MPI_REL_FFLAGS)" \
#	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX15_MPI_REL_LDFLAGS)"
	@echo "Done xx15 RELEASE"

debug:
	@echo
	@echo "Building xx15 DEBUG"
	$(MAKE) xx15 \
	FFLAGS="-c -I$(MODULE_MPI) -I$(MODULE_CHOOSE_SOLVERS) -I$(MODULE_PARAFEM_PETSC) $(XX15_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) -lParaFEM_choose_solvers_D.$(VERSION) -lParaFEM_parafem_petsc_D.$(VERSION) $(XX15_MPI_DBG_LDFLAGS)"
#	$(MAKE) xx15_original \
#	FFLAGS="-c -I$(MODULE_MPI) $(XX15_MPI_DBG_FFLAGS)" \
#	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX15_MPI_DBG_LDFLAGS)"
	@echo "Done xx15 DEBUG"

xx15: xx15.o
	$(FC) xx15.o -o xx15 $(LDFLAGS)

xx15.o: xx15.F90
	$(FC) xx15.F90 $(FFLAGS)

xx15_original: xx15_original.o
	$(FC) xx15_original.o -o xx15_original  $(LDFLAGS)

xx15_original.o: xx15_original.f90
	$(FC) xx15_original.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx15 xx15_original ../../../../bin/xx15 ../../../../bin/xx15_original

relink: xx15 #xx15_original

install:
	cp xx15 ../../../../bin
#	cp xx15_original ../../../../bin

test:
	./test.sh xx15
#	./test.sh xx15_original
