# MAKEFILE src/programs/dev/xx17/makefile
# Author: Louise M. Lever, Lee Margetts, Mark Filipiak
# -----------------------------------------------------------------------
# Compiles xx17 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with DEBUG STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx17 executable as RELEASE edition
#   debug: Build xx17 executable as DEBUG edition
#   install: Copy xx17 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx17 executable
#   relink: Rebuild xx17 from object files and libraries
# -----------------------------------------------------------------------
# Notes:
# -----------------------------------------------------------------------

include ../../../../build/mk_defs.inc
include ../../../../build/$(MACHINE).inc

MODULE_MPI            = ../../../../include/mpi
MODULE_CHOOSE_SOLVERS = ../../../../include/choose_solvers
MODULE_PARAFEM_PETSC  = ../../../../include/parafem_petsc
LIB	= ../../../../lib

.SUFFIXES: .o .F90

all:	release

release:
	@echo
	@echo "Building xx17 RELEASE"
	$(MAKE) xx17 \
	FFLAGS="-c -I$(MODULE_MPI) -I$(MODULE_CHOOSE_SOLVERS) -I$(MODULE_PARAFEM_PETSC) $(XX17_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) -lParaFEM_choose_solvers.$(VERSION) -lParaFEM_parafem_petsc.$(VERSION) $(XX17_MPI_REL_LDFLAGS)"
	@echo "Done xx17 RELEASE"

debug:
	@echo
	@echo "Building xx17 DEBUG"
	$(MAKE) xx17 \
	FFLAGS="-c -I$(MODULE_MPI) -I$(MODULE_CHOOSE_SOLVERS) -I$(MODULE_PARAFEM_PETSC) $(XX17_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) -lParaFEM_choose_solvers_D.$(VERSION) -lParaFEM_parafem_petsc_D.$(VERSION) $(XX17_MPI_DBG_LDFLAGS)"
	@echo "Done xx17 DEBUG"

xx17: xx17.o
	$(FC) -o $@ $+ $(LDFLAGS)

xx17.o: xx17.F90 makefile
# The \"\" are to include "" around the text when it is substituted in
# the program.  The "" are to expand REVISION correctly even if it has
# special characters.  The $$ are to escape the make $ to become a
# shell $.
	$(FC) $(FFLAGS) -D REVISION=\""$$REVISION"\" -o $@ $<


clean:
	rm -f *.o

execlean:
	rm -f xx17 ../../../../bin/xx17

relink: xx17

install:
	cp xx17 ../../../../bin

test:
#	./test.sh xx17
