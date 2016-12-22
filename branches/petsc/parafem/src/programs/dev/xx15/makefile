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
LIB	    = ../../../../lib
EXECUTABLES = xx15_elastic xx15_quadric_linear_hardening xx15_vm_nonlinear_hardening normdiff
.SUFFIXES: .o .F90 .f90

all:	release

release:
	@echo
	@echo "Building xx15 RELEASE"
	$(MAKE) $(EXECUTABLES) \
	FFLAGS="-c -I$(MODULE_MPI) -I$(MODULE_CHOOSE_SOLVERS) -I$(MODULE_PARAFEM_PETSC) $(XX15_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) $(XX15_MPI_REL_LDFLAGS)" \
	LIBS="-lParaFEM_mpi.$(VERSION) -lParaFEM_choose_solvers.$(VERSION) -lParaFEM_parafem_petsc.$(VERSION) $(PETSC_L)"
	@echo "Done xx15 RELEASE"

debug:
	@echo
	@echo "Building xx15 DEBUG"
	$(MAKE) $(EXECUTABLES) \
	FFLAGS="-c -I$(MODULE_MPI) -I$(MODULE_CHOOSE_SOLVERS) -I$(MODULE_PARAFEM_PETSC) $(XX15_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) $(XX15_MPI_DBG_LDFLAGS)" \
	LIBS="-lParaFEM_mpi_D.$(VERSION) -lParaFEM_choose_solvers_D.$(VERSION) -lParaFEM_parafem_petsc_D.$(VERSION) $(PETSC_L)"
	@echo "Done xx15 DEBUG"

xx15_elastic: xx15_elastic.o plasticity_xx15.o
	$(FC) $(LDFLAGS) -o $@ $+ $(LIBS)

xx15_quadric_linear_hardening: xx15_quadric_linear_hardening.o plasticity_xx15.o
	$(FC) $(LDFLAGS) -o $@ $+ $(LIBS)

xx15_vm_nonlinear_hardening: xx15_vm_nonlinear_hardening.o plasticity_xx15.o
	$(FC) $(LDFLAGS) -o $@ $+ $(LIBS)

xx15_original: xx15_original.o plasticity_xx15.o
	$(FC) $(LDFLAGS) -o $@ $+ $(LIBS)

normdiff: normdiff.o
	$(FC) $(LDFLAGS) -o $@ $+ $(LIBS)

# Compiling the xx15 objects needs the plasticity_xx15 module already
# compiled.
xx15_%.o : xx15_%.F90 plasticity_xx15.o
# The \"\" are to include "" around the text when it is substituted in
# the program.  The "" are to expand REVISION correctly even if it has
# special characters.  The $$ are to escape the make $ to become a
# shell $.
	$(FC) $(FFLAGS) -D REVISION=\""$$REVISION"\" -o $@ $< 

plasticity_xx15.o : plasticity_xx15.F90
	$(FC) $(FFLAGS) -o $@ $< 

normdiff.o : normdiff.F90
	$(FC) $(FFLAGS) -o $@ $< 

clean:
	rm -f *.mod *.o

execlean:
	rm -f $(EXECUTABLES)
	cd ../../../../bin; rm -f $(EXECUTABLES)

relink: $(EXECUTABLES)

install:
	cp $(EXECUTABLES) ../../../../bin

test:
#	./test.sh xx15
