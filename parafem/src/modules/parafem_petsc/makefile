# MAKEFILE src/modules/petsc/makefile
# Author: Mark Filipiak
# -----------------------------------------------------------------------
# Compiles PETSc interface files
# Builds as a STATIC library libParaFEM_parafem_petsc.a
# -----------------------------------------------------------------------
# Options:
# (default) : Compile and build library
# install: Copy modules into global include directory
#          Copy library into global lib directory
# -----------------------------------------------------------------------

MPI = ../mpi
# The ParaFEM MPI module is not yet installed, so use the build tree
# directory.
MODULE_MPI = ../mpi

.SUFFIXES: .F90 .o
%.o: %.mod

include ../../../build/mk_defs.inc
include ../../../build/$(MACHINE).inc

PARAFEM_PETSC = parafem_petsc.o

all:	release

release:
	@echo
	@echo "Building MODULE RELEASE (PARAFEM_PETSC)"
	$(MAKE) FFLAGS="-I$(MODULE_MPI) $(LIB_PARAFEM_PETSC_REL_FFLAGS)" lib-release
	@echo "Done MODULE RELEASE (PARAFEM_PETSC)"
debug:
	@echo
	@echo "Building MODULE DEBUG (PARAFEM_PETSC)"
	$(MAKE) FFLAGS="-I$(MODULE_MPI) $(LIB_PARAFEM_PETSC_DBG_FFLAGS)" lib-debug
	@echo "Done MODULE DEBUG (PARAFEM_PETSC)"

lib-release:	$(PARAFEM_PETSC)
	ar -r libParaFEM_parafem_petsc.$(VERSION).a $(PARAFEM_PETSC)

lib-debug:	$(PARAFEM_PETSC)
	ar -r libParaFEM_parafem_petsc_D.$(VERSION).a $(PARAFEM_PETSC)


# There is a dependency on other modules.  The .mod file is not a good
# choice for the prerequisite because it can be xxx.mod or XXX.mod;
# use the .o file instead, even if this means that sometimes there is
# unneeded recompilation.  The compilation of the modules can only be
# done by the makefile in other module libraries, so don't try to
# compile it here (it would be dangerous if a parallel make was used).
# The makefile for all the modules (../makefile) contains the
# dependency between this directory and $(MPI), and that will ensure
# that the other modules (which all go into $(MPI)) are compiled
# first.
parafem_petsc.o: parafem_petsc.F90 \
                 $(MPI)/mp_interface.o \
                 $(MPI)/global_variables.o \
                 $(MPI)/gather_scatter.o
	$(FC) -c $(FFLAGS) $<
$(MPI)/mp_interface.o: ;
$(MPI)/global_variables.o: ;
$(MPI)/gather_scatter.o: ;

clean:
	rm -f *.o *.mod *.a

execlean:
	rm -f ../../../include/parafem_petsc/*.mod
	rm -f ../../../lib/libParaFEM_parafem_petsc*.a

install:
	mkdir -p ../../../include/parafem_petsc
	cp *.mod ../../../include/parafem_petsc
	cp *.a ../../../lib
