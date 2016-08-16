# MAKEFILE src/modules/choose_solvers/makefile
# Author: Mark Filipiak
# -----------------------------------------------------------------------
# Compiles the choose_solvers module
# Builds as a STATIC library libParaFEM_choose_solvers.a
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

CHOOSE_SOLVERS=	choose_solvers.o

all:	release

release:
	@echo
	@echo "Building MODULE RELEASE (CHOOSE_SOLVERS)"
	$(MAKE) FFLAGS="-I$(MODULE_MPI) $(LIB_CHOOSE_SOLVERS_REL_FFLAGS)" lib-release
	@echo "Done MODULE RELEASE (CHOOSE_SOLVERS)"
debug:
	@echo
	@echo "Building MODULE DEBUG (CHOOSE_SOLVERS)"
	$(MAKE) FFLAGS="-I$(MODULE_MPI) $(LIB_CHOOSE_SOLVERS_DBG_FFLAGS)" lib-debug
	@echo "Done MODULE DEBUG (CHOOSE_SOLVERS)"

lib-release:	$(CHOOSE_SOLVERS)
	ar -r libParaFEM_choose_solvers.$(VERSION).a $(CHOOSE_SOLVERS)

lib-debug:	$(CHOOSE_SOLVERS)
	ar -r libParaFEM_choose_solvers_D.$(VERSION).a $(CHOOSE_SOLVERS)


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
choose_solvers.o: choose_solvers.F90 \
                  $(MPI)/global_variables.o
	$(FC) -c $(FFLAGS) $<
$(MPI)/global_variables.o: ;

clean:
	rm -f *.o *.mod *.a

execlean:
	rm -f ../../../include/choose_solvers/*.mod
	rm -f ../../../lib/libParaFEM_choose_solvers*.a

install:
	mkdir -p ../../../include/choose_solvers
	cp *.mod ../../../include/choose_solvers
	cp *.a ../../../lib
