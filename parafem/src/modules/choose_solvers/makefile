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

.SUFFIXES: .F90 .o
%.o: %.mod

include ../../../build/mk_defs.inc
include ../../../build/$(MACHINE).inc

CHOOSE_SOLVERS=	choose_solvers.o

all:	release

release:
	@echo
	@echo "Building MODULE RELEASE (CHOOSE_SOLVERS)"
	$(MAKE) FFLAGS="$(LIB_CHOOSE_SOLVERS_REL_FFLAGS)" lib-release
	@echo "Done MODULE RELEASE (CHOOSE_SOLVERS)"
debug:
	@echo
	@echo "Building MODULE DEBUG (CHOOSE_SOLVERS)"
	$(MAKE) FFLAGS="$(LIB_CHOOSE_SOLVERS_DBG_FFLAGS)" lib-debug
	@echo "Done MODULE DEBUG (CHOOSE_SOLVERS)"

lib-release:	$(CHOOSE_SOLVERS)
	ar -r libParaFEM_choose_solvers.$(VERSION).a $(CHOOSE_SOLVERS)

lib-debug:	$(CHOOSE_SOLVERS)
	ar -r libParaFEM_choose_solvers_D.$(VERSION).a $(CHOOSE_SOLVERS)


.F90.o:
	$(FC) -c $(FFLAGS) choose_solvers.F90

clean:
	rm -f *.o *.mod *.a

execlean:
	rm -f ../../../include/choose_solvers/*.mod
	rm -f ../../../lib/libParaFEM_choose_solvers*.a

install:
	mkdir -p ../../../include/choose_solvers
	cp *.mod ../../../include/choose_solvers
	cp *.a ../../../lib
