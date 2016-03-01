# MAKEFILE src/modules/petsc/makefile
# Author: Mark Filipiak
# -----------------------------------------------------------------------
# Compiles PETSc interface files
# Builds as a STATIC library libParaFEM_petsc.a
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

PETSC=	petsc.o

all:	release

release:
	@echo
	@echo "Building MODULE RELEASE (PETSC)"
	$(MAKE) FFLAGS="$(LIB_PETSC_REL_FFLAGS)" lib-release
	@echo "Done MODULE RELEASE (PETSC)"
debug:
	@echo
	@echo "Building MODULE DEBUG (PETSC)"
	$(MAKE) FFLAGS="$(LIB_PETSC_DBG_FFLAGS)" lib-debug
	@echo "Done MODULE DEBUG (PETSC)"

lib-release:	$(PETSC)
	ar -r libParaFEM_petsc.$(VERSION).a $(PETSC)

lib-debug:	$(PETSC)
	ar -r libParaFEM_petsc_D.$(VERSION).a $(PETSC)


.F90.o:
	$(FC) -c $(FFLAGS) petsc.F90

clean:
	rm -f *.o *.mod *.a

execlean:
	rm -f ../../../include/petsc/*.mod
	rm -f ../../../lib/libParaFEM_petsc*.a

install:
	mkdir -p ../../../include/petsc
	cp *.mod ../../../include/petsc
	cp *.a ../../../lib
