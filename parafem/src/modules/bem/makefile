# MAKEFILE src/modules/bem/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles common "shared" files and specific "bem" files
# Builds as a STATIC library libParaFEM_bem.a for "release"
# Builds as a STATIC library libParaFEM_bem_D.a for "debug"
# -----------------------------------------------------------------------
# Options:
# (default) : Compile and build library
# install: Copy modules into global include directory
#          Copy library into global lib directory
# -----------------------------------------------------------------------
# Notes:
# (1) VPATH required to find "shared" modules when building lib rule
# -----------------------------------------------------------------------

VPATH = ../shared
.SUFFIXES: .f90 .o

include ../../../build/mk_defs.inc
include ../../../build/$(MACHINE).inc

BEM=	bem_lib.o
SHARED=	precision.o

all:	release

release:
	@echo
	@echo "Building MODULE RELEASE (BEM)"
	$(MAKE) FFLAGS="$(LIB_BEM_REL_FFLAGS)" lib-release
	@echo "Done MODULE RELEASE (BEM)"

debug:
	@echo
	@echo "Building MODULE DEBUG (BEM)"
	$(MAKE) FFLAGS="$(LIB_BEM_DBG_FFLAGS)" lib-debug
	@echo "Done MODULE DEBUG (BEM)"

lib-release:	$(SHARED) $(BEM)
	ar -r libParaFEM_bem.$(VERSION).a $(SHARED) $(BEM)

lib-debug:	$(SHARED) $(BEM)
	ar -r libParaFEM_bem_D.$(VERSION).a $(SHARED) $(BEM)

.f90.o:
	$(FC) -c $(FFLAGS) ../shared/precision.f90
	$(FC) -c $(FFLAGS) bem_lib.f90

clean:
	rm -f *.o *.mod *.a

install:
	cp *.mod ../../../include/bem
	cp *.a ../../../lib
