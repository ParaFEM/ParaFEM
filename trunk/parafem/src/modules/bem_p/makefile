# MAKEFILE src/modules/bem_p/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles common "shared" files and specific "bem_p" files
# Builds as a STATIC library libParaFEM_bem_p.a for "release"
# Builds as a STATIC library libParaFEM_bem_p_D.a for "debug"
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

include ../../../include/mk_defs.inc
include ../../../include/$(MACHINE).inc

BEM=	bem_lib_p.o
SHARED=	precision.o global_variables.o

all:	release

release:
	@echo
	@echo "Building MODULE RELEASE (BEM_LIB_P)"
	$(MAKE) FFLAGS=$(LIB_BEM_P_REL_FFLAGS) lib-release
	@echo "Done MODULE RELEASE (BEM_LIB_P)"
debug:
	@echo
	@echo "Building MODULE DEBUG (BEM_LIB_P)"
	$(MAKE) FFLAGS=$(LIB_BEM_P_DBG_FFLAGS) lib-debug
	@echo "Done MODULE DEBUG (BEM_LIB_P)"

lib-release:	$(SHARED) $(BEM)
	ar -r libParaFEM_bem_p.$(VERSION).a $(SHARED) $(BEM)

lib-debug:	$(SHARED) $(BEM)
	ar -r libParaFEM_bem_p_D.$(VERSION).a $(SHARED) $(BEM)

.f90.o:
	$(FC) -c $(FFLAGS) ../shared/precision.f90
	$(FC) -c $(FFLAGS) ../shared/global_variables.f90
	$(FC) -c $(FFLAGS) bem_lib_p.f90

clean:
	rm -f *.o *.mod *.a

install:
	cp *.mod ../../../include
	cp *.a ../../../lib
