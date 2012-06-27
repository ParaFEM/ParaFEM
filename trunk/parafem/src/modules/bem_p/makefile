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

include ../../../build/mk_defs.inc
include ../../../build/$(MACHINE).inc

SHARED=	precision.o \
	global_variables.o \
	timing.o

MPI=	mpi_wrapper.o \
	mp_interface.o \
	input.o \
	maths.o \
	gather_scatter.o \
	bicg.o

BEM=	bem_lib_p.o

all:	release

release:
	@echo
	@echo "Building MODULE RELEASE (BEM_LIB_P)"
	$(MAKE) FFLAGS="$(LIB_BEM_P_REL_FFLAGS)" lib-release
	@echo "Done MODULE RELEASE (BEM_LIB_P)"

debug:
	@echo
	@echo "Building MODULE DEBUG (BEM_LIB_P)"
	$(MAKE) FFLAGS="$(LIB_BEM_P_DBG_FFLAGS)" lib-debug
	@echo "Done MODULE DEBUG (BEM_LIB_P)"

lib-release:	$(SHARED) $(BEM)
	ar -r libParaFEM_bem_p.$(VERSION).a $(SHARED) $(MPI) $(BEM)

lib-debug:	$(SHARED) $(BEM)
	ar -r libParaFEM_bem_p_D.$(VERSION).a $(SHARED) $(MPI) $(BEM)

.f90.o:
	$(FC) -c $(FFLAGS) ../shared/precision.f90
	$(FC) -c $(FFLAGS) ../shared/global_variables.f90
	$(FC) -c $(FFLAGS) ../shared/timing.f90
	$(FC) -c $(FFLAGS) ../mpi/mpi_wrapper.f90
	$(FC) -c $(FFLAGS) ../mpi/mp_interface.f90
	$(FC) -c $(FFLAGS) ../mpi/input.f90
	$(FC) -c $(FFLAGS) ../mpi/maths.f90
	$(FC) -c $(FFLAGS) ../mpi/gather_scatter.f90
	$(FC) -c $(FFLAGS) ../mpi/bicg.f90
	$(FC) -c $(FFLAGS) bem_lib_p.f90

clean:
	rm -f *.o *.mod *.a

execlean:
	rm -f ../../../include/bem_p/*.mod
	rm -f ../../../libParaFEM_bem_p*.a

install:
	cp *.mod ../../../include/bem_p
	cp *.a ../../../lib
