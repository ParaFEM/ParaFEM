# MAKEFILE src/modules/serial/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles common "shared" files and specific "serial" files
# Builds as a STATIC library libParaFEM_serial.a
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

SERIAL=	mp_interface.o \
	input.o
SHARED=	precision.o \
	global_variables.o \
	partition.o \
	steering.o \
	elements.o \
	timing.o \
	plasticity.o \
	geometry.o \
	fluid.o

all:	release

release:
	@echo
	@echo "Building MODULE RELEASE (SERIAL)"
	$(MAKE) FFLAGS="$(LIB_SERIAL_REL_FFLAGS)" lib-release
	@echo "Done MODULE RELEASE (SERIAL)"
debug:
	@echo
	@echo "Building MODULE DEBUG (SERIAL)"
	$(MAKE) FFLAGS="$(LIB_SERIAL_DBG_FFLAGS)" lib-debug
	@echo "Done MODULE DEBUG (SERIAL)"

lib-release:	$(SHARED) $(SERIAL)
	ar -r libParaFEM_serial.$(VERSION).a $(SHARED) $(SERIAL)

lib-debug:	$(SHARED) $(SERIAL)
	ar -r libParaFEM_serial_D.$(VERSION).a $(SHARED) $(SERIAL)


.f90.o:
	$(FC) -c $(FFLAGS) ../shared/precision.f90
	$(FC) -c $(FFLAGS) ../shared/global_variables.f90
	$(FC) -c $(FFLAGS) ../shared/partition.f90
	$(FC) -c $(FFLAGS) ../shared/steering.f90
	$(FC) -c $(FFLAGS) ../shared/elements.f90
	$(FC) -c $(FFLAGS) ../shared/timing.f90
	$(FC) -c $(FFLAGS) ../shared/plasticity.f90
	$(FC) -c $(FFLAGS) ../shared/geometry.f90
	$(FC) -c $(FFLAGS) ../shared/fluid.f90
	$(FC) -c $(FFLAGS) mp_interface.f90
	$(FC) -c $(FFLAGS) input.f90

clean:
	rm -f *.o *.mod *.a

execlean:
	rm -f ../../../include/serial/*.mod
	rm -f ../../../lib/libParaFEM_serial*.a

install:
	cp *.mod ../../../include/serial
	cp *.a ../../../lib
