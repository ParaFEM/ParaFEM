# MAKEFILE src/programs/4th_ed/gaitfem/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles gaitfem program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_serial.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_serial_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build gaitfem executable as RELEASE edition
#   debug: Build gaitfem executable as DEBUG edition
#   install: Copy gaitfem executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local gaitfem executable
#   relink: Rebuild gaitfem from object files and libraries
# -----------------------------------------------------------------------
# Notes:
# -----------------------------------------------------------------------

include ../../../include/mk_defs.inc
include ../../../include/$(MACHINE).inc

MODULE	= ../../../include
LIB	= ../../../lib
.SUFFIXES: .o .f90

all:	release

release:
	@echo
	@echo "Building gaitfem RELEASE"
	$(MAKE) gaitfem \
	FFLAGS="-c -I$(MODULE) $(GAITFEM_SERIAL_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_serial.$(VERSION) $(GAITFEM_SERIAL_REL_LDFLAGS)"
	@echo "Done gaitfem RELEASE"

debug:
	@echo
	@echo "Building gaitfem DEBUG"
	$(MAKE) gaitfem \
	FFLAGS="-c -I$(MODULE) $(GAITFEM_SERIAL_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_serial_D.$(VERSION) $(GAITFEM_SERIAL_DBG_LDFLAGS)"
	@echo "Done gaitfem DEBUG"

gaitfem: gaitfem.o
	$(FC) gaitfem.o -o gaitfem $(LDFLAGS)

gaitfem.o: gaitfem.f90
	$(FC) gaitfem.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f gaitfem

relink: gaitfem

install:
	cp gaitfem ../../../../bin
