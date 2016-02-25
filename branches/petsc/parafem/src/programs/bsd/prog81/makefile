# MAKEFILE src/programs/bsd/prog81/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles prog81 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_bem.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_bem_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build prog81 executable as RELEASE edition
#   debug: Build prog81 executable as DEBUG edition
#   install: Copy prog81 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local prog81 executable
#   relink: Rebuild prog81 from object files and libraries
# -----------------------------------------------------------------------
# Notes:
# -----------------------------------------------------------------------

include ../../../../build/mk_defs.inc
include ../../../../build/$(MACHINE).inc

MODULE	= ../../../../include/bem
LIB	= ../../../../lib
.SUFFIXES: .o .f90

all:	release

release:
	@echo
	@echo "Building prog81 RELEASE"
	$(MAKE) prog81 \
	FFLAGS="-c -I$(MODULE) $(PROG81_BEM_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_bem.$(VERSION) $(PROG81_BEM_REL_LDFLAGS)"
	@echo "Done prog81 RELEASE"

debug:
	@echo
	@echo "Building prog81 DEBUG"
	$(MAKE) prog81 \
	FFLAGS="-c -I$(MODULE) $(PROG81_BEM_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_bem_D.$(VERSION) $(PROG81_BEM_DBG_LDFLAGS)"
	@echo "Done prog81 DEBUG"

prog81: prog81.o
	$(FC) prog81.o -o prog81 $(LDFLAGS)

prog81.o: prog81.f90
	$(FC) prog81.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f prog81 ../../../../bin/prog81

relink: prog81

install:
	cp prog81 ../../../../bin
