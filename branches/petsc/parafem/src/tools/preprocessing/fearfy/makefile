# MAKEFILE src/tools/preprocessing/fearfy/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles fearfy program executables.
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build fearfy executables as RELEASE edition
#   debug: Build fearfy executables as DEBUG edition
#   install: Copy fearfy executables into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local fearfy executables
#   relink: Rebuild fearfy from object files and libraries
# -----------------------------------------------------------------------
# Notes:
# -----------------------------------------------------------------------

include ../../../../build/mk_defs.inc
include ../../../../build/$(MACHINE).inc

MODULE	= ../../../../include/mpi
LIB	= ../../../../lib

.SUFFIXES: .o .f90

all:	release

release:
	@echo
	@echo "Building fearfy RELEASE"
	$(MAKE) fearfy \
	FFLAGS="-c -I$(MODULE) $(FEARFY_REL_FFLAGS)"
	@echo "Done fearfy RELEASE"

debug:
	@echo
	@echo "Building fearfy DEBUG"
	$(MAKE) fearfy \
	FFLAGS="-c -I$(MODULE) $(FEARFY_DBG_FFLAGS)"
	@echo "Done fearfy DEBUG"

fearfy: fearfy.o
	$(FC) fearfy.o -o fearfy $(LDFLAGS)

fearfy.o: fearfy.f90
	$(FC) fearfy.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f fearfy ../../../../bin/fearfy

relink: fearfy 

install:
	cp fearfy ../../../../bin
