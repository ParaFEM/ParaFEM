# MAKEFILE src/tools/preprocessing/dinomesh/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles dinomesh program executables.
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build dinomesh executables as RELEASE edition
#   debug: Build dinomesh executables as DEBUG edition
#   install: Copy dinomesh executables into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local dinomesh executables
#   relink: Rebuild dinomesh from object files and libraries
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
	@echo "Building dinomesh RELEASE"
	$(MAKE) dinomesh \
	FFLAGS="-c -I$(MODULE) $(DINOMESH_REL_FFLAGS)"
	@echo "Done dinomesh RELEASE"

debug:
	@echo
	@echo "Building dinomesh DEBUG"
	$(MAKE) dinomesh \
	FFLAGS="-c -I$(MODULE) $(DINOMESH_DBG_FFLAGS)"
	@echo "Done dinomesh DEBUG"

dinomesh: dinomesh.o
	$(FC) dinomesh.o -o dinomesh $(LDFLAGS)

dinomesh.o: dinomesh.f90
	$(FC) dinomesh.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f dinomesh ../../../../bin/dinomesh

relink: dinomesh 

install:
	cp dinomesh ../../../../bin
