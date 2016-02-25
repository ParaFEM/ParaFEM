# MAKEFILE src/tools/preprocessing/iotest/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles iotest program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build iotest executable as RELEASE edition
#   debug: Build iotest executable as DEBUG edition
#   install: Copy iotest executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local iotest executable
#   relink: Rebuild iotest from object files and libraries
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
	@echo "Building iotest RELEASE"
	$(MAKE) iotest \
	FFLAGS="-c -I$(MODULE) $(IOTEST_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(IOTEST_MPI_REL_LDFLAGS)"
	@echo "Done iotest RELEASE"

debug:
	@echo
	@echo "Building iotest DEBUG"
	$(MAKE) iotest \
	FFLAGS="-c -I$(MODULE) $(IOTEST_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(IOTEST_MPI_DBG_LDFLAGS)"
	@echo "Done iotest DEBUG"

iotest: iotest.o
	$(FC) iotest.o -o iotest $(LDFLAGS)

iotest.o: iotest.f90
	$(FC) iotest.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f iotest ../../../../bin/iotest

relink: iotest

install:
	cp iotest ../../../../bin
