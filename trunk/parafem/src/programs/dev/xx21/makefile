# MAKEFILE src/programs/dev/xx21/makefile
# Author: Lee Margetts
# -----------------------------------------------------------------------
# Compiles xx21 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build xx21 executable as RELEASE edition
#   debug: Build xx21 executable as DEBUG edition
#   install: Copy xx21 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local xx21 executable
#   relink: Rebuild xx21 from object files and libraries
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
	@echo "Building XX21 RELEASE"
	$(MAKE) xx21 \
	FFLAGS="-c -I$(MODULE) $(XX21_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(XX21_MPI_REL_LDFLAGS)"
	@echo "Done xx21 RELEASE"

debug:
	@echo
	@echo "Building xx21 DEBUG"
	$(MAKE) xx21 \
	FFLAGS="-c -I$(MODULE) $(XX21_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(XX21_MPI_DBG_LDFLAGS)"
	@echo "Done xx21 DEBUG"

xx21: xx21.o
	$(FC) xx21.o -o xx21 $(LDFLAGS)

xx21.o: xx21.f90
	$(FC) xx21.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f xx21 ../../../../bin/xx21

relink: xx21

install:
	cp xx21 ../../../../bin

test:
	./test.sh xx21
	
