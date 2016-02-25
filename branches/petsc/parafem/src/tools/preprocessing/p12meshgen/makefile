# MAKEFILE src/tools/preprocessing/p12meshgen/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p12meshgen program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p12meshgen executable as RELEASE edition
#   debug: Build p12meshgen executable as DEBUG edition
#   install: Copy p12meshgen executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p12meshgen executable
#   relink: Rebuild p12meshgen from object files and libraries
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
	@echo "Building p12meshgen RELEASE"
	$(MAKE) p12meshgen \
	FFLAGS="-c -I$(MODULE) $(P12MESHGEN_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P12MESHGEN_MPI_REL_LDFLAGS)"
	@echo "Done p12meshgen RELEASE"

debug:
	@echo
	@echo "Building p12meshgen DEBUG"
	$(MAKE) p12meshgen \
	FFLAGS="-c -I$(MODULE) $(P12MESHGEN_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P12MESHGEN_MPI_DBG_LDFLAGS)"
	@echo "Done p12meshgen DEBUG"

p12meshgen: p12meshgen.o
	$(FC) p12meshgen.o -o p12meshgen $(LDFLAGS)

p12meshgen.o: p12meshgen.f90
	$(FC) p12meshgen.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p12meshgen ../../../../bin/p12meshgen

relink: p12meshgen

install:
	cp p12meshgen ../../../../bin
