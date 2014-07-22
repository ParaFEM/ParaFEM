# MAKEFILE src/tools/preprocessing/p12meshgenbin/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p12meshgenbin program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p12meshgenbin executable as RELEASE edition
#   debug: Build p12meshgenbin executable as DEBUG edition
#   install: Copy p12meshgenbin executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p12meshgenbin executable
#   relink: Rebuild p12meshgenbin from object files and libraries
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
	@echo "Building p12meshgenbin RELEASE"
	$(MAKE) p12meshgenbin \
	FFLAGS="-c -I$(MODULE) $(P12MESHGENBIN_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P12MESHGENBIN_MPI_REL_LDFLAGS)"
	@echo "Done p12meshgenbin RELEASE"

debug:
	@echo
	@echo "Building p12meshgenbin DEBUG"
	$(MAKE) p12meshgenbin \
	FFLAGS="-c -I$(MODULE) $(P12MESHGENBIN_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P12MESHGENBIN_MPI_DBG_LDFLAGS)"
	@echo "Done p12meshgenbin DEBUG"

p12meshgenbin: p12meshgenbin.o
	$(FC) p12meshgenbin.o -o p12meshgenbin $(LDFLAGS)

p12meshgenbin.o: p12meshgenbin.f90
	$(FC) p12meshgenbin.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p12meshgenbin ../../../../bin/p12meshgenbin

relink: p12meshgenbin

install:
	cp p12meshgenbin ../../../../bin
