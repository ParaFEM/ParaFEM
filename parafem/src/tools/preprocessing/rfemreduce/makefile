# MAKEFILE src/tools/preprocessing/rfemreduce/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles rfemreduce program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build rfemreduce executable as RELEASE edition
#   debug: Build rfemreduce executable as DEBUG edition
#   install: Copy rfemreduce executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local rfemreduce executable
#   relink: Rebuild rfemreduce from object files and libraries
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
	@echo "Building rfemreduce RELEASE"
	$(MAKE) rfemreduce \
	FFLAGS="-c -I$(MODULE) $(RFEMREDUCE_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(RFEMREDUCE_MPI_REL_LDFLAGS)"
	@echo "Done rfemreduce RELEASE"

debug:
	@echo
	@echo "Building rfemreduce DEBUG"
	$(MAKE) rfemreduce \
	FFLAGS="-c -I$(MODULE) $(RFEMREDUCE_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(RFEMREDUCE_MPI_DBG_LDFLAGS)"
	@echo "Done rfemreduce DEBUG"

rfemreduce: rfemreduce.o
	$(FC) rfemreduce.o -o rfemreduce $(LDFLAGS) 

rfemreduce.o: rfemreduce.f90
	$(FC) rfemreduce.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f rfemreduce ../../../../bin/rfemreduce

relink: rfemreduce

install:
	$(INSTALL_CMD) rfemreduce ../../../../bin
