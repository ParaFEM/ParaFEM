# MAKEFILE src/programs/5th_ed/p121_5/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p121_5 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p121_5 executable as RELEASE edition
#   debug: Build p121_5 executable as DEBUG edition
#   install: Copy p121_5 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p121_5 executable
#   relink: Rebuild p121_5 from object files and libraries
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
	@echo "Building P121_5 RELEASE"
	$(MAKE) p121_5 \
	FFLAGS="-c -I$(MODULE) $(P121_5_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P121_5_MPI_REL_LDFLAGS)"
	@echo "Done p121_5 RELEASE"

debug:
	@echo
	@echo "Building p121_5 DEBUG"
	$(MAKE) p121_5 \
	FFLAGS="-c -I$(MODULE) $(P121_5_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P121_5_MPI_DBG_LDFLAGS)"
	@echo "Done p121_5 DEBUG"

p121_5: p121_5.o
	$(FC) p121_5.o -o p121_5 $(LDFLAGS)

p121_5.o: p121_5.f90
	$(FC) p121_5.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p121_5 ../../../../bin/p121_5

relink: p121_5

install:
	cp p121_5 ../../../../bin
