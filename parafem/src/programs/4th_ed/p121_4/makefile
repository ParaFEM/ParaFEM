# MAKEFILE src/programs/4th_ed/p121_4/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles p121_4 program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build p121_4 executable as RELEASE edition
#   debug: Build p121_4 executable as DEBUG edition
#   install: Copy p121_4 executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local p121_4 executable
#   relink: Rebuild p121_4 from object files and libraries
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
	@echo "Building p121_4 RELEASE"
	$(MAKE) p121_4 \
	FFLAGS="-c -I$(MODULE) $(P121_4_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(P121_4_MPI_REL_LDFLAGS)"
	@echo "Done p121_4 RELEASE"

debug:
	@echo
	@echo "Building p121_4 DEBUG"
	$(MAKE) p121_4 \
	FFLAGS="-c -I$(MODULE) $(P121_4_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(P121_4_MPI_DBG_LDFLAGS)"
	@echo "Done p121_4 DEBUG"

p121_4: p121_4.o
	$(FC) p121_4.o -o p121_4 $(LDFLAGS)

p121_4.o: p121_4.f90
	$(FC) p121_4.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f p121_4 ../../../../bin/p121_4

relink: p121_4

install:
	cp p121_4 ../../../../bin
