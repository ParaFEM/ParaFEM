# MAKEFILE src/tools/postprocessing/endian_convert/makefile
# Authors: Louise M. Lever
#          Llion Evans
# -----------------------------------------------------------------------
# Compiles ttbttr program executable as RELEASE or DEBUG
# Links with RELEASE STATIC library libParaFEM_mpi.$(VERSION).a
#    or with RELEASE STATIC library libParaFEM_mpi_D.$(VERSION).a
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build upvw2dis executable as RELEASE edition
#   debug: Build upvw2dis executable as DEBUG edition
#   install: Copy upvw2dis executable into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local upvw2dis executable
#   relink: Rebuild upvw2dis from object files and libraries
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
	@echo "Building endian_convert RELEASE"
	$(MAKE) endian_convert \
	FFLAGS="-c -I$(MODULE) $(ENDIAN_CONVERT_MPI_REL_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi.$(VERSION) $(ENDIAN_CONVERT_MPI_REL_LDFLAGS)"
	@echo "Done endian_convert RELEASE"

debug:
	@echo
	@echo "Building endian_convert DEBUG"
	$(MAKE) endian_convert \
	FFLAGS="-c -I$(MODULE) $(ENDIAN_CONVERT_MPI_DBG_FFLAGS)" \
	LDFLAGS="-L$(LIB) -lParaFEM_mpi_D.$(VERSION) $(ENDIAN_CONVERT_MPI_DBG_LDFLAGS)"
	@echo "Done endian_convert DEBUG"

endian_convert: endian_convert.o
	$(FC) endian_convert.o -o endian_convert $(LDFLAGS)

endian_convert.o: endian_convert.f90
	$(FC) endian_convert.f90 $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f endian_convert ../../../../bin/endian_convert

relink: endian_convert

install:
	cp endian_convert ../../../../bin
