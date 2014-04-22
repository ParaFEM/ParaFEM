# MAKEFILE src/tools/preprocessing/partitioner/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles pf2metin and metout2pf program executables.
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build d2metin/metout2pf executables as RELEASE edition
#   debug: Build pf2metin/metout2pf executables as DEBUG edition
#   install: Copy pf2metin/metout2pf executables into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local pf2metin/metout2pf executables
#   relink: Rebuild pf2metin/metout2pf from object files and libraries
# -----------------------------------------------------------------------
# Notes:
# -----------------------------------------------------------------------

include ../../../../build/mk_defs.inc
include ../../../../build/$(MACHINE).inc

MODULE	= ../../../../include/mpi
LIB	= ../../../../lib

.SUFFIXES: .o .c

all:	release

release:
	@echo
	@echo "Building pf2metin and metout2pf RELEASE"
	$(MAKE) pf2metin \
	FFLAGS="-c -I$(MODULE) $(PARTITIONER_REL_FFLAGS)"
	$(MAKE) metout2pf \
	FFLAGS="-c -I$(MODULE) $(PARTITIONER_REL_FFLAGS)"
	@echo "Done pf2metin and metout2pf RELEASE"

debug:
	@echo
	@echo "Building pf2metin and metout2pf DEBUG"
	$(MAKE) pf2metin \
	FFLAGS="-c -I$(MODULE) $(PARTITIONER_DBG_FFLAGS)"
	$(MAKE) metout2pf \
	FFLAGS="-c -I$(MODULE) $(PARTITIONER_DBG_FFLAGS)"
	@echo "Done pf2metin and metout2pf DEBUG"

pf2metin: pf2metin.o
	$(CC) pf2metin.o -o pf2metin

metout2pf: metout2pf.o
	$(CC) metout2pf.o -o metout2pf

pf2metin.o: pf2metin.c
	$(CC) pf2metin.c $(FFLAGS)

metout2pf.o: metout2pf.c
	$(CC) metout2pf.c $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f pf2metin metout2pf ../../../../bin/pf2metin ../../../../bin/metout2pf

relink: pf2metin metout2pf

install:
	cp pf2metin ../../../../bin
	cp metout2pf ../../../../bin
