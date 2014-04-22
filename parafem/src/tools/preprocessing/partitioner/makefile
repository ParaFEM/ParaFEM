# MAKEFILE src/tools/preprocessing/partitioner/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles pf2metin and metout2pf program executables.
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build d2metin/metoutpf executables as RELEASE edition
#   debug: Build pf2metin/metoutpf executables as DEBUG edition
#   install: Copy pf2metin/metoutpf executables into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local pf2metin/metoutpf executables
#   relink: Rebuild pf2metin/metoutpf from object files and libraries
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
	@echo "Building pf2metin and metoutpf RELEASE"
	$(MAKE) pf2metin \
	FFLAGS="-c -I$(MODULE) $(PARTITIONER_REL_FFLAGS)"
	$(MAKE) metoutpf \
	FFLAGS="-c -I$(MODULE) $(PARTITIONER_REL_FFLAGS)"
	@echo "Done pf2metin and metoutpf RELEASE"

debug:
	@echo
	@echo "Building pf2metin and metoutpf DEBUG"
	$(MAKE) pf2metin \
	FFLAGS="-c -I$(MODULE) $(PARTITIONER_DBG_FFLAGS)"
	$(MAKE) metoutpf \
	FFLAGS="-c -I$(MODULE) $(PARTITIONER_DBG_FFLAGS)"
	@echo "Done pf2metin and metoutpf DEBUG"

pf2metin: pf2metin.o
	$(CC) pf2metin.o -o pf2metin

metoutpf: metoutpf.o
	$(CC) metoutpf.o -o metoutpf

pf2metin.o: pf2metin.c
	$(CC) pf2metin.c $(FFLAGS)

metoutpf.o: metoutpf.c
	$(CC) metoutpf.c $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f pf2metin metoutpf ../../../../bin/pf2metin ../../../../bin/metoutpf

relink: pf2metin metoutpf

install:
	cp pf2metin ../../../../bin
	cp metoutpf ../../../../bin
