# MAKEFILE src/tools/preprocessing/partitioner/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles d2metin and metout2d program executables.
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Build d2metin/metout2d executables as RELEASE edition
#   debug: Build d2metin/metout2d executables as DEBUG edition
#   install: Copy d2metin/metout2d executables into global bin directory
#   clean: Remove .o object files
#   execlean: Remove the local d2metin/metout2d executables
#   relink: Rebuild d2metin/metout2d from object files and libraries
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
	@echo "Building d2metin and metout2d RELEASE"
	$(MAKE) d2metin \
	FFLAGS="-c -I$(MODULE) $(PARTITIONER_REL_FFLAGS)"
	$(MAKE) metout2d \
	FFLAGS="-c -I$(MODULE) $(PARTITIONER_REL_FFLAGS)"
	@echo "Done d2metin and metout2d RELEASE"

debug:
	@echo
	@echo "Building d2metin and metout2d DEBUG"
	$(MAKE) d2metin \
	FFLAGS="-c -I$(MODULE) $(PARTITIONER_DBG_FFLAGS)"
	$(MAKE) metout2d \
	FFLAGS="-c -I$(MODULE) $(PARTITIONER_DBG_FFLAGS)"
	@echo "Done d2metin and metout2d DEBUG"

d2metin: d2metin.o
	$(CC) d2metin.o -o d2metin

metout2d: metout2d.o
	$(CC) metout2d.o -o metout2d

d2metin.o: d2metin.c
	$(CC) d2metin.c $(FFLAGS)

metout2d.o: metout2d.c
	$(CC) metout2d.c $(FFLAGS)

clean:
	rm -f *.o

execlean:
	rm -f d2metin metout2d ../../../../bin/d2metin ../../../../bin/metout2d

relink: d2metin metout2d

install:
	cp d2metin ../../../../bin
	cp metout2d ../../../../bin
