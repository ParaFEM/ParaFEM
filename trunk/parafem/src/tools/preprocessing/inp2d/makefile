# MAKEFILE src/tools/preprocessing/inp2d/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Options:
# (default) : install
# install: Copy shell and awk script into bin directory
# -----------------------------------------------------------------------

include ../../../../build/mk_defs.inc
include ../../../../build/$(MACHINE).inc

all:	install

install:
	@echo
	@echo "Installing TOOL RELEASE (INP2D)"
	cp inp2d.awk ../../../../bin
	cp inp2d ../../../../bin
	chmod +x ../../../../bin/inp2d
	@echo "Done TOOL RELEASE (INP2D)"

release:
	@echo
	@echo "Nothing to do for this tool"

debug:
	@echo
	@echo "Nothing to do for this tool"

clean:
	@echo
	@echo "Nothing to do for this tool"

execlean:
	@rm -rf ../../../../bin/inp2d.awk
	@rm -rf ../../../../bin/inp2d
