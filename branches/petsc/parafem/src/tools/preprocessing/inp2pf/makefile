# MAKEFILE src/tools/preprocessing/inp2pf/makefile
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
	@echo "Installing TOOL RELEASE (INP2PF)"
	cp inp2pf.awk ../../../../bin
	cp inp2pf ../../../../bin
	chmod +x ../../../../bin/inp2pf
	@echo "Done TOOL RELEASE (INP2PF)"

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
	@rm -rf ../../../../bin/inp2pf.awk
	@rm -rf ../../../../bin/inp2pf
