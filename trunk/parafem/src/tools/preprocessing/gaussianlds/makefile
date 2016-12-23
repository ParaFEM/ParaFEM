# MAKEFILE src/tools/preprocessing/gaussianlds/makefile
# Author: Llion M Evans
# -----------------------------------------------------------------------
# Options:
# (default) : install
# install: Copy py script into bin directory
# -----------------------------------------------------------------------

include ../../../../build/mk_defs.inc
include ../../../../build/$(MACHINE).inc

all:	install

install:
	@echo
	@echo "Installing TOOL RELEASE (GAUSSIANLDS)"
	cp gaussianlds.py ../../../../bin
	chmod +x ../../../../bin/gaussianlds.py
	@echo "Done TOOL RELEASE (GAUSSIANLDS)"

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
	@rm -rf ../../../../bin/gaussianlds.py
