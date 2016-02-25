# MAKEFILE src/tools/postprocessing/pf2ensi/makefile
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
	@echo "Installing TOOL RELEASE (PF2ENSI)"
	cp pf2ensi.geo.awk ../../../../bin
	cp pf2ensi.var.awk ../../../../bin
	cp pf2ensi ../../../../bin
	chmod +x ../../../../bin/pf2ensi
	@echo "Done TOOL RELEASE (PF2ENSI)"

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
	@rm -rf ../../../../bin/pf2ensi.geo.awk
	@rm -rf ../../../../bin/pf2ensi.var.awk
	@rm -rf ../../../../bin/pf2ensi
