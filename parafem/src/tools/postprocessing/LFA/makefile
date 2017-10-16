# MAKEFILE src/tools/postprocessing/LFA/makefile
# Author: Llion Evans
# -----------------------------------------------------------------------
# Copies LFA postprocessing programs as RELEASE or DEBUG
# -----------------------------------------------------------------------
# Available Rules:
#   (default) : execute the release rule
#
#   release: Does nothing
#   debug: Does nothing
#   install: Copy LFA scripts into global bin directory
#   execlean: Remove the local LFA script
#   relink: 
# -----------------------------------------------------------------------
# Notes:
# -----------------------------------------------------------------------

all:	release

release:
	@echo
	@echo "Done LFA RELEASE"

debug:
	@echo
	@echo "Done LFA DEBUG"

execlean:
	rm -f LFA_calc_k.py ../../../../bin/LFA_calc_k.py

relink: LFA_calc_k.py

install:
	cp LFA_calc_k.py ../../../../bin
	chmod -f 755 ../../../../bin/LFA_calc_k.py
