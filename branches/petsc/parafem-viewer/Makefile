# Top level makefile for ParaFEM-Viewer 

# update version number for runtime releases
VER=0.5

all:
	base -gen_proc express -exit
	$(MAKE) -f express.mk

# runtime and package rules supplied by $(MACHINE).mak
# NOTE: Windows (nmake) requires !
!include $(MACHINE).mak
