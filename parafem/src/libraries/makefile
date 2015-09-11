# MAKEFILE src/libraries/makefile
# Author: Louise M. Lever
# -----------------------------------------------------------------------
# Compiles external libraries required by some ParaFEM programs.
# -----------------------------------------------------------------------
# Options:
# (default/all)  : Compile and install all libraries (release)
# release-arpack : Compile and install the ARPACK library
# release-metis  : Compile and install the METIS library and tools
# -----------------------------------------------------------------------

include ../../build/mk_defs.inc
include ../../build/$(MACHINE).inc

all:	release-arpack release-metis

release-arpack:
	@echo
	@echo "Building and Installing LIBRARY RELEASE (ARPACK)"
	@cd ARPACK && $(MAKE) lib && cd ..
	@echo "Done LIBRARY RELEASE (ARPACK)"

debug-arpack:
	@echo
	@echo "****************************************************************"
	@echo "** There is currently no DEBUG option for the arpack library  **"
	@echo "** Building the RELEASE version instead                       **"
	@echo "****************************************************************"
	$(MAKE) release-arpack

release-metis:
	@echo
	@echo "Building LIBRARY RELEASE (METIS)"
	@cd metis-4.0.3 && $(MAKE) default && cd ..
	@echo "Done LIBRARY RELEASE (METIS)"

debug-metis:
	@echo
	@echo "****************************************************************"
	@echo "** There is currently no DEBUG option for the metis library   **"
	@echo "** Building the RELEASE version instead                       **"
	@echo "****************************************************************"
	$(MAKE) release-metis

release-gaf77:
	@echo
	@echo "Building LIBRARY RELEASE (GAF77)"
	@cd gaf77 && $(MAKE) FC=$(EXT_LIB_GAF77_FC) FFLAGS="$(EXT_LIB_GAF77_FFLAGS)" LDFLAGS="$(EXT_LIB_GAF77_LDFLAGS)" && cd ..
	@echo "Done LIBRARY RELEASE (GAF77)"

debug-gaf77:
	@echo
	@echo "****************************************************************"
	@echo "** There is currently no DEBUG option for the gaf77 library   **"
	@echo "** Building the RELEASE version instead                       **"
	@echo "****************************************************************"
	$(MAKE) release-gaf77

clean:
	@cd ARPACK && $(MAKE) clean && cd ..
#	@cd metis-4.0.3 && $(MAKE) clean && cd ..

execlean:
	@rm -f ../../lib/libarpack*.a
	@rm -f ../../lib/libmetis*.a
	@rm -f gaf77/*.o
	@rm -f gaf77/*.a
#	@cd metis-4.0.3 && $(MAKE) realclean && cd ..

#install:
#	@echo
#	@echo "Installing: "
#	$(MAKE) install-mpi install-serial install-bem_lib install-bem_lib_p

#install-mpi:
#	@cd mpi && $(MAKE) install & cd ..

