# ParaFEM Top-Level Makefile
# Author: Louise M. Lever (louise.lever@manchester.ac.uk)
# Version: 0.1.3

# Set Latest ParaFEM VERSION Number
include include/mk_defs.inc

all:
	@echo
	@echo "PARAFEM BUILD - Version $(VERSION)"
	@echo
	@echo "Please specify one of the following build options:"
	@echo
	@echo "  [xt]             - Builds RELEASE library and programs for HECToR using MPI"
	@echo "  [ppc970]         - Builds RELEASE library and programs for PPC970 using MPI"
	@echo "  [sp6]            - Builds RELEASE library and programs for SP6 using MPI"
	@echo
	@echo "  [xt-bem]         - Builds RELEASE library and programs for HECToR using MPI"
	@echo "  [ppc970-bem]     - Builds RELEASE library and programs for PPC970 using MPI"
	@echo "  [sp6-bem]        - Builds RELEASE library and programs for SP6 using MPI"
	@echo
	@echo "  [xt-bem-mpi]     - Builds RELEASE library and programs for HECToR using MPI"
	@echo "  [ppc970-bem-mpi] - Builds RELEASE library and programs for PPC970 using MPI"
	@echo "  [sp6-bem-mpi]    - Builds RELEASE library and programs for SP6 using MPI"
	@echo
	@echo "  [debug-xt]             - Builds DEBUG library and programs for HECToR using MPI"
	@echo "  [debug-ppc970]         - Builds DEBUG library and programs for PPC970 using MPI"
	@echo "  [debug-sp6]            - Builds DEBUG library and programs for SP6 using MPI"
	@echo
	@echo "  [debug-xt-bem]         - Builds DEBUG library and programs for HECToR using MPI"
	@echo "  [debug-ppc970-bem]     - Builds DEBUG library and programs for PPC970 using MPI"
	@echo "  [debug-sp6-bem]        - Builds DEBUG library and programs for SP6 using MPI"
	@echo
	@echo "  [debug-xt-bem-mpi]     - Builds DEBUG library and programs for HECToR using MPI"
	@echo "  [debug-ppc970-bem-mpi] - Builds DEBUG library and programs for PPC970 using MPI"
	@echo "  [debug-sp6-bem-mpi]    - Builds DEBUG library and programs for SP6 using MPI"
	@echo
	@echo "Or select one or more of the following options:"
	@echo
	@echo "  [tools-xt]      - Builds tools for HECToR (xt) and installs to bin directory"
	@echo "  [tools-ppc970]  - Builds tools for PPC970 and installs to bin directory"
	@echo "  [tools-sp6]     - Builds tools for SP6 and installs to bin directory"
	@echo
	@echo "  [clean]    - Cleans out all module *.mod, object *.o and library lib*.a files from all directories"
	@echo "  [execlean] - Cleans out all executable programs"
	@echo

# RELEASE build rules
xt:
	$(MAKE) modules programs OPT=mpi MACHINE=xt
ppc970:
	$(MAKE) modules programs OPT=mpi MACHINE=ppc970
sp6:
	$(MAKE) modules programs OPT=mpi MACHINE=sp6
xt-bem:
	$(MAKE) modules programs OPT=bem MACHINE=xt
ppc970-bem:
	$(MAKE) modules programs OPT=bem MACHINE=ppc970
sp6-bem:
	$(MAKE) modules programs OPT=bem MACHINE=sp6
xt-bem-mpi:
	$(MAKE) modules programs OPT=bem-mpi MACHINE=xt
ppc970-bem-mpi:
	$(MAKE) modules programs OPT=bem-mpi MACHINE=ppc970
sp6-bem-mpi:
	$(MAKE) modules programs OPT=bem-mpi MACHINE=sp6

# DEBUG build rules
debug-xt:
	$(MAKE) modules-debug programs-debug OPT=mpi MACHINE=xt
debug-ppc970:
	$(MAKE) modules-debug programs-debug OPT=mpi MACHINE=ppc970
debug-sp6:
	$(MAKE) modules-debug programs-debug OPT=mpi MACHINE=sp6
debug-xt-bem:
	$(MAKE) modules-debug programs-debug OPT=bem MACHINE=xt
debug-ppc970-bem:
	$(MAKE) modules-debug programs-debug OPT=bem MACHINE=ppc970
debug-sp6-bem:
	$(MAKE) modules-debug programs-debug OPT=bem MACHINE=sp6
debug-xt-bem-mpi:
	$(MAKE) modules-debug programs-debug OPT=bem-mpi MACHINE=xt
debug-ppc970-bem-mpi:
	$(MAKE) modules-debug programs-debug OPT=bem-mpi MACHINE=ppc970
debug-sp6-bem-mpi:
	$(MAKE) modules-debug programs-debug OPT=bem-mpi MACHINE=sp6

# Module RELEASE and DEBUG rules
modules:
	@echo "Building $(OPT) Modules and Library Release"
	@cd src/modules && \
	$(MAKE) -f make$(OPT).$(MACHINE) release && \
	$(MAKE) -f make$(OPT).$(MACHINE) install && \
	cd ../..

modules-debug:
	@echo "Building $(OPT) Modules and Library Debug"
	@cd src/modules && \
	$(MAKE) -f make$(OPT).$(MACHINE) debug && \
	$(MAKE) -f make$(OPT).$(MACHINE) install && \
	cd ../..

# Programs RELEASE and DEBUG rules
programs:
	@echo "Building p121 $(OPT) Release for $(MACHINE)"
	@cd src/programs/4th_ed/p121 && \
	$(MAKE) -f makeprograms.$(MACHINE) release && \
	$(MAKE) -f makeprograms.$(MACHINE) install && \
	cd ../../../..
	@echo "Building p126 $(OPT) Release for $(MACHINE)"
	@cd src/programs/4th_ed/p126 && \
	$(MAKE) -f makeprograms.$(MACHINE) release && \
	$(MAKE) -f makeprograms.$(MACHINE) install && \
	cd ../../../..
	@echo "Building p129 $(OPT) Release for $(MACHINE)"
	@cd src/programs/4th_ed/p129 && \
	$(MAKE) -f makeprograms.$(MACHINE) release && \
	$(MAKE) -f makeprograms.$(MACHINE) install && \
	cd ../../../..
	@echo "Building p1211 $(OPT) Release for $(MACHINE)"
	@cd src/programs/5th_ed/p1211 && \
	$(MAKE) -f makeprograms.$(MACHINE) release && \
	$(MAKE) -f makeprograms.$(MACHINE) install && \
	cd ../../../..
	@echo "Building p1212 $(OPT) Release for $(MACHINE)"
	@cd src/programs/5th_ed/p1212 && \
	$(MAKE) -f makeprograms.$(MACHINE) release && \
	$(MAKE) -f makeprograms.$(MACHINE) install && \
	cd ../../../..

programs-debug:
	@echo "Building p121 $(OPT) Release for $(MACHINE)"
	@cd src/programs/4th_ed/p121 && \
	$(MAKE) -f makeprograms.$(MACHINE) debug && \
	$(MAKE) -f makeprograms.$(MACHINE) install && \
	cd ../../../..
	@echo "Building p126 $(OPT) Release for $(MACHINE)"
	@cd src/programs/4th_ed/p126 && \
	$(MAKE) -f makeprograms.$(MACHINE) debug && \
	$(MAKE) -f makeprograms.$(MACHINE) install && \
	cd ../../../..
	@echo "Building p129 $(OPT) Release for $(MACHINE)"
	@cd src/programs/4th_ed/p129 && \
	$(MAKE) -f makeprograms.$(MACHINE) debug && \
	$(MAKE) -f makeprograms.$(MACHINE) install && \
	cd ../../../..
	@echo "Building p1211 $(OPT) Release for $(MACHINE)"
	@cd src/programs/5th_ed/p1211 && \
	$(MAKE) -f makeprograms.$(MACHINE) debug && \
	$(MAKE) -f makeprograms.$(MACHINE) install && \
	cd ../../../..
	@echo "Building p1212 $(OPT) Release for $(MACHINE)"
	@cd src/programs/5th_ed/p1212 && \
	$(MAKE) -f makeprograms.$(MACHINE) debug && \
	$(MAKE) -f makeprograms.$(MACHINE) install && \
	cd ../../../..

# MACHINE specific Tool Builds and Installs
tools-xt:
	$(MAKE) tools scripts OPT=mpi MACHINE=xt

tools-ppc970:	tool-sg12mg
	$(MAKE) tools scripts OPT=mpi MACHINE=ppc970

tools-sp6:	tool-sg12mg
	$(MAKE) tools scripts OPT=mpi MACHINE=sp6

# Tool and Script Builds and Installs
tools:	tool-sg12mg

scripts:	tool-inp2d

tool-inp2d:
	@echo "Installing the inp2d script in bin"
	@cp src/tools/preprocessing/inp2d/inp2d.awk bin
	@cp src/tools/preprocessing/inp2d/inp2d bin
	@chmod +x bin/inp2d

tool-sg12mg:
	@echo "Building sg12mg Release for $(MACHINE)"
	@cd src/tools/preprocessing/ed4mesher && \
	$(MAKE) -f makeprograms.$(MACHINE) release && \
	$(MAKE) -f makeprograms.$(MACHINE) install && \
	cd ../../../..

# Clean rules
clean:
	@echo "Cleaning modules build directories"
	cd src/modules && make -f makempi.xt clean && cd ../..
	@echo "Cleaning program build directories"
	cd src/programs/4th_ed/p121 && make -f makeprograms.xt clean && cd ../../..
	cd src/programs/4th_ed/p126 && make -f makeprograms.xt clean && cd ../../..
	cd src/programs/4th_ed/p129 && make -f makeprograms.xt clean && cd ../../..
	cd src/programs/5th_ed/p1211 && make -f makeprograms.xt clean && cd ../../..
	cd src/programs/5th_ed/p1212 && make -f makeprograms.xt clean && cd ../../..
	@echo "Cleaning global include directory"
	cd include && rm -f *.mod && cd ..
	@echo "Cleaning global lib directory"
	cd lib && rm -f *.a && cd ..

execlean:
	@echo "Cleaning program build directories"
	cd src/programs/4th_ed/p121 && make -f makeprograms.xt execlean && cd ../../..
	cd src/programs/4th_ed/p126 && make -f makeprograms.xt execlean && cd ../../..
	cd src/programs/4th_ed/p129 && make -f makeprograms.xt execlean && cd ../../..
	cd src/programs/5th_ed/p1211 && make -f makeprograms.xt execlean && cd ../../..
	cd src/programs/5th_ed/p1212 && make -f makeprograms.xt execlean && cd ../../..
	@echo "Cleaning global bin directory"
	cd bin && rm -f p121 p126 p129 p1211 p1212 && cd ..
