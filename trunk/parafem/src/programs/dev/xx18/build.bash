#!/bin/bash

# Don't use the -debug flag.  That would make the comparison between
# ParaFEM and PETSc unfair.

if [[ $# == 1 && $1 == 'xx18' ]]; then
    export build=xx18 log=build-xx18.log
elif [[ $# == 1 && $1 == 'modules' ]]; then
    export build=modules log=build-modules.log
else
    echo 'Usage:
build.bash xx18 to build xx18 only
build.bash modules to re-build the ParaFEM modules and p12meshgen'
    exit 2 # incorrect usage
fi

# The following two commands are needed until the new modules are set
# to be the defaults on 24 February 2016.
source $EPCC_PE_RELEASE/nov2015
# This is needed because cray-tpsl is used by cray-petsc without
# actually loading the module.
#module load cray-tpsl
module load cray-tpsl-64

#module load cray-petsc
module load cray-petsc-64

# Split in two because /work was slow.  Seems OK on /home to compile modules and xx18 all at once.

(
    export PARAFEM_HOME=$(readlink --canonicalize $PWD/../../../..)
    if [[ $build == 'xx18' ]]; then
	sleep 1 # TDS and /home are sleepy
	./make-parafem-xx18 MACHINE=xc30 --no-libs --no-tools -mpi -xx > $log 2>&1
	sleep 1 # TDS and /home are sleepy
	mkdir -p $HOME/modulefiles/xx18
	cat > $HOME/modulefiles/xx18/0.1 <<EOF
#%Module
#
# Module xx18
#

set PARAFEM_DIR $PARAFEM_HOME
setenv PARAFEM_DIR \$PARAFEM_DIR

proc ModulesHelp { } {
    puts stderr "xx18: ParaFEM p121 with PETSc"
    puts stderr {=======================

Maintained by: Mark Filipiak, EPCC
}
}
prepend-path PATH \$PARAFEM_DIR/bin

module-whatis "xx18: ParaFEM p121 with PETSc"

EOF
	sleep 1 # TDS and /home are sleepy
    elif [[ $build == 'modules' ]]; then
	./make-parafem-xx18 MACHINE=xc30 clean execlean &> clean.log
	sleep 1 # TDS and /home are sleepy
	./make-parafem-xx18 MACHINE=xc30 --no-libs --no-tools -mpi --only-modules > $log 2>&1
	sleep 1 # TDS and /home are sleepy
	./make-parafem-xx18 MACHINE=xc30 --no-libs -mpi --only-tools -preproc >> $log 2>&1
	sleep 1 # TDS and /home are sleepy
    fi
)
