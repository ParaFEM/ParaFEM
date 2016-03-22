#!/bin/bash

version="0.1"

# Don't use the -debug flag.  That would make the comparison between
# ParaFEM and PETSc unfair.

if [[ $# == 1 && $1 == 'xx' ]]; then
    export build=xx log=xx.log
elif [[ $# == 1 && $1 == 'modules' ]]; then
    export build=modules log=modules.log
else
    echo 'Usage:
build.bash xx to build xx15 and xx18 only
build.bash modules to re-build the ParaFEM modules and p12meshgen'
    exit 2 # incorrect usage
fi

mkdir -p $HOME/modulefiles/parafem
export modulefile=$HOME/modulefiles/parafem/$version

#module load cray-petsc
module load cray-petsc-64

# There is clock skew between TDS and /home so sleeps are needed.

(
    export PARAFEM_HOME=$(readlink --canonicalize $PWD)
    if [[ $build == 'xx' ]]; then
	sleep 1 # TDS and /home are sleepy
	./make-parafem MACHINE=xc30 --no-libs --no-tools -mpi -parafem_petsc --install-mpi-parafem_petsc-only -xx > $log 2>&1
	sleep 1 # TDS and /home are sleepy
	cat > $modulefile <<EOF
#%Module
#
# Module parafem
#

set PARAFEM_DIR $PARAFEM_HOME
setenv PARAFEM_DIR \$PARAFEM_DIR

proc ModulesHelp { } {
    puts stderr "xx15 and xx18: ParaFEM with PETSc"
    puts stderr {=================================

Maintained by: Mark Filipiak, EPCC
}
}
prepend-path PATH \$PARAFEM_DIR/bin

module-whatis "xx15 and xx18: ParaFEM with PETSc"

EOF
	sleep 1 # TDS and /home are sleepy
    elif [[ $build == 'modules' ]]; then
	sleep 1 # TDS and /home are sleepy
	./make-parafem MACHINE=xc30 clean execlean &> clean.log
	rm -f $modulefile
	sleep 1 # TDS and /home are sleepy
	./make-parafem MACHINE=xc30 --no-libs --no-tools -mpi -parafem_petsc --install-mpi-parafem_petsc-only --only-modules > $log 2>&1
	sleep 1 # TDS and /home are sleepy
	./make-parafem MACHINE=xc30 --no-libs -mpi -parafem_petsc --install-mpi-parafem_petsc-only --only-tools -preproc >> $log 2>&1
	sleep 1 # TDS and /home are sleepy
    fi
)
