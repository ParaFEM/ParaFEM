#!/bin/bash

version="0.1"

# Don't use the -debug flag.  That would make the comparison between
# ParaFEM and PETSc unfair.

if [[ $# == 1 && $1 == 'xx' ]]; then
    export build=xx log=xx.log
elif [[ $# == 1 && $1 == 'modules' ]]; then
    export build=modules log=modules.log
elif [[ $# == 1 && $1 == 'clean' ]]; then
    export build=clean log=clean.log
elif [[ $# == 1 && $1 == 'execlean' ]]; then
    export build=execlean log=execlean.log
else
    echo 'Usage:
build.bash clean to clean everything
build.bash modules to re-build the ParaFEM modules and p12meshgen
build.bash xx to build xx15 and xx18 only'
    exit 2 # incorrect usage
fi

mkdir -p $HOME/modulefiles/parafem
export modulefile=$HOME/modulefiles/parafem/$version

# When the default modules are up to date, then only this line is
# needed.
#module load cray-petsc-64

# To get the 14 April 2016 Archer modules the following lines are
# needed, including an explicit load of cray-tpsl-64.  Switch back to
# a simple 'module load cray-petsc-64' when the defaults are changed.
module swap modules modules/3.2.10.3
module swap cray-mpich cray-mpich/7.3.2
module swap craype craype/2.5.3
module swap cce cce/8.4.5
module swap cray-libsci cray-libsci/16.03.1
module load cray-tpsl-64/16.03.1
module load cray-petsc-64/3.6.3.0

module list -t 2>&1 | head -n 1 > $log
module list -t 2>&1 | tail -n +2 | sort >> $log
echo >> $log

# There is clock skew between TDS and /home so sleeps are needed.

(
    export PARAFEM_HOME=$(readlink --canonicalize $PWD)
    if [[ $build == 'xx' ]]; then
	sleep 1 # TDS and /home are sleepy
	./make-parafem MACHINE=xc30 --no-libs --no-tools -mpi -parafem_petsc --install-mpi-parafem_petsc-only -xx >> $log 2>&1
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
	./make-parafem MACHINE=xc30 --no-libs --no-tools -mpi -parafem_petsc --install-mpi-parafem_petsc-only --only-modules >> $log 2>&1
	sleep 1 # TDS and /home are sleepy
	./make-parafem MACHINE=xc30 --no-libs -mpi -parafem_petsc --install-mpi-parafem_petsc-only --only-tools -preproc >> $log 2>&1
	sleep 1 # TDS and /home are sleepy
    elif [[ $build == 'clean' ]]; then
	sleep 1 # TDS and /home are sleepy
	./make-parafem MACHINE=xc30 clean execlean >> $log 2>&1
	rm -f $modulefile
	sleep 1 # TDS and /home are sleepy
    elif [[ $build == 'execlean' ]]; then
	sleep 1 # TDS and /home are sleepy
	./make-parafem MACHINE=xc30 execlean >> $log 2>&1
	rm -f $modulefile
	sleep 1 # TDS and /home are sleepy
    fi
)
