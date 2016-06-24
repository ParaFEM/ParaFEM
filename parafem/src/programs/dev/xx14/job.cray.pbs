#!/bin/bash --login
#
#$Id$
#
#        Resource: ARCHER (Cray XC30 (24-core per node))
#    Batch system: PBSPro_select
#
#PBS -j oe
#PBS -A ecse0505
#PBS -l select=8
#PBS -l walltime=00:20:0
#PBS -N l8

module add PrgEnv-cray perftools craype-hugepages8M atp stat
module list
export XT_SYMMETRIC_HEAP_SIZE=700m
export ATP_ENABLED=1
export PGAS_MEMINFO_DISPLAY=1
#export CRAY_PGAS_USE_DMAPP_QUEUE=y
export MPICH_MPIIO_STATS=1

# resolve all symlinks to absolute paths
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)

# Switch to current working directory
cd $PBS_O_WORKDIR

# Run the parallel program
#aprun -n 192 -N 24 -S 12 -d 1 -T testABM.x
aprun -n 192 -N 24 -S 12 -d 1 -T xx14.x p121_large
#aprun -n 30000 -N 24 -S 12 -d 1 -T xx14noio-dm0.5.x+pat p121_large 
#sleep 600
#pid=`ps | grep aprun | head -1 | awk '{print $1}'`
#stat-cl -c $pid
