#!/bin/sh
# inputs:
# 1 - name of the executable file, xx14.x or xx14noio.x
# 2 - queue, short or standard
project=e347
user=mexas
bld_dir=/home/$project/$project/$user/parafem/parafem/bin
work_dir=/work/$project/$project/$user
prog=$1
queue=$2

if [ $prog = "xx14.x" ]
then
	job_script=zpfem
elif [ $prog = "xx14noio.x" ]
then
	job_script=zscal
else
	echo $prog does not exist
	exit 1
fi
# check the file exists under /work
if [ -x $work_dir/$prog ]
then
	differ=`diff -q $bld_dir/$prog $work_dir/$prog`
else
	differ="$prog does not exist under /work"
fi
# if the program has been updated, or does not exist
# under /work at all, copy it there and run
if [ -z "$differ" ]
then
	echo $prog not changed
else
	echo $differ
	cd $work_dir
	rsync -qca $bld_dir/$prog . &&  qsub -q $queue $job_script && qstat -u $USER
fi
