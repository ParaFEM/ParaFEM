#!/bin/sh
project=e347
user=mexas
bld_dir=/home/$project/$project/$user/parafem/parafem/bin
work_dir=/work/$project/$project/$user
prog=xx14noio.x
job_script=zscal

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
	rsync -qca $bld_dir/$prog .
	qsub -q short $job_script
	qstat -u $USER
fi
