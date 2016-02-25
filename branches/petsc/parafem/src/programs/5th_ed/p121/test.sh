#!/bin/bash

if [ $# -ne 1 ]; then
	echo "Usage: test.sh <driver>"; exit 1
fi

do_test () {
	
	t_file="${1}"
	t_pattern="${2}"
	t_counter="${3}"
	t_value=$(grep "${t_pattern}" "${demo_file_path}/${t_file}")
	
	echo -n "** Test: ${t_counter}... "
	#echo "Looking at: ${tmp_dir}/demo/${t_file}"
	#echo "Comparing with "${demo_file_path}/${t_file}"
	#echo "Looking for: $t_pattern"
	#echo "The reference value is: $t_value"

	# The test
	grep "${t_value}" ${tmp_dir}/demo/${t_file} >> ${log_file}
	
	if [ $? -eq 0 ]; then
		
		echo "[PASSED]" | tee -a ${log_file}
		e_code=0
	else
		echo "[FAILED]" | tee -a ${log_file}
		echo "* Generated output in ${t_file} did not match reference value" | tee -a ${log_file}
        echo "* Expected: ${t_value}" | tee -a ${log_file}
        echo "* Found: $(grep "${t_pattern}" ${tmp_dir}/demo/${t_file})" | tee -a ${log_file}
        e_code=1
	fi
	return ${e_code}	 
}

# Start here
driver_name=$1
echo "*** Started testing ${driver_name}" | tee ${log_file}

# Set up logging
log_file=$(pwd)/${driver_name}-test.log
rm -f ${log_file}

# Load the details of the paths and tests to be carried out
source ./tests
demo_file_path=$(readlink -e ${demo_file_path}) # We want an absolute path

# Set up a temporary working directory
tmp_dir=$(mktemp -d)
cp -a "${demo_file_path}" "${tmp_dir}/"

# Run the executable - If this fails, we know we're in trouble...
# Execute from the top-level bin dir.
# FIXME - We should able to execute in place, but it this causes issues finding the input files

echo -n "** Executing ${driver_name}... "
cd ../../../../bin
cpu_cores=$(expr $(nproc) / 2)
mpirun --mca btl self,sm -np ${cpu_cores} ${driver_name} ${tmp_dir}/demo/${driver_name}_demo 2>&1 >> ${log_file}

if [ $? -eq 0 ]; then
	# Executed correctly, check output
	#echo "* ${driver_name} Executed OK" | tee -a ${log_file}
	echo "[PASSED]" | tee -a ${log_file}
	for ((i=0; i < ${#test_file[@]}; i++)) do
		do_test "${test_file[$i]}" "${test_pattern[$i]}" $i		
		exit_code=$(expr $exit_code + $?)
	done
else
	echo "[FAILED]" | tee -a ${log_file}
	echo "* ${driver_name} failed to execute" | tee -a ${log_file}
    exit_code=2
fi

if [ ${exit_code} -eq 0 ]; then
	echo "*** ${driver_name} All tests: [PASSED]" | tee -a ${log_file}
    rm ${log_file}
else
    echo "*** ${driver_name} At least one test: [FAILED]" | tee -a ${log_file}
    echo "*** See log for details: ${log_file}"
fi  
rm -rf "${tmp_dir}"


exit ${exit_code}

