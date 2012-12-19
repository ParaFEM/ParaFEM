#!/bin/bash --login

#$ -N parafem_test_gpu

# Set the job time
#$ -l h_rt=1:00:0

# Set the account to charge to (change this to your account)
#$ -A ge199

# Shift to the directory that the job was submitted from
#$ -cwd

# Send environment with script (needed to get code modules)
#$ -V

# Select queue
###$ -q tesla@gpu3
#$ -q tesla

# Merge standard output and error
#$ -j y

##export CUDA_VISIBLE_DEVICES="1,2"

./xx3 ed4-gpu1-c3d20






