#!/bin/bash

## Example SLURM batch script to use 1 host cpu and 1 gpu on Bede 

## To run, type:                $ sbatch job.bede
## To check queue status, type: $ squeue -u <username>
## To cancel job, type:         $ scancel <JOBID>

#SBATCH --account=<ACCOUNTID>
#SBATCH --time=0:5:0

#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gres=gpu:1

nvidia-smi

bede-mpirun --bede-par 1ppg -np 1 ~/parafem/parafem/bin/xx3 xx3_small

