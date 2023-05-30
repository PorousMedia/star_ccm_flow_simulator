#!/bin/bash

## This portion of the script is to request for a computer allocation on an hpc. Slurm is like a language itself. It can get complex but below is generally what is needed.
## See https://hpcc.okstate.edu/pete-tutorial.html for basic tutorials.

#SBATCH -p express  # This is the partition
#SBATCH -t 1:00:00 # Wall time
#SBATCH --nodes=1 # Number CPU
#SBATCH --ntasks-per-node=32 # This is useful for distributed work
#SBATCH --mail-user=olubukola.ishola@okstate.edu # For getting notifications on job status
#SBATCH --mail-type=ALL # For getting notifications on job status 

## This portion of the script is the series of tasks to execute on the computer when allocated. It is the bash scripting used in linux computers.

module load star-ccm/15.02.007

scontrol show hostname $SLURM_NODELIST > machinefile_${SLURM_JOB_ID}

sim_file="Star.sim"
pod_key=" "

starccm+ -power -mpi intel -machinefile machinefile_${SLURM_JOB_ID} -rsh ssh -podkey $pod_key -licpath 1999@flex.cd-adapco.com -np $SLURM_NTASKS $sim_file -batch flow_simulator.java

rm machinefile_${SLURM_JOB_ID}