#!/bin/bash

#SBATCH --time=00-3:00:00
#SBATCH --cpus-per-task=32

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --output=slurm-grid-%j.log

echo "Allocated workers: $SLURM_CPUS_PER_TASK"
echo "Allocated node:    $SLURM_JOB_NODELIST"

echo "Started at:        $(date -d @$SLURM_JOB_START_TIME)"
echo "Expected end:      $(date -d @$SLURM_JOB_END_TIME)"

echo " "
julia -t$SLURM_CPUS_PER_TASK misc/grid/manager.jl
echo " "

echo "Done"
