#!/bin/bash

#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=3G
#SBATCH --output=slurm-grid-%j.log

echo "Workers: $SLURM_CPUS_PER_TASK"
echo "Time: $SBATCH_TIMELIMIT"
echo "Node: $SLURM_NODEID"

julia -t$SLURM_CPUS_PER_TASK misc/grid/manager.jl

echo "Done"
