#!/bin/bash

# CHANGE
#SBATCH --time=00-8:00:00
#SBATCH --cpus-per-task=35
#SBATCH --mem-per-cpu=2800M

# LEAVE
#SBATCH --export=ALL
#SBATCH --job-name=AGNIgrid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=slurm-grid-%j.log


echo "Home dir:          $HOME"
source $HOME/.bashrc

echo "Julia:             $(which julia)"
echo "TMPDIR:            $TMPDIR"
echo "Allocated workers: $SLURM_CPUS_PER_TASK"
echo "Allocated node:    $SLURM_JOB_NODELIST"

echo "Started at:        $(date -d @$SLURM_JOB_START_TIME)"
echo "Expected end:      $(date -d @$SLURM_JOB_END_TIME)"

echo " "
srun julia -t$SLURM_CPUS_PER_TASK misc/grid/manager.jl
echo " "

echo "Done"
