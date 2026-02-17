#!/bin/bash

# CHANGE
#SBATCH --time=05-02:00:00
#SBATCH --cpus-per-task=80
#SBATCH --mem-per-cpu=2750M

# LEAVE
#SBATCH --export=ALL
#SBATCH --job-name=AGNIgrid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=slurm-%x-%j-stdout.log
#SBATCH --error=slurm-%x-%j-stderr.log


echo "Home dir:          $HOME"
source $HOME/.bashrc

echo "Julia:             $(which julia)"
echo "TMPDIR:            $TMPDIR"
echo "Allocated workers: $SLURM_CPUS_PER_TASK"
echo "Allocated node:    $SLURM_JOB_NODELIST"

echo "Started at:        $(date)"
echo "Expected end:      $(date -d @$SLURM_JOB_END_TIME)    [EPOCH=$SLURM_JOB_END_TIME]"

echo " "
srun julia -t$SLURM_CPUS_PER_TASK misc/grid/manager.jl
echo " "

echo "Done"
