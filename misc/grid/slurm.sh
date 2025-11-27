#!/bin/bash

#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=3G
#SBATCH --output=slurm-grid-stdout-%j.log
#SBATCH --output=slurm-grid-stderr-%j.log

julia misc/grid/manager.jl
