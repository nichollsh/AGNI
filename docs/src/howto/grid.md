# Running a grid of models

AGNI is not explicitly parallelised. However, there is functionality to run a grid
of models using the script located at `misc/grid/worker.jl`. This script operates a
single worker (of potentially many). All parameter configuration should be done by editing
the `worker.jl` file directly.

## Running workers manually

For example, to run worker ID=1 and allocate two workers to the whole grid:
```console
julia --project=. misc/grid/worker.jl 1 2
```

## Running a managed parallel grid

By allocating multiple workers and running them simultaneously using the manager script
at `misc/grid/manager.jl`, simulations can be parallelised. The number of workers is
defined by the number of threads with which the manager script is run.

For example, run the manager with four threads (and thus four workers):
```console
julia -t4 misc/grid/manager.jl
```

## Running on an HPC cluster with Slurm

The manager script can also be executed on a HPC cluster using the Slurm workload
manager. Define the number of workers (and the allotted time) inside `slurm.sh` by
setting `cpus-per-task=XX`.

For example:
```console
sbatch --export=ALL misc/grid/slurm.sh
```

## Consolidating results

Once all workers have finished, combine results using the `consolidate.jl` script.
This generates `consolidated_*` files (2 CSV, 1 NetCDF) in the root folder of the grid.

For example, assuming the grid root is at `out/`:
```console
julia misc/grid/consolidate.jl out/
```
