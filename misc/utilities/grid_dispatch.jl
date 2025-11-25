#!/usr/bin/env -S julia --color=yes --startup-file=no

# This script is used to dispatch multiple workers of `grid_worker.jl`
# No grid-configuration is done in this file - all done in the worker script.
# You must run this script with the same number of threads as there are workers.

using Distributed
using Base.Threads
using Logging

@info "Started dispatch script"

# Path to script
const EXEC_WORKER::String = joinpath(dirname(@__FILE__),"grid_worker.jl")
if !isfile(EXEC_WORKER)
    @error "Could not locate grid_worker script at $grid_worker"
    exit(1)
end

# Path to julia
const EXEC_JULIA::String = joinpath(Sys.BINDIR,"julia")

# Path to AGNI root folder
const AGNI_DIR::String = abspath(dirname(@__FILE__),"..","..")

# Get number of workers specified in grid_worker
grid_script = readchomp(EXEC_WORKER)
num_work::Int = parse(Int, split(readuntil(EXEC_WORKER,"# %NUM_WORKERS"),"=")[end])
@info "Identified $num_work workers"

# Check we have the right number of threads
if nthreads() < num_work
    @error "Only $(nthreads()) threads available, but $num_work workers allocated"
    exit(1)
end

# Loop over workers to dispatch them
@threads for id in 1:num_work
    # sleep if worker>1
    if id > 1
        sleep(5.0 * id)
    end

    # Start worker...
    @info "Dispatching worker $id in thread $id..."
    run(
        pipeline(`$EXEC_JULIA --project=$AGNI_DIR $EXEC_WORKER $id`;
                    stdout=devnull
                ),
        wait=true
    )
end

@info "All workers completed!"
