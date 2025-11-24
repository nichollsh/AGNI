#!/usr/bin/env -S julia --color=yes --startup-file=no

# To be run from within the AGNI root directory
# Dispatch multiple workers of `grid_worker.jl`, each as a separate process
# No grid-configuration is done in this file. It's all done in `grid_worker.jl`

using Distributed
using Base.Threads
using Logging

@info "Started dispatch script"

# Path to script
const GRID_WORKER::String = joinpath(dirname(@__FILE__),"grid_worker.jl")
if !isfile(GRID_WORKER)
    @error "Could not locate grid_worker script at $grid_worker"
    exit(1)
end

# Path to julia
const JULIA::String = joinpath(Sys.BINDIR,"julia")

# Path to AGNI root folder
const ROOT_DIR::String = abspath(dirname(@__FILE__),"..","..")

# Get number of workers specified in grid_worker
grid_script = readchomp(GRID_WORKER)
num_workers::Int = parse(Int, split(readuntil(GRID_WORKER,"# NUM_WORKERS"),"=")[end])
@info "Identified $num_workers workers"

# Check we have the right number of threads
if nthreads() < num_workers
    @error "Only $(nthreads()) threads are available!"
    exit(1)
end

# Loop over workers to dispatch them
@threads for id in 1:num_workers
    @info "Dispatching worker $id in thread $id..."
    run(`$JULIA --project=$ROOT_DIR $GRID_WORKER wk_$id`, devnull, devnull, wait=true)
end
