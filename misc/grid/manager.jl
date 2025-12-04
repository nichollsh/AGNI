#!/usr/bin/env -S julia --color=yes --startup-file=no

# This script is used to dispatch multiple workers of `worker.jl`
#   No grid-configuration is done in this file - all done in the worker script.
#   You must run this script with the same number of threads as there are workers.

using Distributed
using Base.Threads
using Dates

println("Started manager script at $(now())")

# Path to script
const EXEC_WORKER::String = joinpath(dirname(@__FILE__),"worker.jl")
if !isfile(EXEC_WORKER)
    println(stderr,"Could not locate worker script at '$grid_worker'")
    exit(1)
end

# Path to julia
const EXEC_JULIA::String = joinpath(Sys.BINDIR,"julia")

# Path to AGNI root folder
const AGNI_DIR::String = abspath(dirname(@__FILE__),"..","..")

# Get number of workers specified in grid_worker
grid_script = readchomp(EXEC_WORKER)
num_work::Int = nthreads()
println("Requested $(nthreads()) threads, so will allocate $num_work workers")

# Loop over workers to dispatch them
@threads for id in 1:num_work
    # sleep if worker>1
    if id > 1
        sleep(10.0 + 2.0 * id)
    end

    # Start worker...
    println("Dispatching worker $id in thread $id...")
    run(
        pipeline(`$EXEC_JULIA --project=$AGNI_DIR $EXEC_WORKER $id $num_work`;
                    stdout=devnull
                ),
        wait=true
    )
end

println("All workers completed at $(now())")
println("Manager script exiting")
