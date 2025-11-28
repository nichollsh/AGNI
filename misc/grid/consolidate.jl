#!/usr/bin/env -S julia --color=yes --startup-file=no

# Consolidate worker output data into single data files
# First argument must provide path to grid output folder

using Glob
using DelimitedFiles

# Get output dir from ARGS
output_dir::String = ""
if length(ARGS) == 1
    output_dir = ARGS[1]
else
    println(stderr,"Invalid arguments: $ARGS")
    exit(1)
end

if !isdir(output_dir)
    println(stderr,"Not a directory: $indir")
    exit(1)
end

work_dirs = glob("wk_*", output_dir)
if isempty(work_dirs)
    println(stderr,"No worker folders found in $indir")
    exit(1)
end
println("Found worker directories: $work_dirs")

# Dictionary of dictionaries - one per worker
dfs = Dict{String,Dict}()
for dir in work_dirs
    fpath = joinpath(dir,"result_table.csv")
    (data, head) = readdlm(fpath, ',', Float64; header=true)

    # Read into dictionary
    df = Dict{String,Array}()
    for (i,h) in enumerate(head)
       df[h] = data[:,i]
    end

    # Store by worker index integer
    dfs[parse(Int,split(dir,"_")[end])] = df
end

# Consolidate results into a single dictionary (each with 2D array)
results_table = Dict{String,Array}()
for wk in sort(keys(dfs))
    println("Adding worker $wk")
end


# combined = reduce((a,b)->vcat(a,b; cols=:union), dfs)
# if outpath === nothing
#     outpath = joinpath(indir, "combined.csv")
# end
# mkpath(dirname(outpath))
# CSV.write(outpath, combined)
# @info "Wrote $(nrow(combined)) rows, $(ncol(combined)) cols" outpath
