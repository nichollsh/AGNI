#!/usr/bin/env -S julia --color=yes --startup-file=no

# Consolidate worker output data into single data files
# First argument must provide path to grid output folder

using Glob
using CSV
using DataFrames
using NCDatasets
using Printf

# Get output dir from ARGS
output_dir::String = ""
if length(ARGS) == 1
    output_dir = ARGS[1]
else
    println(stderr,"Invalid arguments: $ARGS")
    exit(1)
end
if isdir(output_dir)
    println("Found target directory: $output_dir")
else
    println(stderr,"Not a directory: $output_dir")
    exit(1)
end

work_dirs = glob("wk_*", output_dir)
if isempty(work_dirs)
    println(stderr,"No worker folders found in $indir")
    exit(1)
end
num_work = length(work_dirs)
println("Found $num_work worker sub-directories")

# Combined results dataframe
println("Combining result tables...")
dfs_table = DataFrame[] # not sorted
for iwk in 1:num_work
    dir = joinpath(output_dir, "wk_$iwk")
    print("    $iwk/$num_work ")
    f  = joinpath(dir, "result_table.csv")
    if isfile(f)
        println("    reading $f")
        df = CSV.read(f, DataFrame; normalizenames=true, missingstring=["", "NA", "NaN"])
        push!(dfs_table, df)
    else
        println("    skipped; could not find $f")
    end
end
combined = reduce((a,b)->vcat(a,b; cols=:union), dfs_table)
outpath = joinpath(output_dir, "consolidated_table.csv")
rm(outpath, force=true)
CSV.write(outpath, combined)
println("    wrote $(nrow(combined))x$(ncol(combined)) table to '$outpath'")
println(" ")

println("Statistics...")
num_tot = nrow(combined)
num_suc = nrow( filter(row -> row["succ"]  >  0.1, combined) )
num_fai = nrow( filter(row -> row["succ"]  < -0.1, combined) )
num_vis = num_suc + num_fai
@printf("    total:   %7d \n",num_tot)
@printf("    visited: %7d \n",num_vis )
@printf("    success: %7d (%.1f%% of visited)\n",num_suc,100*num_suc/num_vis )
@printf("    failure: %7d (%.1f%% of visited)\n",num_fai,100*num_fai/num_vis )
println(" ")

# Combined fluxes dataframe
println("Combining emission fluxes...")
dfs_emits = DataFrame[] # not sorted
for iwk in 1:num_work
    dir = joinpath(output_dir, "wk_$iwk")
    print("    $iwk/$num_work ")
    f  = joinpath(dir, "result_emits.csv")
    if isfile(f)
        println("    reading $f")
        df = CSV.read(f, DataFrame; normalizenames=true, missingstring=["", "NA", "NaN"])
        push!(dfs_emits, df)
    else
        println("    skipped; could not find $f")
    end
end
combined = reduce((a,b)->vcat(a,b; cols=:union), dfs_emits)

# manually read column names as wavelengths
head_first = split(readchomp(joinpath(work_dirs[1], "result_emits.csv")),"\n")[1]
head_emits = String["index"]
append!(head_emits, split(head_first,",")[2:end])

outpath = joinpath(output_dir, "consolidated_emits.csv")
rm(outpath, force=true)
CSV.write(outpath, combined, header=head_emits)
println("    wrote $(nrow(combined))x$(ncol(combined)) emits to '$outpath'")
println(" ")

# Combined NetCDF file
println("Combining NetCDF profiles...")
dfs_profs = Dict{String,Array}[] # not sorted
for iwk in 1:num_work
    dir = joinpath(output_dir, "wk_$iwk")
    print("    $iwk/$num_work ")
    f  = joinpath(dir, "result_profs.nc")

    if isfile(f)
        println("    reading $f")
        ds = Dataset(f) # open
        df = Dict([(k,ds[k][:,:]) for k in ("t","p","r")]) # read T,P,R arrays
        close(ds) # close
        push!(dfs_profs, df)
    else
        println("    skipped; could not find $f")
    end
end

outpath = joinpath(output_dir, "consolidated_profs.nc")
rm(outpath, force=true)
ds = Dataset(outpath,"c")

ds.attrib["description"]  = "AGNI grid consolidated profiles (t-p-r)"
ds.attrib["hostname"]     = gethostname()
ds.attrib["username"]     = ENV["USER"]

(nlev, gridsize) = size(dfs_profs[1]["t"])
defDim(ds, "nlev_c",   nlev)
defDim(ds, "gridsize", gridsize)
println("    found nlev_c = $nlev")

var_p = defVar(ds, "p", Float64, ("nlev_c","gridsize",) ) # saved in python dimension order
var_t = defVar(ds, "t", Float64, ("nlev_c","gridsize",) )
var_r = defVar(ds, "r", Float64, ("nlev_c","gridsize",) )

for iwk in 1:num_work  # for each worker directory
    # get mask of indicies that this worker performed
    mask_work = dfs_table[iwk][!,"index"][:]

    print("    $iwk/$num_work ")
    println("    writing profile")

    # loop over these indices
    for i in mask_work
        for j in 1:nlev
            # println("Get profile value for index=$i at level=$j")
            var_p[j,i] = dfs_profs[iwk]["p"][j,i]
            var_t[j,i] = dfs_profs[iwk]["t"][j,i]
            var_r[j,i] = dfs_profs[iwk]["r"][j,i]
        end
    end
end
close(ds)

println("    wrote $(gridsize)x$(nlev) t-p-r profs to '$outpath'")
println(" ")


println("Done!")
