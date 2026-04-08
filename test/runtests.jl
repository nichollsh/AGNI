#!/usr/bin/env -S julia --color=yes --startup-file=no

# Run this function from inside the AGNI root folder
# e.g. as `julia --project=. test/runtests.jl`

# To run with coverage reporting, add `--code-coverage` flag

# Get AGNI root directory
ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))

# Activate environment
ENV["GKSwstype"] = "100"
import Pkg
Pkg.activate(ROOT_DIR)

# Include libraries
using LoggingExtras
using NCDatasets
using AGNI
using Glob
using Test

@info "Begin AGNI tests"

# Configure
SLOW_TESTS = ["integration", "chemistry", "deep_heating", "kzz", "spectrum"]

# Prepare
RES_DIR         = joinpath(ROOT_DIR,"res/")
OUT_DIR         = joinpath(ROOT_DIR,"out/")
TEST_DIR        = joinpath(ROOT_DIR,"test/")
p_top           = 1e-8
nlev_centre     = 100
radius          = 1.0e7    # metres
gravity         = 10.0      # m s-2
total  = 0
failed = 0

# which test suite to run?
suite::String = "all"
if length(ARGS)>0
    suite = strip(ARGS[1])
    if suite == "0"
        suite = "fast"
    end
end
@info "Requested test suite '$suite'"

rm(OUT_DIR,force=true,recursive=true)
if !isdir(OUT_DIR) && !isfile(OUT_DIR)
    mkdir(OUT_DIR)
end

rtol   = 1e-3

# Test module imported
@test isdefined(AGNI.atmosphere, :setup!)

# Find test names
test_names = sort([replace(split(basename(f), ".jl")[1], "test_"=>"")
                        for f in glob("test_*.jl", TEST_DIR)])

# Select tests
if suite == "none"
    # no tests
    @info "No tests selected. Exiting."
    exit()

elseif suite == "fast"
    # exclude slow tests
    for slow in SLOW_TESTS
        filter!(t -> !occursin(slow, t), test_names)
    end

elseif suite in test_names
    # run only the specified test suite
    test_names = [suite]

else
    if suite != "all"
        @error "Test suite '$suite' not found"
        exit(1)
    end
end

# Collect tests
test_files = String[]
for test_name in test_names
    push!(test_files, joinpath(TEST_DIR, "test_$test_name.jl"))
end
@info "Collected tests: $(join(test_names, ", "))"

# Configure logging to show only error messages during tests
LoggingExtras.global_logger(Logging.SimpleLogger(Logging.Error))

# Run tests
for test_file in test_files
    @info "Running '$(test_file)'"
    include(test_file)
end

exit(0)
