#!/usr/bin/env -S julia --color=yes --startup-file=no
# Run this function from inside the `test/` folder
# e.g. as `julia --project=.. runtests.jl `

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
using Test

@info "Begin AGNI tests"

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
@info "Using test suite '$suite'"

rm(OUT_DIR,force=true,recursive=true)
if !isdir(OUT_DIR) && !isfile(OUT_DIR)
    mkdir(OUT_DIR)
end

rtol   = 1e-3

# Test module imported
@test isdefined(AGNI.atmosphere, :setup!)
if suite == "none"
    @info "No tests selected. Exiting."
    exit()
end

# Other tests
LoggingExtras.global_logger(Logging.SimpleLogger(Logging.Warn))
include("test_consts.jl")
include("test_phys.jl")
if suite != "fast"
    include("test_integration.jl")
end
