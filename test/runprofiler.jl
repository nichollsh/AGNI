#!/usr/bin/env -S julia --color=yes --startup-file=no
# Run this tool from inside the AGNI root folder
# e.g. as `julia --project test/runprofiler.jl`

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)), "../"))

# Activate environment
ENV["GKSwstype"] = "100"
import Pkg
Pkg.activate(ROOT_DIR)

using Profile
using AGNI

# Help message
function _print_help()
    println("AGNI profiler tool")
    println("")
    println("Profiles a single AGNI run using a TOML config and writes an HTML flame report.")
    println("")
    println("Usage:")
    println("  julia --project test/runprofiler.jl [config_path] [report_dir]")
    println("")
    println("Defaults:")
    println("  config_path = test/test.toml")
    println("  report_dir  = profile_report/")
end

# Check that library is installed
function _load_profiler_library()
    try
        @eval using StatProfilerHTML
        return true
    catch err
        @error "StatProfilerHTML.jl is required. Install it once with: julia --project -e 'using Pkg; Pkg.add(\"StatProfilerHTML\")'" exception=(err, catch_backtrace())
        return false
    end
end

# Main function
function do_profiling(cfg_path, report_dir)

    isfile(cfg_path) || error("Config not found: $cfg_path")
    mkpath(dirname(report_dir))

    cfg = AGNI.open_config(cfg_path)
    output_dir = cfg["files"]["output_dir"]

    @info "Profiling AGNI.run_from_config using config: $cfg_path"
    @info "Model outputs will be written to: $output_dir"
    @info "Profiler report directory: $report_dir"

    Profile.clear()
    Profile.Allocs.clear()
    try
        Profile.@profile begin
            succ = AGNI.run_from_config(cfg)
            succ || error("AGNI run_from_config returned false")
        end
    catch ex
        ex isa InterruptException || rethrow(ex)
        @info "Interrupted; writing profile for samples collected so far."
    end

    # Write report
    StatProfilerHTML.statprofilehtml(path=report_dir)

    #
end

# If executed directly...
if abspath(PROGRAM_FILE) == @__FILE__

    # Parse arguments
    if any(arg -> arg in ("-h", "--help"), ARGS)
        _print_help()
        exit(0)
    end

    # Check that profiler library is available
    if !_load_profiler_library()
        exit(1)
    end

    # Get paths
    cfg_path = length(ARGS) >= 1 ? abspath(ARGS[1]) : joinpath(ROOT_DIR, "test", "test.toml")
    report_dir = length(ARGS) >= 2 ? abspath(ARGS[2]) : joinpath(ROOT_DIR, "profile_report")

    # Main profiling function
    do_profiling(cfg_path, report_dir)
end
