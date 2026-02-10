#!/usr/bin/env -S julia --color=yes --startup-file=no
# Run this function from inside the AGNI root folder
# e.g. as `julia --project=. tests/get_coverage.jl `

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))

# Activate environment
ENV["GKSwstype"] = "100"
import Pkg
Pkg.activate(ROOT_DIR)

# https://github.com/JuliaCI/Coverage.jl

# Include libraries
using Coverage

# process '*.cov' files
coverage = process_folder() # defaults to src/; alternatively, supply the folder name as argument
coverage = append!(coverage, process_folder("deps"))  # useful if you want to analyze more than just src/

# process '*.info' files, if you collected them
coverage = merge_coverage_counts(coverage, filter!(
    let prefixes = (joinpath(pwd(), "src", ""),
                    joinpath(pwd(), "deps", ""))
        c -> any(p -> startswith(c.filename, p), prefixes)
    end,
    LCOV.readfolder("test")))

# Get total coverage for all Julia files
covered_lines, total_lines = get_summary(coverage)
coverage_pct = round(covered_lines / total_lines * 100, digits=1)
@info "Total coverage: $covered_lines of $total_lines ($coverage_pct% coverage)"

# Write total coverage to single file
open("coverage.total", "w") do io
    write(io, "$coverage_pct")
end

# Write to files that CI can read
LCOV.writefile("coverage.info", coverage)
export_codecov_json(coverage, "coverage.json")

# Or process a single file
# @show get_summary(process_file(joinpath("src", "AGNI.jl")))
