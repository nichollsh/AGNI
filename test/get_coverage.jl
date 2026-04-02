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
using Printf
using Dates


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


# Calculate per-file statistics
file_stats = Dict{String, Dict{String, Any}}()
for c in coverage
    if !haskey(file_stats, c.filename)
        file_stats[c.filename] = Dict(
            "covered" => 0,
            "total" => 0,
            "uncovered_lines" => Int[],
            "covered_lines" => Int[]
        )
    end

    # c.coverage is a Vector{Union{Nothing, Int64}}
    for (line_num, cov_count) in enumerate(c.coverage)
        if cov_count !== nothing
            file_stats[c.filename]["total"] += 1
            if cov_count > 0
                file_stats[c.filename]["covered"] += 1
                push!(file_stats[c.filename]["covered_lines"], line_num)
            else
                push!(file_stats[c.filename]["uncovered_lines"], line_num)
            end
        end
    end
end

# Calculate overall statistics
total_covered = sum(stats["covered"] for stats in values(file_stats))
total_lines = sum(stats["total"] for stats in values(file_stats))
overall_pct = total_lines > 0 ? total_covered / total_lines * 100 : 0.0

# Sort files by coverage percentage (ascending) and then by name
sorted_files = sort(collect(file_stats),
    by = x -> (x[2]["total"] > 0 ? x[2]["covered"]/x[2]["total"] : 0.0, x[1]))

# Generate markdown report
output_file = joinpath(ROOT_DIR, "coverage.md")
open(output_file, "w") do io
    println(io, "# AGNI Code Coverage Report")
    println(io, "")
    println(io, "Generated: $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))")
    println(io, "")

    # Overall summary
    println(io, "## Overall Coverage")
    println(io, "")
    println(io, @sprintf("**Total Coverage: %.1f%%** (%d / %d lines)",
        overall_pct, total_covered, total_lines))
    println(io, "")

    # Per-file summary table
    println(io, "## Coverage by File")
    println(io, "")
    println(io, "| File | Coverage | Lines Covered | Total Lines | Uncovered Lines |")
    println(io, "|------|----------|---------------|-------------|-----------------|")

    for (filename, stats) in sorted_files
        pct = stats["total"] > 0 ? stats["covered"] / stats["total"] * 100 : 0.0
        short_name = replace(filename, ROOT_DIR * "/" => "")

        # Format uncovered line ranges
        uncov_lines = stats["uncovered_lines"]
        if length(uncov_lines) > 1
            line_summary = ""
            line_group = Int[]
            for i in eachindex(uncov_lines)[2:end]
                # find groups of lines
                if uncov_lines[i] == uncov_lines[i-1] + 1
                    # continue group
                    push!(line_group, uncov_lines[i])
                else
                    if length(line_group) > 2
                        line_summary *= string(line_group[1], "-", line_group[end], ", ")
                        line_group = Int[]
                    else
                        line_summary *= string(uncov_lines[i-1], ", ")
                    end
                end
            end
        elseif length(uncov_lines) == 1
            line_summary = string(uncov_lines[1])
        else
            line_summary = "None"
        end
        line_summary = rstrip(line_summary, (',', ' '))

        # Color code: <50% = 🔴, 50-80% = 🟡, >80% = 🟢
        icon = pct >= 80 ? "🟢" : (pct >= 50 ? "🟡" : "🔴")

        println(io, @sprintf("| %s | %s %.1f%% | %d | %d | %s |",
            short_name, icon, pct, stats["covered"], stats["total"], line_summary))
    end

    println(io, "")

    # Detailed breakdown of low-coverage files (<50%)
    println(io, "## Files Needing Attention (< 50% coverage)")
    println(io, "")

    low_coverage_files = [(f, s) for (f, s) in sorted_files
        if s["total"] > 0 && s["covered"] / s["total"] < 0.5]

    if length(low_coverage_files) == 0
        println(io, "✅ All files have >= 50% coverage!")
    else
        for (filename, stats) in low_coverage_files
            pct = stats["covered"] / stats["total"] * 100
            short_name = replace(filename, ROOT_DIR * "/" => "")
            println(io, "### $(short_name)")
            println(io, "")
            println(io, @sprintf("- **Coverage:** %.1f%% (%d / %d lines)",
                pct, stats["covered"], stats["total"]))
            println(io, "")
        end
    end

    # Quick wins: files with 0% coverage but small size
    println(io, "## Quick Wins (0% coverage, < 100 lines)")
    println(io, "")

    quick_wins = [(f, s) for (f, s) in sorted_files
        if s["total"] > 0 && s["covered"] == 0 && s["total"] < 100]

    if length(quick_wins) == 0
        println(io, "No small untested files found.")
    else
        println(io, "These files are currently untested but are small enough for quick test additions:")
        println(io, "")
        println(io, "| File | Total Lines |")
        println(io, "|------|-------------|")

        for (filename, stats) in sort(quick_wins, by=x->x[2]["total"])
            short_name = replace(filename, ROOT_DIR * "/" => "")
            println(io, @sprintf("| %s | %d |", short_name, stats["total"]))
        end
    end

    println(io, "")
    println(io, "---")
    println(io, "")
    println(io, "_Generated by `test/coverage_to_markdown.jl`_")
end

println("Coverage report written to: $output_file")
println(@sprintf("Overall coverage: %.1f%% (%d / %d lines)",
    overall_pct, total_covered, total_lines))
