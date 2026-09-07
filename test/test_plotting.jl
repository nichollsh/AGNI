# Tests for src/interface/plotting.jl
# Covers:
#   - plotting._savefig_safe(): the shared savefig helper used by every plot_*()
#     function, guarding against an empty filename and a missing output
#     directory (the insufficient-disk-space branch is not exercised here,
#     since it can only be hit by a genuinely near-full filesystem).

using Test
using AGNI
using Plots
using Logging

@testset "plotting" begin

    @testset "savefig_safe" begin
        plt = plot([1.0, 2.0, 3.0], [1.0, 4.0, 9.0])

        # empty filename: a no-op, used by plot_*() functions when the
        # caller only wants the Plots.Plot object back, not a saved file
        logs, _ = Test.collect_test_logs() do
            plotting._savefig_safe(plt, "")
        end
        @test isempty([l for l in logs if l.level >= Logging.Warn])

        tmpdir = mktempdir()

        # nonexistent output directory: warns, does not create the file
        missing_dir_path = joinpath(tmpdir, "does_not_exist", "fig.png")
        logs, _ = Test.collect_test_logs() do
            plotting._savefig_safe(plt, missing_dir_path)
        end
        @test !isfile(missing_dir_path)
        @test any(occursin("Directory does not exist", l.message)
                        for l in logs if l.level == Logging.Warn)

        # happy path: directory exists and there's ample disk space, so the
        # figure is actually written to disk
        fig_path = joinpath(tmpdir, "fig.png")
        @test !isfile(fig_path)
        plotting._savefig_safe(plt, fig_path)
        @test isfile(fig_path)

        rm(tmpdir; force=true, recursive=true)
    end
end
