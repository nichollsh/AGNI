using Test
using AGNI

const SE = AGNI.solver.solve_energy

@testset "solver_energy" begin

    # only CODE_SUC reports convergence, and must not also request a plot
    @testset "classify_status_success" begin
        converged, level, message, should_plot = SE._classify_status(SE.CODE_SUC, 7)
        @test converged
        @test level == :info
        @test occursin("7 steps", message)
        @test !should_plot
    end

    # every failure code must report non-converged at :warn level
    @testset "classify_status_failure_codes" begin
        cases = [
            (SE.CODE_ITE, true),
            (SE.CODE_SIN, true),
            (SE.CODE_TIM, false),
            (SE.CODE_NAN, true),
            (SE.CODE_CFG, false),
            (SE.CODE_OBJ, true),
            (SE.CODE_STP, true),
            (SE.CODE_HYD, true),
        ]
        for (code, expect_plot) in cases
            converged, level, message, should_plot = SE._classify_status(code, 99)
            @test !converged
            @test level == :warn
            @test !isempty(message)
            @test should_plot == expect_plot
        end

        # ensure that the failure codes all have their own unique messages
        messages = [SE._classify_status(code, 0)[3] for (code, _) in cases]
        @test length(unique(messages)) == length(messages)
    end

    # edge case which matches none of status code types
    @testset "classify_status_default_fallback" begin
        converged, level, message, should_plot = SE._classify_status(SE.CODE_99, 3)
        @test !converged
        @test level == :warn
        @test should_plot
        @test occursin("other", message)
    end
end
