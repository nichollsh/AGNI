using Test
using AGNI
using Logging

const SE = AGNI.solver.solve_energy
ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)), "../"))
OUT_DIR  = joinpath(ROOT_DIR, "out/")

# Cheap, fully-allocated greygas atmosphere for directly unit testing _fev!/
# _calc_jac_res! (now standalone functions), without needing a real SOCRATES
# solve.
function _cheap_atmos(; nlev_centre::Int=16)
    atmos = atmosphere.Atmos_t()
    ok = atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                            "greygas",
                            1000.0, 1.0, 0.0, 60.0,
                            500.0,
                            10.0, 6.37e6,
                            nlev_centre, 100.0, 1e-6,
                            Dict("H2O" => 1.0), "";
                            real_gas=false,
                            thermo_functions=false,
                            flag_rayleigh=false,
                            flag_cloud=false)
    ok || error("Failed to setup test atmosphere")
    atmosphere.allocate!(atmos, ""; check_safe_gas=false) ||
        error("Failed to allocate test atmosphere")
    return atmos
end

# Solution-array guess built from the atmosphere's current temperatures,
# sized correctly for the given sol_type (cell-centres only for sol_type=1,
# plus a trailing surface/skin temperature for sol_type in [2,3,4]).
function _x_for(atmos::atmosphere.Atmos_t, sol_type::Int64)
    arr_len = sol_type == 1 ? atmos.nlev_c : atmos.nlev_c + 1
    x = zeros(Float64, arr_len)
    for i in 1:atmos.nlev_c
        x[i] = atmos.tmp[i]
    end
    if sol_type >= 2
        x[end] = atmos.tmp[atmos.nlev_c] + 1.0
    end
    return x
end

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

    # _fev! residual formulas differ by sol_type; each is checked against the
    # exact formula from the source (not just "is populated"), since these
    # branches (sol_type 1/2/4) are never exercised by the SOCRATES-based
    # integration test, which always uses sol_type=3.
    @testset "fev_residual_formulas" begin
        step_ok = Ref(true)
        code = Ref(SE.CODE_99)

        # sol_type=1: residual is exactly the flux-divergence array, no
        # surface term appended
        atmos1 = _cheap_atmos()
        x1 = _x_for(atmos1, 1)
        resid1 = zeros(Float64, length(x1))
        ok1 = SE._fev!(atmos1, x1, resid1, 1, true, false, true,
                        true, true, true, true, true, 1.0, step_ok, code)
        @test ok1
        @test resid1 == atmos1.flux_dif
        atmosphere.deallocate!(atmos1)

        # sol_type=2: last element uses the conductive skin flux instead of
        # the flux-divergence value a sol_type=1/3/4 formula would give there
        atmos2 = _cheap_atmos()
        x2 = _x_for(atmos2, 2)
        resid2 = zeros(Float64, length(x2))
        ok2 = SE._fev!(atmos2, x2, resid2, 2, true, false, true,
                        true, true, true, true, true, 1.0, step_ok, code)
        @test ok2
        @test resid2[1:end-1] == atmos2.flux_dif[1:end]
        @test isapprox(resid2[end], atmos2.flux_tot[end] - energy.skin_flux(atmos2); atol=1e-10)
        # discrimination guard: the skin-flux term must actually differ from
        # the flux_dif value a wrong (sol_type=1-style) formula would give
        @test !isapprox(resid2[end], atmos2.flux_dif[end]; atol=1e-6)
        atmosphere.deallocate!(atmos2)

        # sol_type=4: last element compares outgoing LW flux against a target
        atmos4 = _cheap_atmos()
        x4 = _x_for(atmos4, 4)
        resid4 = zeros(Float64, length(x4))
        atmos4.target_olr = 123.456  # arbitrary, distinguishable target
        ok4 = SE._fev!(atmos4, x4, resid4, 4, true, false, true,
                        true, true, true, true, true, 1.0, step_ok, code)
        @test ok4
        @test resid4[1:end-1] == atmos4.flux_dif[1:end]
        @test isapprox(resid4[end], 123.456 - atmos4.flux_u_lw[1]; atol=1e-10)
        atmosphere.deallocate!(atmos4)
    end

    # edge case: a non-finite residual (forced via a NaN target_olr) must be
    # caught and reported as CODE_NAN, not silently propagated
    @testset "fev_nan_residual_failure" begin
        atmos = _cheap_atmos()
        x = _x_for(atmos, 4)
        resid = zeros(Float64, length(x))
        atmos.target_olr = NaN
        step_ok = Ref(true)
        code = Ref(SE.CODE_99)
        logs, ok = Test.collect_test_logs() do
            SE._fev!(atmos, x, resid, 4, true, false, true,
                        true, true, true, true, true, 1.0, step_ok, code)
        end
        @test !ok
        @test code[] == SE.CODE_NAN
        @test any(occursin("NaNs and/or Infs", l.message) for l in logs if l.level == Logging.Warn)
        atmosphere.deallocate!(atmos)
    end

    # _calc_jac_res!: finite-difference branches not exercised by the
    # integration test (which always uses central=true, order=2, and
    # perturbs every level).
    @testset "calc_jac_res_finite_difference_branches" begin
        step_ok = Ref(true)
        code = Ref(SE.CODE_99)

        # which[i]=false: that column must be left untouched (still zero),
        # while a column that IS updated must be nonzero
        atmos = _cheap_atmos()
        x = _x_for(atmos, 1)
        arr_len = length(x)
        jacob = zeros(Float64, arr_len, arr_len)
        resid = zeros(Float64, arr_len)
        which = fill(true, arr_len)
        which[1] = false
        ok = SE._calc_jac_res!(atmos, x, jacob, resid, true, 2, which,
                                1, true, false, true, true, true, true, true, true,
                                1.0, 5.0, step_ok, code)
        @test ok
        @test all(jacob[:, 1] .== 0.0)
        @test any(jacob[:, 2] .!= 0.0)
        atmosphere.deallocate!(atmos)

        # central difference, 4th order: exercises the +-2*fd_s evaluations
        # and the 4th-order central-difference jacobian formula
        atmos4c = _cheap_atmos()
        x4c = _x_for(atmos4c, 1)
        arr_len4c = length(x4c)
        jacob4c = zeros(Float64, arr_len4c, arr_len4c)
        resid4c = zeros(Float64, arr_len4c)
        ok4c = SE._calc_jac_res!(atmos4c, x4c, jacob4c, resid4c, true, 4, fill(true, arr_len4c),
                                1, true, false, true, true, true, true, true, true,
                                1.0, 5.0, step_ok, code)
        @test ok4c
        @test any(jacob4c .!= 0.0)
        atmosphere.deallocate!(atmos4c)

        # forward (non-central) difference, both 2nd and 4th order
        for order in (2, 4)
            atmosf = _cheap_atmos()
            xf = _x_for(atmosf, 1)
            arr_lenf = length(xf)
            jacobf = zeros(Float64, arr_lenf, arr_lenf)
            residf = zeros(Float64, arr_lenf)
            okf = SE._calc_jac_res!(atmosf, xf, jacobf, residf, false, order, fill(true, arr_lenf),
                                    1, true, false, true, true, true, true, true, true,
                                    1.0, 5.0, step_ok, code)
            @test okf
            @test any(jacobf .!= 0.0)
            atmosphere.deallocate!(atmosf)
        end
    end

    # solve_energy!() failure/config paths not exercised by the SOCRATES
    # integration test (which always converges via method=1/Newton-Raphson
    # with default ls_method/conv_type). modplot=0 keeps these cheap by
    # skipping in-loop plotting; the post-loop plot on CODE_ITE still runs
    # once, since that's part of what's being covered here.
    @testset "solve_energy_failure_and_method_paths" begin
        # CODE_ITE: max_steps reached before convergence
        atmos_ite = _cheap_atmos()
        logs, converged = Test.collect_test_logs() do
            SE.solve_energy!(atmos_ite; sol_type=1, max_steps=1, modplot=0, save_frames=false)
        end
        @test !converged
        @test any(occursin("maximum iterations", l.message) for l in logs if l.level == Logging.Warn)
        atmosphere.deallocate!(atmos_ite)

        # CODE_TIM: max_runtime exceeded before the first step can complete
        atmos_tim = _cheap_atmos()
        logs, converged = Test.collect_test_logs() do
            SE.solve_energy!(atmos_tim; sol_type=1, max_runtime=0.0, modplot=0, save_frames=false)
        end
        @test !converged
        @test any(occursin("maximum time", l.message) for l in logs if l.level == Logging.Warn)
        atmosphere.deallocate!(atmos_tim)

        # CODE_CFG via an invalid linesearch algorithm choice
        atmos_ls = _cheap_atmos()
        logs, converged = Test.collect_test_logs() do
            SE.solve_energy!(atmos_ls; sol_type=1, ls_method=99, max_steps=2, modplot=0, save_frames=false)
        end
        @test !converged
        @test any(occursin("failure (configuration)", l.message) for l in logs if l.level == Logging.Warn)
        atmosphere.deallocate!(atmos_ls)

        # CODE_CFG via an invalid convergence-metric choice (a different
        # source line than the linesearch case above)
        atmos_cv = _cheap_atmos()
        converged_cv = SE.solve_energy!(atmos_cv; sol_type=1, conv_type=99, max_steps=2,
                                            modplot=0, save_frames=false)
        @test !converged_cv
        atmosphere.deallocate!(atmos_cv)

        # method=2 (Gauss-Newton) and method=4 (Jacobi-preconditioned
        # Newton): alternative solver methods, never selected by any other
        # test's config
        for method in (2, 4)
            atmos_m = _cheap_atmos()
            SE.solve_energy!(atmos_m; sol_type=1, method=method, max_steps=1,
                                modplot=0, save_frames=false)
            atmosphere.deallocate!(atmos_m)
        end

        # partial Jacobian update (perturb_all=false, step > 2): only reached
        # once the solver has taken a few steps and isn't near convergence;
        # the integration test never exercises this since it always uses
        # perturb_all=true
        atmos_pj = _cheap_atmos()
        SE.solve_energy!(atmos_pj; sol_type=1, perturb_all=false, max_steps=4,
                            modplot=0, save_frames=false)
        @test atmos_pj.is_solved
        atmosphere.deallocate!(atmos_pj)

        # method=3 (Levenberg-Marquardt): documents a pre-existing bug (an
        # undefined `dtd` variable in this branch) rather than fixing it -
        # pins the current (crashing) behaviour so a future change is deliberate
        atmos_lm = _cheap_atmos()
        @test_throws UndefVarError SE.solve_energy!(atmos_lm; sol_type=1, method=3,
                                                        max_steps=1, modplot=0, save_frames=false)
        atmosphere.deallocate!(atmos_lm)
    end
end
