using Test
using AGNI
using Logging
using LoggingExtras
using LinearAlgebra

const SE = AGNI.solver.solve_energy
ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)), "../"))
OUT_DIR  = joinpath(ROOT_DIR, "out/")

# Cheap, fully-allocated greygas atmosphere for directly unit testing _fev!/
# _calc_jac_res! (now standalone functions), without needing a real SOCRATES
# solve.
function _cheap_atmos(; nlev_centre::Int=16,
                        tmp_floor::Float64=atmosphere.CFG_tmp_floor,
                        hill_radius::Float64=atmosphere.CFG_hill_radius)
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
                            flag_cloud=false,
                            tmp_floor=tmp_floor,
                            hill_radius=hill_radius)
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
                                1.0, 5.0, step_ok, code, 1e-4, 1e-5, false)
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
                                1.0, 5.0, step_ok, code, 1e-4, 1e-5, false)
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
                                    1.0, 5.0, step_ok, code, 1e-4, 1e-5, false)
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

        # method=2 (Gauss-Newton): alternative solver method, never selected
        # by any other test's config
        atmos_m = _cheap_atmos()
        SE.solve_energy!(atmos_m; sol_type=1, method=2, max_steps=1,
                            modplot=0, save_frames=false)
        @test atmos_m.is_solved
        atmosphere.deallocate!(atmos_m)

        # invalid options are configuration errors, rejected before solving:
        #    method=4 (the former Jacobi-preconditioned Newton, which was a
        #    duplicate of Gauss-Newton) no longer exists, the non-monotone
        #    memory must be at least one, and the Broyden refresh interval
        #    must be non-negative
        for (kw, msg) in (((method=4,),           "Invalid solver method"),
                          ((ls_memory=0,),        "Invalid linesearch memory"),
                          ((broyden_refresh=-1,), "Invalid Broyden refresh"))
            atmos_bad = _cheap_atmos()
            logs, converged = Test.collect_test_logs() do
                SE.solve_energy!(atmos_bad; sol_type=1, max_steps=2, modplot=0,
                                    save_frames=false, kw...)
            end
            @test !converged
            @test !atmos_bad.is_solved
            @test any(occursin(msg, l.message) for l in logs if l.level == Logging.Warn)
            atmosphere.deallocate!(atmos_bad)
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

        # method=3 (Levenberg-Marquardt): previously raised UndefVarError (an
        # undefined `dtd` variable); a single step must now run and report that
        # the iteration limit was reached
        atmos_lm = _cheap_atmos()
        logs, converged = Test.collect_test_logs() do
            SE.solve_energy!(atmos_lm; sol_type=1, method=3, max_steps=1, modplot=0,
                                save_frames=false)
        end
        @test !converged
        @test atmos_lm.is_solved
        @test any(occursin("maximum iterations", l.message) for l in logs if l.level == Logging.Warn)
        atmosphere.deallocate!(atmos_lm)
    end

    # Convergence is only accepted on steps which were not rejected, except for when
    # they are flagged as unbound
    @testset "convergence_requires_accepted_step" begin
        conv = (conv_atol=1e30, conv_rtol=0.0, max_steps=3, modplot=0, save_frames=false)

        # baseline: an accepted step converges immediately
        atmos_ok = _cheap_atmos()
        converged = with_logger(MinLevelLogger(current_logger(), Logging.Error+1)) do
            SE.solve_energy!(atmos_ok; sol_type=1, conv...)
        end
        @test converged
        @test atmos_ok.is_converged
        atmosphere.deallocate!(atmos_ok)

        # every step is rejected by the surface-temperature collapse check (the fixed
        #    500 K surface is within tmp_pad+1 of the floor), so the solver must not report
        #    convergence, despite the criterion being satisfied. It instead runs out of steps.
        atmos_rej = _cheap_atmos(; tmp_floor=495.0)
        logs, converged = Test.collect_test_logs() do
            SE.solve_energy!(atmos_rej; sol_type=1, conv...)
        end
        @test !converged
        @test !atmos_rej.is_converged
        @test any(occursin("maximum iterations", l.message) for l in logs if l.level == Logging.Warn)
        atmosphere.deallocate!(atmos_rej)

        # layers above a Hill radius placed within the atmosphere are unbound on every step,
        #    but the solver can still converge, and warns that it has done so
        base = _cheap_atmos()
        atmosphere.calc_layer_props!(base)
        r_hill = 0.5 * (base.rl[end] + base.rl[1])
        atmosphere.deallocate!(base)
        atmos_unb = _cheap_atmos(; hill_radius=r_hill)
        logs, converged = Test.collect_test_logs() do
            SE.solve_energy!(atmos_unb; sol_type=1, conv...)
        end
        @test converged
        @test any(.!atmos_unb.layer_isbound)
        @test any(occursin("not bound by gravity", l.message) for l in logs if l.level == Logging.Warn)
        atmosphere.deallocate!(atmos_unb)

        # unbound layers must not make a step acceptable when it was also rejected for
        #    another reason (here, the surface-temperature collapse on every step)
        atmos_both = _cheap_atmos(; tmp_floor=495.0, hill_radius=r_hill)
        logs, converged = Test.collect_test_logs() do
            SE.solve_energy!(atmos_both; sol_type=1, conv...)
        end
        @test !converged
        @test any(occursin("maximum iterations", l.message) for l in logs if l.level == Logging.Warn)
        @test !any(occursin("not bound by gravity", l.message) for l in logs if l.level == Logging.Warn)
        atmosphere.deallocate!(atmos_both)
    end

    # Shared settings for the tests below: tolerances tight enough that each method
    #    must make real progress from the fixture's isothermal initial state
    solve_kw = (modplot=0, save_frames=false, conv_atol=1e-3, conv_rtol=1e-4, max_steps=100)
    function _solve(; sol_type::Int64=1, kw...)
        atmos = _cheap_atmos()
        logs, ok = Test.collect_test_logs(; min_level=Logging.Info) do
            SE.solve_energy!(atmos; sol_type=sol_type, solve_kw..., kw...)
        end
        tmp  = copy(atmos.tmp)
        nrt  = atmos.num_rt_eval
        atmosphere.deallocate!(atmos)
        return ok, tmp, nrt, [string(l.message) for l in logs]
    end

    # The Levenberg-Marquardt step solves (b'b + λD) dx = -b'r. As λ → 0 it must equal
    #    the Newton step, its D-scaled length must shrink monotonically as λ grows, and for
    #    large λ it must turn towards the scaled gradient direction -D⁻¹b'r/λ.
    @testset "levenberg_marquardt_step_interpolates_newton_and_gradient_directions" begin
        J = [4.0 1.0 0.0; 1.0 3.0 1.0; 0.0 2.0 5.0]    # well-conditioned, non-symmetric
        r = [1.0, -2.0, 0.5]
        btb = J' * J
        btr = J' * r
        d   = diag(btb)
        dx  = zeros(3)
        fac = zeros(3, 3)

        # λ = 0 reproduces the Newton step exactly (up to rounding)
        SE._lm_step!(dx, btb, btr, 0.0, d, fac)
        dx_newton = -(J \ r)
        @test isapprox(dx, dx_newton; rtol=1e-10)

        # D-scaled step length decreases monotonically with λ
        lens = Float64[]
        for λ in (1e-4, 1e-2, 1.0, 1e2, 1e4)
            SE._lm_step!(dx, btb, btr, λ, d, fac)
            push!(lens, sqrt(sum(d .* dx .^ 2)))
        end
        @test all(diff(lens) .< 0.0)

        # large λ: step tends to the scaled gradient direction
        λ = 1e8
        SE._lm_step!(dx, btb, btr, λ, d, fac)
        dx_grad = -btr ./ d ./ λ
        @test isapprox(dx, dx_grad; rtol=1e-6)
        # discrimination guard: the gradient and Newton directions differ markedly here
        @test !isapprox(normalize(dx_grad), normalize(dx_newton); atol=0.1)

        # error contract: an undamped singular system cannot be solved, while any positive
        #    damping with a positive scaling makes it solvable
        S = [1.0 2.0; 2.0 4.0]
        @test_throws SingularException SE._lm_step!(zeros(2), S, [1.0, 1.0], 0.0, [1.0, 1.0], zeros(2, 2))
        dx_s = zeros(2)
        SE._lm_step!(dx_s, S, [1.0, 1.0], 1.0, [1.0, 4.0], zeros(2, 2))
        @test all(isfinite, dx_s)
    end

    # Levenberg-Marquardt (method 3) must converge on the fixture, to the same profile as
    #    Newton-Raphson in the deep levels, which are strongly coupled to the fixed surface.
    #    (Levels near the top are optically thin, so weakly constrained by the residuals.)
    @testset "levenberg_marquardt_converges_to_newton_solution" begin
        ok_lm, tmp_lm, _, logs_lm = _solve(sol_type=1, method=3)
        ok_nr, tmp_nr, _, _       = _solve(sol_type=1, method=1)
        @test ok_lm && ok_nr
        @test isapprox(tmp_lm[end-3:end], tmp_nr[end-3:end]; atol=0.5)
        @test any(occursin("Lm", m) for m in logs_lm)   # LM steps were taken
    end

    # When steps are clipped by dx_max, the backtracking linesearch (ls_method=2) compares
    #    trials with the cost of the full step, and so rejects steps which genuinely reduce
    #    the cost. The non-monotone Armijo linesearch (ls_method=3) must accept them, and so
    #    converge to the same solution with far fewer RT evaluations.
    @testset "nonmonotone_armijo_accepts_clipped_steps_rejected_by_backtracking" begin
        ok3, tmp3, nrt3, _ = _solve(sol_type=1, ls_method=3, dx_max=35.0)
        ok2, tmp2, nrt2, _ = _solve(sol_type=1, ls_method=2, dx_max=35.0)
        @test ok3 && ok2
        @test nrt3 < 0.5 * nrt2       # measured ratio is about 0.3
        @test isapprox(tmp3, tmp2; atol=0.1)

        # free surface temperature (sol_type=3), and monotone Armijo (ls_memory=1) as an
        #    edge case of the non-monotone condition
        ok3s, tmp3s, _, _ = _solve(sol_type=3, ls_method=3)
        ok1s, tmp1s, _, _ = _solve(sol_type=3, ls_method=3, ls_memory=1)
        @test ok3s && ok1s
        @test isapprox(tmp3s, tmp1s; atol=0.1)
    end

    # Broyden updates replace most finite-difference Jacobians, so the solver must converge
    #    with fewer RT evaluations than when finite-differencing at every step. With
    #    broyden_refresh=0 every step is finite-differenced, which must exactly reproduce
    #    the solver without Broyden updates.
    @testset "broyden_updates_reduce_rt_evaluations" begin
        ok_f, tmp_f, nrt_f, logs_f = _solve(sol_type=3, ls_method=3, broyden=false)
        ok_b, tmp_b, nrt_b, logs_b = _solve(sol_type=3, ls_method=3, broyden=true)
        ok_0, tmp_0, nrt_0, logs_0 = _solve(sol_type=3, ls_method=3, broyden=true,
                                            broyden_refresh=0)
        @test ok_f && ok_b && ok_0
        @test nrt_b < nrt_f
        @test isapprox(tmp_b, tmp_f; atol=0.1)
        @test any(occursin("Br", m) for m in logs_b)     # Broyden steps were taken
        @test !any(occursin("Br", m) for m in logs_f)
        @test !any(occursin("Br", m) for m in logs_0)
        @test nrt_0 == nrt_f
        @test isapprox(tmp_0, tmp_f; rtol=1e-12)
    end

    # Regression: with hundreds of levels, the running product of LU pivots underflows to
    #    zero, so the former determinant-based test reported this well-conditioned Jacobian
    #    as singular and the solver stopped at its first step. Singularity must instead be
    #    judged from the factorisation itself.
    @testset "singular_jacobian_check_is_not_fooled_by_determinant_underflow" begin
        atmos = _cheap_atmos(nlev_centre=400)    # 400 levels triggers the underflow
        x = _x_for(atmos, 1)
        n = length(x)
        J = zeros(Float64, n, n)
        r = zeros(Float64, n)
        SE._calc_jac_res!(atmos, x, J, r, true, 2, fill(true, n),
                            1, true, false, true, true, true, true, true, true,
                            1.0, 5.0, Ref(true), Ref(SE.CODE_99), 1e-4, 1e-5, false)

        # discrimination guard: the determinant underflows, so the former test fails...
        @test abs(det(J)) < floatmin() * 10.0
        # ...although |det J| is mathematically far above that threshold...
        @test logabsdet(J)[1] > log(floatmin() * 10.0) + 100.0
        # ...and the matrix is far from singular to working precision
        @test cond(J) < 1e10
        @test !SE._is_singular(J)

        # the solver must now proceed beyond the first step
        logs, _ = Test.collect_test_logs() do
            SE.solve_energy!(atmos; sol_type=1, max_steps=2, modplot=0, save_frames=false)
        end
        warns = [l.message for l in logs if l.level == Logging.Warn]
        @test any(occursin("maximum iterations", m) for m in warns)
        @test !any(occursin("singular jacobian", m) for m in warns)
        atmosphere.deallocate!(atmos)

        # error contract: exactly singular and non-finite matrices are still detected
        @test SE._is_singular(zeros(Float64, 3, 3))
        @test SE._is_singular([1.0 2.0; 2.0 4.0])
        @test SE._is_singular([1.0 Inf; 0.0 1.0])
        @test !SE._is_singular([2.0 1.0; 1.0 3.0])
    end
end
