# Tests for src/solver/linesearch.jl: step-length control for Newton-type solvers.
#
# Invariants and contract clauses exercised in this file:
#   - armijo: accepted steps satisfy φ(α) <= φ_ref + c1 α dφ0; the full step is tried first;
#             non-descent directions are refused; non-finite merits are backtracked from;
#             invalid parameters throw ArgumentError; failure returns the best point found.
#   - nonmonotone_ref: maximum over the last M merits; M=1 recovers the monotone condition.
#   - broyden_update!: the secant equation B s = y holds after the update; the update has
#             rank one; it vanishes when B already satisfies the secant equation.
#
# All problems are analytic, so no SOCRATES data is needed.

using Test
using LinearAlgebra

import AGNI.solver.linesearch: armijo, nonmonotone_ref, broyden_update!

# Scalar test problem r(x) = atan(x), with root x=0. Undamped Newton diverges from
#    |x0| > ~1.39, which is the classic example of why globalisation is needed.
_r(x)  = atan(x)
_dr(x) = 1.0 / (1.0 + x^2)

# Damped Newton iteration on r(x)=atan(x) with merit |r|, using the Armijo linesearch
function _armijo_newton(x0::Float64, nmax::Int)
    x = x0
    merits = [abs(_r(x))]
    alphas = Float64[]
    for _ in 1:nmax
        d = -_r(x) / _dr(x)
        φ(α) = abs(_r(x + α * d))
        α, _, _, ok = armijo(φ, abs(_r(x)), -abs(_r(x)))
        ok || break
        x += α * d
        push!(alphas, α)
        push!(merits, abs(_r(x)))
        merits[end] < 1e-12 && break
    end
    return x, merits, alphas
end

# Acceptance rule of the backtracking linesearch in solve_energy! (ls_method=2):
#    a trial is accepted if its cost is below ls_increase times the cost of the full step.
_legacy_accepts(c_trial, c_full, ls_increase) = c_trial <= c_full * ls_increase

@testset "linesearch" begin

    # Contract: Armijo-damped Newton converges on atan(x) from x0=1.5, where undamped
    #    Newton diverges, and the merit decreases monotonically along the way.
    @testset "armijo_globalises_newton_where_undamped_newton_diverges" begin
        x0 = 1.5    # beyond the ~1.39 radius of convergence of undamped Newton

        # undamped Newton overshoots to a point further from the root
        x1_undamped = x0 - _r(x0) / _dr(x0)
        @test abs(x1_undamped) > abs(x0)

        x, merits, alphas = _armijo_newton(x0, 30)
        @test abs(x) < 1e-10
        @test all(diff(merits) .< 0.0)          # monotone decrease
        @test alphas[1] < 1.0                   # the full first step was rejected
        @test isapprox(alphas[end], 1.0; atol=1e-12)  # full steps are taken near the root
    end

    # Discrimination guard for the legacy rule: from x0=3, the first backtracked step
    #    (α=0.495) increases the merit, yet the legacy rule accepts it with ls_increase=1.1
    #    (as used by PROTEUS in early coupling loops). Armijo must reject it.
    @testset "armijo_rejects_cost_increase_accepted_by_legacy_rule" begin
        x0 = 3.0
        d  = -_r(x0) / _dr(x0)
        c0     = abs(_r(x0))
        c_full = abs(_r(x0 + 0.99 * d))
        c_half = abs(_r(x0 + 0.99 * 0.5 * d))

        @test c_full > 1.1 * c0                       # legacy linesearch is triggered
        @test c_half > c0                             # trial is an increase in cost...
        @test _legacy_accepts(c_half, c_full, 1.1)    # ...which the legacy rule accepts
        @test !(c_half <= c0 + 1e-4 * 0.495 * (-c0))  # ...but violates Armijo

        φ(α) = abs(_r(x0 + α * d))
        α, φα, _, ok = armijo(φ, c0, -c0)
        @test ok
        @test φα < c0
        @test φα <= c0 * (1.0 - 1e-4 * α) + 1e-15
    end

    # Edge case: on a linear system the Newton step is exact, so the full step must be
    #    accepted after a single merit evaluation.
    @testset "armijo_accepts_exact_newton_step_on_linear_system" begin
        A  = [4.0 1.0; 2.0 3.0]
        bv = [1.0, -2.0]
        x0 = [10.0, -7.0]
        d  = -(A \ (A * x0 - bv))
        φ(α) = norm(A * (x0 .+ α .* d) .- bv)
        φ0 = norm(A * x0 .- bv)
        α, φα, nfev, ok = armijo(φ, φ0, -φ0)
        @test ok
        @test nfev == 1
        @test isapprox(α, 1.0; atol=1e-15)
        @test φα < 1e-12 * φ0
    end

    # Error contract: invalid parameters throw, and a non-descent direction is refused
    #    without evaluating the merit.
    @testset "armijo_refuses_invalid_parameters_and_ascent_directions" begin
        φ(α) = (1.0 - α)^2
        @test_throws ArgumentError armijo(φ, 1.0, -2.0; c1=0.0)
        @test_throws ArgumentError armijo(φ, 1.0, -2.0; c1=1.0)
        @test_throws ArgumentError armijo(φ, 1.0, -2.0; α0=-1.0)

        α, φα, nfev, ok = armijo(φ, 1.0, +2.0)
        @test !ok
        @test nfev == 0
        @test isapprox(α, 0.0; atol=0.0)
        @test isapprox(φα, 1.0; atol=0.0)
    end

    # Edge case: trials where the merit cannot be evaluated (non-finite) must be backtracked
    #    from; if no finite decrease exists, failure returns α=0 and the original merit.
    @testset "armijo_backtracks_from_nonfinite_merits" begin
        # merit is undefined beyond α=0.3, and a smooth descent below it
        φ(α) = (α > 0.3) ? Inf : (1.0 - α)^2
        α, φα, nfev, ok = armijo(φ, 1.0, -2.0)
        @test ok
        @test 0.0 < α <= 0.3
        @test φα <= 1.0 + 1e-4 * α * (-2.0)
        @test nfev >= 2

        # merit is never defined: must give up at α_min, not loop forever
        ψ(α) = NaN
        α2, φα2, nfev2, ok2 = armijo(ψ, 1.0, -2.0; α_min=1e-3)
        @test !ok2
        @test isapprox(α2, 0.0; atol=0.0)
        @test isapprox(φα2, 1.0; atol=0.0)
        @test 1 < nfev2 <= 30
    end

    # Contract: the non-monotone reference is the maximum over the last M merits. A trial
    #    which increases the merit is accepted relative to that reference, but not under
    #    the monotone condition (M=1).
    @testset "nonmonotone_reference_permits_bounded_increase" begin
        hist = [10.0, 2.0, 3.0]
        @test isapprox(nonmonotone_ref(hist, 3), 10.0; atol=0.0)
        @test isapprox(nonmonotone_ref(hist, 2), 3.0; atol=0.0)
        @test isapprox(nonmonotone_ref(hist, 1), 3.0; atol=0.0)
        @test isapprox(nonmonotone_ref(hist, 99), 10.0; atol=0.0)   # memory longer than history

        # merit along the line rises from 3 to 5 at the full step
        φ(α) = 3.0 + 2.0 * α
        _, _, _, ok_mono = armijo(φ, 3.0, -3.0; φ_ref=nonmonotone_ref(hist, 1), α_min=0.5)
        _, φα, _, ok_nm  = armijo(φ, 3.0, -3.0; φ_ref=nonmonotone_ref(hist, 3))
        @test !ok_mono
        @test ok_nm
        @test φα > 3.0          # an increase was accepted

        @test_throws ArgumentError nonmonotone_ref(hist, 0)
        @test_throws ArgumentError nonmonotone_ref(Float64[], 3)
    end

    # Contract: Broyden's update satisfies the secant equation with a rank-one change,
    #    and leaves B unchanged when it already maps s to y.
    @testset "broyden_update_satisfies_secant_equation_with_rank_one_change" begin
        B0 = [2.0 0.5 0.0; -1.0 3.0 0.2; 0.3 0.0 1.5]
        s  = [0.1, -0.2, 0.05]
        y  = [0.4, -0.9, 0.3]
        B  = copy(B0)
        @test broyden_update!(B, s, y)
        @test isapprox(B * s, y; atol=1e-14)
        # discrimination guard: the original matrix did not satisfy the secant equation
        @test !isapprox(B0 * s, y; atol=1e-3)
        sv = svdvals(B .- B0)
        @test sv[2] < 1e-14 * sv[1]          # rank-one change

        # consistent secant pair: no change
        B2 = copy(B0)
        @test broyden_update!(B2, s, B0 * s)
        @test isapprox(B2, B0; atol=1e-15)

        # zero step: update refused, matrix untouched
        B3 = copy(B0)
        @test !broyden_update!(B3, zeros(3), y)
        @test isapprox(B3, B0; atol=0.0)
        @test !broyden_update!(B3, s, [NaN, 0.0, 0.0])
    end
end
