using Test
using AGNI

# Access the guillot submodule through setpt
const guillot = AGNI.setpt.guillot

# Guillot (2010) analytical grey-gas T(p) profile module.
# Constants inside the module:
#   κ_th = 7e-2 * 10 = 0.7  m2/kg  (LW opacity)
#   κ_vs = 4e-3 * 10 = 0.04 m2/kg  (SW opacity)
#   γ = κ_vs/κ_th = 0.04/0.7       (opacity ratio)

@testset "guillot" begin

    # -------------
    # Optical depth: τ = p * κ_th / g
    # -------------
    @testset "eval_tau" begin
        p_test = [1e3,   1e5,      1e7     ]  # pressure [Pa]
        g_test = [10.0,  10.0,     25.0    ]  # gravity [m s-2]
        v_expt = [70.0,  7000.0,   280000.0]  # τ = p * 0.7 / g [dimensionless]
        for i in eachindex(p_test)
            @test isapprox(guillot.eval_tau(p_test[i], g_test[i]), v_expt[i]; rtol=1e-10)
        end
        # τ must be linear in pressure at fixed gravity
        @test isapprox(guillot.eval_tau(2e5, 10.0), 2 * guillot.eval_tau(1e5, 10.0); rtol=1e-10)
    end

    # -------------
    # Stellar temperatures: geometric relations
    #   eval_Tirr(Tstar, Rstar, sep) = Tstar * sqrt(Rstar / sep)
    #   eval_Teqm(Tstar, Rstar, sep) = Tstar * sqrt(Rstar / (2*sep))
    # So Tirr / Teqm == sqrt(2) exactly.
    # -------------
    @testset "stellar_temps" begin
        Tstar = 5778.0    # solar effective temperature [K]
        Rstar = 6.957e8   # solar radius [m]
        sep   = 1.496e11  # Earth-Sun distance [m]

        Tirr = guillot.eval_Tirr(Tstar, Rstar, sep)
        Teqm = guillot.eval_Teqm(Tstar, Rstar, sep)

        @test Tirr > 0.0
        @test Teqm > 0.0
        @test Teqm < Tirr
        # exact geometric relationship between the two temperatures
        @test isapprox(Tirr / Teqm, sqrt(2.0); rtol=1e-10)
    end

    # -------------
    # T^4 functions must return positive values for physical inputs
    # Note: eval_T4_avg uses E2(γτ) which diverges at τ=0 (τ>0 only).
    # eval_T4_cos does not use E1/E2 so τ=0 is fine.
    # -------------
    @testset "T4_positive" begin
        τ_vals_avg = [0.01, 0.1, 1.0, 10.0, 100.0]  # τ > 0 for avg (E1 singularity at 0)
        τ_vals_cos = [0.0,  0.1, 1.0, 10.0, 100.0]  # τ = 0 is fine for cos
        Tint = 100.0    # internal temperature [K]
        Teqm = 800.0    # equilibrium temperature [K]
        Tirr = 1000.0   # irradiation temperature [K]
        θ    = 45.0     # zenith angle [degrees]
        for τ in τ_vals_avg
            @test guillot.eval_T4_avg(τ, Tint, Teqm) > 0.0
        end
        for τ in τ_vals_cos
            @test guillot.eval_T4_cos(τ, Tint, Tirr, θ) > 0.0
        end
    end

    # -------------
    # T^4 must increase monotonically with optical depth (deeper = hotter)
    # for the average profile at sensible parameters
    # -------------
    @testset "T4_monotone" begin
        τ_lo, τ_hi = 0.1, 10.0
        Tint, Teqm, Tirr, θ = 200.0, 1200.0, 1500.0, 30.0
        @test guillot.eval_T4_avg(τ_hi, Tint, Teqm) > guillot.eval_T4_avg(τ_lo, Tint, Teqm)
        @test guillot.eval_T4_cos(τ_hi, Tint, Tirr, θ) > guillot.eval_T4_cos(τ_lo, Tint, Tirr, θ)
    end

    # -------------
    # Profile helpers: correct length and all-positive temperatures
    # -------------
    @testset "profiles" begin
        grav  = 10.0   # [m s-2]
        Tint  = 100.0  # [K]
        Teqm  = 900.0  # [K]
        Tirr  = 1200.0 # [K]
        θ     = 45.0   # [degrees]
        p_arr = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6]  # pressure grid [Pa]

        T_avg = guillot.calc_profile_avg(p_arr, grav, Tint, Teqm)
        @test length(T_avg) == length(p_arr)
        @test all(T_avg .> 0.0)
        @test issorted(T_avg)  # temperature increases towards higher pressures

        T_cos = guillot.calc_profile_cos(p_arr, grav, Tint, Tirr, θ)
        @test length(T_cos) == length(p_arr)
        @test all(T_cos .> 0.0)
        @test issorted(T_cos)
    end

end
