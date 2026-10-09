using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)), "../"))
OUT_DIR = joinpath(ROOT_DIR, "out/")

# helper to create a simple grey gas atmosphere for testing diagnostics and energy utilities
function _make_greygas_atmos(; instellation::Float64=1200.0,
                               nlev_c::Int64=30,
                               p_surf::Float64=10.0,
                               p_top::Float64=1e-6)
    atmos = atmosphere.Atmos_t()
    ok = atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                            "greygas",
                            instellation, 1.0, 0.0, 0.0,
                            350.0,
                            10.0, 1.0e7,
                            nlev_c, p_surf, p_top,
                            Dict("N2" => 1.0), "";
                            real_gas=false,
                            thermo_functions=false,
                            flag_rayleigh=false,
                            benchmark_rt=true,
                            flag_cloud=false)
    ok || error("Failed to setup test atmosphere")
    atmosphere.allocate!(atmos, ""; check_safe_gas=false) || error("Failed to allocate test atmosphere")
    return atmos
end

@testset "diagnostics" begin
    atmos = _make_greygas_atmos()
    setpt.isothermal!(atmos, 400.0)
    atmosphere.calc_layer_props!(atmos)

    @testset "convective_variables" begin
        # test finding convective zone from mask
        fill!(atmos.mask_c, false)
        atmos.mask_c[4:9] .= true
        p_top, p_bot = diagnostics.estimate_convective_zone(atmos)
        @test p_top == atmos.pl[4]
        @test p_bot == atmos.pl[9]

        # test Rayleigh number estimation
        atmos.layer_kc .= 0.5
        atmos.layer_ρ .= 1.2
        atmos.layer_cp .= 1.0e3
        atmos.w_conv .= 2.0
        atmos.λ_conv .= 3.0
        diagnostics.estimate_Ra!(atmos)
        κ = phys.calc_therm_diffus(0.5, 1.2, 1.0e3)
        expected_Ra = (2.0 * 3.0 / κ)^(1.0 / phys.βRa)
        @test isapprox(atmos.diagnostic_Ra[1], expected_Ra; rtol=rtol)

        # test convective timescale estimation
        atmos.w_conv[1] = 1.0
        atmos.λ_conv[1] = 10.0
        diagnostics.estimate_timescale_conv!(atmos)
        @test isapprox(atmos.timescale_conv[1], 10.0; rtol=rtol)
        @test all(isfinite.(atmos.timescale_conv))
    end

    @testset "radiative_variables" begin
        # test radiative timescale estimation
        energy.calc_fluxes!(atmos, radiative=true)
        diagnostics.estimate_timescale_rad!(atmos)

        @test all(isfinite.(atmos.timescale_rad))
        @test all(atmos.timescale_rad .> 0.0)

        @test atmos.num_rt_eval == 2
        @test atmos.tim_rt_eval > 0.0
    end

    @testset "energy_utils_safe" begin
        # test handling NaN and Inf in energy fluxes
        arr = [1.0, Inf, NaN, -Inf]
        energy._make_finite!(arr, -5.0)
        @test isapprox(arr, [1.0, -5.0, -5.0, -5.0]; rtol=rtol) # converts to -5.0

        # test TKE-scheme exchange coefficient evaluation
        cd1 = energy.eval_exchange_coeff(2.0, 1.0e-2)
        cd2 = energy.eval_exchange_coeff(1.0e-8, 1.0e-2)
        @test isfinite(cd1) && cd1 > 0.0
        @test isfinite(cd2) && cd2 > 0.0

        # test calculating conductive-skin CBL flux
        atmos.tmp_magma = 1600.0
        atmos.tmp_surf = 1000.0
        atmos.skin_k = 3.0
        atmos.skin_d = 0.2
        fsk = energy.skin_flux(atmos)
        @test isapprox(fsk, 9000.0; rtol=rtol)

        # test calculating CBL skin depth
        @test isapprox(energy.skin_depth(atmos, fsk), 0.2; rtol=rtol)
        @test isapprox(energy.skin_depth(atmos, 1e20), 1.0e-6; rtol=rtol)
        @test isapprox(energy.skin_depth(atmos, 1e-20), 1.0e6; rtol=rtol)

        # sensible heat from TKE scheme, does not touch other terms
        atmos.flux_sens = 1.0
        atmos.is_out_lw = true
        atmos.is_out_sw = true
        fill!(atmos.flux_tot, 2.0)
        fill!(atmos.flux_cdct, 3.0)
        @test energy.reset_fluxes!(atmos)
        @test atmos.flux_sens == 0.0
        @test !atmos.is_out_lw
        @test !atmos.is_out_sw
        @test all(atmos.flux_tot .== 0.0)
        @test all(atmos.flux_cdct .== 0.0)

        # heating-rate helper
        atmos.flux_tot .= collect(0.0:(atmos.nlev_l-1))
        @test energy.calc_hrates!(atmos)
        expected_hr1 = (atmos.a[1] / atmos.layer_cp[1]) *
                        (atmos.flux_tot[2] - atmos.flux_tot[1]) /
                        (atmos.pl[2] - atmos.pl[1]) * 86400.0
        @test isapprox(atmos.heating_rate[1], expected_hr1; rtol=rtol)
    end
end

# Exobase diagnostic (src/state/diagnostics.jl: estimate_exobase!, _mfp_over_H).
@testset "exobase" begin

    function _isothermal_n2(T::Float64; p_top::Float64, nlev_c::Int64)
        atmos = _make_greygas_atmos(nlev_c=nlev_c, p_top=p_top)
        setpt.isothermal!(atmos, T)
        atmosphere.calc_layer_props!(atmos)
        return atmos
    end

    # A column reaching 1e-14 bar contains the N2 exobase (near 8e-12 bar for g = 10).
    # The exobase is reported at the first layer centre with l/H >= 1, so the analytic
    # crossing m g / (sqrt(2) σ), with g at the exobase radius, lies between that layer and
    # the one below it. p / g there does not depend on temperature, while the radius does.
    @testset "exobase_pressure_over_gravity_is_independent_of_temperature" begin
        res = Dict{Float64,Tuple{Float64,Float64}}()
        atmos = nothing
        for T in (300.0, 1500.0)
            atmos = _isothermal_n2(T; p_top=1e-14, nlev_c=100)
            @test diagnostics.estimate_exobase!(atmos)
            @test atmos.exobase_in_domain
            @test atmos.p[1] < atmos.exobase_p < atmos.p[end]
            @test isapprox(atmos.exobase_tmp, T; rtol=1e-9)

            # analytic crossing must sit between the exobase layer and the layer below it
            i     = findfirst(p -> isapprox(p, atmos.exobase_p; rtol=1e-12), atmos.p)
            gas   = atmos.gas_dat["N2"]
            σ     = pi * gas.particle_d^2
            g_exo = atmos.grav_surf * (atmos.rp / atmos.exobase_r)^2
            p_an  = gas.particle_m * g_exo / (sqrt(2.0) * σ)
            @test atmos.p[i] <= p_an < atmos.p[i+1]
            res[T] = (atmos.exobase_p / g_exo, atmos.exobase_r)
        end

        # log-uniform grid: p / g at the two temperatures may differ by at most one level
        i    = findfirst(p -> isapprox(p, atmos.exobase_p; rtol=1e-12), atmos.p)
        dlnp = log(atmos.p[i+1] / atmos.p[i])
        @test abs(log(res[300.0][1] / res[1500.0][1])) < dlnp
        @test res[1500.0][2] > 1.05 * res[300.0][2]   # hotter column is more extended

        # discrimination guard: a scale height taken with the surface gravity would put the
        # 1500 K crossing above the bracketing layers (gravity ratio ~1.3 > level spacing)
        gas = atmos.gas_dat["N2"]
        p_surf_g = gas.particle_m * atmos.grav_surf / (sqrt(2.0) * pi * gas.particle_d^2)
        @test !(atmos.p[i] <= p_surf_g < atmos.p[i+1])
    end

    # A column stopping at 1e-6 bar does not reach the exobase
    @testset "shallow_column_finds_exobase_above_top" begin
        deep    = _isothermal_n2(300.0; p_top=1e-14, nlev_c=100)
        shallow = _isothermal_n2(300.0; p_top=1e-6,  nlev_c=30)
        diagnostics.estimate_exobase!(deep)
        @test !diagnostics.estimate_exobase!(shallow)
        @test !shallow.exobase_in_domain
        @test shallow.exobase_p < shallow.p[2]
    end

    # Doubling every collision cross-section halves the exobase pressure; with no gas in
    # a layer the mean free path is unbounded rather than NaN.
    @testset "exobase_pressure_scales_inversely_with_cross_section" begin
        atmos = _isothermal_n2(300.0; p_top=1e-14, nlev_c=100)
        gas   = atmos.gas_dat["N2"]
        σ_ref = pi * gas.particle_d^2
        diagnostics.estimate_exobase!(atmos)
        i_ref = findfirst(p -> isapprox(p, atmos.exobase_p; rtol=1e-12), atmos.p)

        # doubled cross-section: the exact crossing (with σ_ref and g at the new exobase
        # radius) halves, and must lie between the new exobase layer and the one below it
        gas.particle_d *= sqrt(2.0)
        diagnostics.estimate_exobase!(atmos)
        i     = findfirst(p -> isapprox(p, atmos.exobase_p; rtol=1e-12), atmos.p)
        g_exo = atmos.grav_surf * (atmos.rp / atmos.exobase_r)^2
        p_an  = gas.particle_m * g_exo / (sqrt(2.0) * σ_ref)
        @test atmos.p[i] <= 0.5 * p_an < atmos.p[i+1]
        @test i < i_ref   # exobase moved to lower pressure

        # discrimination guard: a σ^(-1/2) law, or no σ dependence, would put the crossing
        # outside the bracketing layers (each factor > one level spacing of ~1.42)
        @test !(atmos.p[i] <= p_an / sqrt(2.0) < atmos.p[i+1])
        @test !(atmos.p[i] <= p_an < atmos.p[i+1])

        atmos.gas_vmr["N2"][1] = 0.0
        ratio, sigma, _ = diagnostics._mfp_over_H(atmos, 1)
        @test isinf(ratio)
        @test sigma < 1e-300
    end
end
