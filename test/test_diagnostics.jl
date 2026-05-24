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
        @test arr == [1.0, -5.0, -5.0, -5.0] # converts to -5.0

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
