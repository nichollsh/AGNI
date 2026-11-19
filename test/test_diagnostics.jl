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

    # Energy-transport timescale diagnostics (diagnostics.estimate_timescale_rad!,
    #   diagnostics.estimate_timescale_cond!). Invariants checked: the radiative estimate
    #   reproduces Showman, Cho & Menou (2010) Eq. 19; τ_rad ∝ T⁻³ and grows with Δp;
    #   τ_cond ∝ Δz²; τ_cond/τ_rad = 4σT³Δz/k by hydrostatic balance; the zero-conductivity
    #   limit returns +Inf rather than NaN or a negative time.
    @testset "timescales" begin
        σ_codata = 5.670374419e-8   # W m-2 K-4, written out so the pins do not reuse phys.σSB

        # The tests above overwrite cp, k and density with arbitrary values, so restore
        #   physically consistent layer properties for the isothermal 400 K N2 column
        setpt.isothermal!(atmos, 400.0)
        atmosphere.calc_layer_props!(atmos)
        pl_saved  = copy(atmos.pl)
        tmp_saved = copy(atmos.tmp)
        i_bot     = atmos.nlev_c

        @testset "conductive_to_radiative_ratio" begin
            # With Δp = ρ g Δz the heat capacity cancels in the ratio, leaving
            #   τ_cond/τ_rad = 4σT³Δz/k. This links the two diagnostics through the
            #   hydrostatic integrator rather than re-deriving either formula. The ρ, Δz
            #   and Δp/g from the integrator agree to better than 1e-3 in every layer
            diagnostics.estimate_timescale_rad!(atmos)
            diagnostics.estimate_timescale_cond!(atmos)
            ratio    = atmos.timescale_cond ./ atmos.timescale_rad
            expected = 4 .* σ_codata .* atmos.tmp.^3 .* atmos.layer_thick ./ atmos.layer_kc
            @test all(isapprox.(ratio, expected; rtol=1e-3))

            # Gas conduction is far slower than blackbody cooling in every layer, from the
            #   thin surface-refined bottom layer to the km-thick upper layers
            @test all(ratio .> 1e3)

            # Discrimination guard: a conductive time linear in Δz would shift the ratio by
            #   a factor Δz (16 m to 9 km on this grid), failing the pin in every layer
            ratio_wrong = (atmos.layer_ρ .* atmos.layer_cp .* atmos.layer_thick ./
                            atmos.layer_kc) ./ atmos.timescale_rad
            @test !any(isapprox.(ratio_wrong, expected; rtol=1e-3))
        end

        @testset "radiative_timescale_scales_mass" begin
            # In an isothermal column with constant cp, τ_rad ∝ Δp/g (the layer's mass per
            #   unit area), so τ_rad g/Δp is uniform even though Δp spans 7 decades
            diagnostics.estimate_timescale_rad!(atmos)
            dp     = atmos.pl[2:end] .- atmos.pl[1:end-1]
            scaled = atmos.timescale_rad .* atmos.g ./ dp
            @test all(isapprox.(scaled, scaled[1]; rtol=1e-12))
            @test maximum(dp) / minimum(dp) > 1e5
            @test all(isfinite.(atmos.timescale_rad)) && all(atmos.timescale_rad .> 0.0)

            # Edge case: the grid is refined towards the surface, so the bottom layer has
            #   higher p but smaller Δp than the one above. τ_rad must follow Δp, not p.
            @test atmos.p[end] > atmos.p[end-1]
            @test dp[end] < dp[end-1]
            @test atmos.timescale_rad[end] < atmos.timescale_rad[end-1]
        end

        @testset "radiative_timescale_scales_temperature" begin
            # Doubling T at fixed cp, Δp and g must shorten τ_rad by exactly 2³ = 8.
            #   Layer properties are not recomputed, so cp stays fixed.
            diagnostics.estimate_timescale_rad!(atmos)
            τ_ref = copy(atmos.timescale_rad)
            atmos.tmp .*= 2.0
            diagnostics.estimate_timescale_rad!(atmos)
            ratio = τ_ref ./ atmos.timescale_rad
            @test all(isapprox.(ratio, 8.0; rtol=1e-10))

            # Discrimination guard: a σT⁴-based formula would give a ratio of 16
            @test all(abs.(ratio .- 16.0) .> 1.0)
            atmos.tmp .= tmp_saved
        end

        @testset "radiative_timescale_reproduces_terrestrial" begin
            # Showman, Cho & Menou (2010) Eq. 19 and the text that follows it: a 0.2 bar
            #   layer at 300 K on an Earth-like planet (cp = 1000 J/kg/K, g = 9.81 m/s2)
            #   relaxes on a timescale comparable to an advection time of ~3 days.
            #   The bottom layer is used because it can hold Δp = 2e4 Pa below 10 bar.
            g_saved  = copy(atmos.g)
            cp_saved = copy(atmos.layer_cp)
            atmos.pl[i_bot]       = atmos.pl[i_bot+1] - 2.0e4
            atmos.g[i_bot]        = 9.81
            atmos.tmp[i_bot]      = 300.0
            atmos.layer_cp[i_bot] = 1000.0
            diagnostics.estimate_timescale_rad!(atmos)
            τ_rad = atmos.timescale_rad[i_bot]

            # Hand-calculated: 2e4/9.81 * 1000 / (4 σ 300³) = 3.3291e5 s = 3.85 days
            @test isapprox(τ_rad, 3.329090e5; rtol=1e-5)
            @test 1.0 * 86400.0 < τ_rad < 10.0 * 86400.0

            # Discrimination guard: dropping the factor 4 gives 1.33e6 s, and using σT⁴ in
            #   place of 4σT³ gives 4.4e3 s; both lie far outside the tolerance
            @test !isapprox(τ_rad, 1.331636e6; rtol=0.1)
            @test !isapprox(τ_rad, 4.438786e3; rtol=0.1)

            atmos.g        .= g_saved
            atmos.layer_cp .= cp_saved
            atmos.pl       .= pl_saved
            atmos.tmp      .= tmp_saved
        end

        @testset "conductive_timescale_matches_diffusion" begin
            # Air near 1 bar and 300 K (k = 0.026 W/m/K, ρ = 1.16 kg/m3, cp = 1005 J/kg/K)
            #   across a 100 m layer. Δz = 100 m is chosen so that Δz² and Δz differ by 100x.
            atmos.layer_kc[i_bot]    = 0.026
            atmos.layer_ρ[i_bot]     = 1.16
            atmos.layer_cp[i_bot]    = 1005.0
            atmos.layer_thick[i_bot] = 100.0
            diagnostics.estimate_timescale_cond!(atmos)
            τ_cond = atmos.timescale_cond[i_bot]

            # Hand-calculated: 1.16 * 1005 * 100² / 0.026 = 4.4838e8 s, about 14 years
            @test isapprox(τ_cond, 4.483846e8; rtol=1e-5)

            # Discrimination guard: a timescale linear in Δz would give 4.5e6 s
            @test !isapprox(τ_cond, 4.483846e6; rtol=0.1)

            # Doubling the thickness quadruples the diffusion time (linear would double it)
            atmos.layer_thick[i_bot] = 200.0
            diagnostics.estimate_timescale_cond!(atmos)
            @test isapprox(atmos.timescale_cond[i_bot] / τ_cond, 4.0; rtol=1e-10)
        end

        @testset "conductive_timescale_is_infinite" begin
            # Error contract: zero conductivity is the no-conduction limit and gives +Inf.
            #   Negative or NaN conductivity is invalid input and gives NaN, the "not
            #   computed" value, rather than a believable Inf or a negative time.
            atmosphere.calc_layer_props!(atmos)
            atmos.layer_kc[1] = 0.0
            atmos.layer_kc[2] = -1.0
            atmos.layer_kc[3] = NaN
            diagnostics.estimate_timescale_cond!(atmos)
            @test isinf(atmos.timescale_cond[1]) && atmos.timescale_cond[1] > 0.0
            @test isnan(atmos.timescale_cond[2])
            @test isnan(atmos.timescale_cond[3])

            # Layers with physical conductivity are unaffected: finite and positive
            @test all(isfinite.(atmos.timescale_cond[4:end]))
            @test all(atmos.timescale_cond[4:end] .> 0.0)
        end

        @testset "convective_timescale_is_zero" begin
            # Contract: τ_conv = λ/w is positive and finite at convective edges, and exactly
            #   zero where MLT finds the profile stable (w = λ = 0). The profile is
            #   superadiabatic (∇ = 0.4 > ∇_ad ≈ 0.29 for N2) below ~1.8 bar and isothermal
            #   at 200 K above, so both regimes are present in one column.
            prof(p) = max(400.0 * (p / atmos.pl[end])^0.4, 200.0)
            atmos.tmp  .= prof.(atmos.p)
            atmos.tmpl .= prof.(atmos.pl)
            atmosphere.calc_layer_props!(atmos)
            energy.convection!(atmos)
            diagnostics.estimate_timescale_conv!(atmos)
            conv = atmos.mask_c[1:atmos.nlev_c]
            @test any(conv) && !all(conv)   # guard: both regimes are exercised
            @test all(isfinite.(atmos.timescale_conv[conv]))
            @test all(atmos.timescale_conv[conv] .> 0.0)
            @test all(isapprox.(atmos.timescale_conv[.!conv], 0.0; atol=1e-300))

            # Edge case: a vanishing velocity with a non-zero mixing length is floored at
            #   eps(Float32), giving a large but finite time rather than Inf or NaN
            atmos.λ_conv[1] = 100.0
            atmos.w_conv[1] = 1e-12
            diagnostics.estimate_timescale_conv!(atmos)
            @test isapprox(atmos.timescale_conv[1], 100.0 / eps(Float32); rtol=1e-12)

            setpt.isothermal!(atmos, 400.0)
        end

        # Leave the shared atmosphere in a consistent state for the tests below
        atmosphere.calc_layer_props!(atmos)
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
