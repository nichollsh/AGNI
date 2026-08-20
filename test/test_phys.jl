using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
RES_DIR         = joinpath(ROOT_DIR,"res/")
OUT_DIR         = joinpath(ROOT_DIR,"out/")
TEST_DIR        = joinpath(ROOT_DIR,"test/")
const lookup_mmw = AGNI.formulae._lookup_mmw
const lookup_count_atoms = AGNI.formulae._lookup_count_atoms
const lookup_colour = AGNI.style._lookup_colour
const lookup_liquid_rho = AGNI.density._lookup_liquid_rho

@testset "phys" begin
    # atom counting
    @test AGNI.formulae.count_atoms("H2O") == lookup_count_atoms["H2O"]
    @test AGNI.formulae.count_atoms("CO2") == lookup_count_atoms["CO2"]

    # all molecules
    for molec in AGNI.consts.vols_standard
        @test length(AGNI.formulae.count_atoms(molec)) > 0
        @test AGNI.formulae.get_mmw(molec) > 0.0
    end

    # same_atoms
    @test AGNI.formulae.same_atoms(Dict("H"=>2, "O"=>1), Dict("O"=>1, "H"=>2))

    # mean molecular weight
    @test isapprox(AGNI.formulae.get_mmw("H2O"), lookup_mmw["H2O"]; rtol=1e-12)

    # pretty name replaces digits (subscript unicode); ensure result differs
    pn = AGNI.style.pretty_name("H2O")
    @test pn != "H2O"
    @test !occursin("2", pn)

    # pretty colour for known gas
    @test AGNI.style.pretty_colour("H2O") == lookup_colour["H2O"]

    # ideal density positive
    rho = AGNI.density._rho_ideal(300.0, 1e5, lookup_mmw["CO2"])
    @test isfinite(rho) && rho > 0.0

    # planck positive
    p = AGNI.phys.evaluate_planck(500.0, 300.0)
    @test p > 0.0

    # gravity approximate
    g = AGNI.phys.grav_accel(5.972e24, 6.371e6)
    @test isapprox(g, 9.81; atol=1.0)

    # liquid density table and fallback
    @test AGNI.density.liquid_rho("H2O") == lookup_liquid_rho["H2O"]
    @test AGNI.density.liquid_rho("UNKNOWN") == AGNI.consts.BIGFLOAT

    # thermal diffusivity
    α = AGNI.phys.calc_therm_diffus(0.6, 1000.0, 1000.0)
    @test isapprox(α, 6e-7; rtol=1e-12)

    # equilibrium and skin temperatures
    Teq = AGNI.phys.calc_Teq(1361.0, 0.3)
    @test Teq > 0.0
    Tskin = AGNI.phys.calc_Tskin(1361.0, 0.3)
    @test isapprox(Tskin, Teq * (0.5^0.25); rtol=1e-12)


    # -------------
    # Test thermodynamic lookup data validity
    # -------------
    ideal_H2O::species.Gas_t = species.load_gas("$RES_DIR/thermodynamics/", "H2O", true, false)
    @test !ideal_H2O.fail


    # -------------
    # Test heat capacity lookup tables
    # -------------
    @testset "cp" begin
        t_test  = [10.0,  500.0, 1000.0, 2000.0, 3000.0]     # Tested values of temperature
        v_expt  = [4.975, 35.22, 41.27 , 51.20 , 55.74 ]     # Expected values of cp [J mol-1 K-1]
        v_obs   = zero(t_test)
        test_pass = true
        for i in 1:5
            v_obs[i] = species.get_Cp(ideal_H2O, t_test[i]) * ideal_H2O.mmw # get value and convert units
            test_pass &= isapprox(v_expt[i], v_obs[i]; rtol=1e-3)
        end
        if !test_pass
            @error "Expected values = $(v_expt) J mol-1 K-1\n Modelled values = $(v_obs) J mol-1 K-1"
        end
        @test test_pass
    end

    # -------------
    # Test ideal gas equation of state
    # -------------
    @testset "ideal_EOS" begin
        t_test = [200.0,  300.0, 500.0,   1273.0,  3200.0] # Tested values of temperature [K]
        p_test = [1e0,    1e3,   1e5,     1e7,     1e8]    # Tested values of pressure [Pa]
        v_expt = [1.0833532e-5, 7.2223549e-3, 4.33341295e-1, 1.7020475e1, 6.7709577e1]  # Expected rho [kg m-3]
        v_obs  = zero(p_test)
        test_pass = true
        for i in 1:5
            v_obs[i] = density.calc_rho_gas(t_test[i], p_test[i], ideal_H2O)
            test_pass &= isapprox(v_expt[i], v_obs[i]; rtol=1e-3, atol=1e-12)
        end

        if !test_pass
            @error "Expected values = $(v_expt) kg m-3\n Modelled values = $(v_obs) kg m-3"
        end
        @test test_pass
    end

    # -------------
    # Test AQUA equation of state (phs_method=2)
    # This method switches to ideal gas in the condensed region
    # -------------
    @testset "AQUA_EOS_PHSMETHOD2" begin
        aqua_H2O::species.Gas_t = species.load_gas("$RES_DIR/thermodynamics/", "H2O", true, true)
        t_test = [100.0,  200.0, 500.0,   1273.0,  4000.0] # Tested values of temperature [K]
        p_test = [1e1,    1e3,   1e5,     1e7,     1e8]    # Tested values of pressure [Pa]
        v_expt = [0.00021667064761346432, 0.010833532380673217, 0.43521435354778626, 17.02216182825162, 51.22864662875654]
        v_obs  = zero(p_test)
        test_pass = true
        for i in 1:5
            v_obs[i] = density.calc_rho_gas(t_test[i], p_test[i], aqua_H2O; phs_method=2)
            test_pass &= isapprox(v_expt[i], v_obs[i]; rtol=1e-3)
        end

        if !test_pass
            @error "Expected values = $(v_expt) kg m-3\n Modelled values = $(v_obs) kg m-3"
        end
        @test test_pass
    end


    # -------------
    # Test AQUA equation of state (phs_method=3)
    # This method uses a shifted log10pressure to evaluate the density
    # -------------
    @testset "AQUA_EOS_PHSMETHOD3" begin
        aqua_H2O::species.Gas_t = species.load_gas("$RES_DIR/thermodynamics/", "H2O", true, true)
        t_test = [200.0,  300.0, 500.0,   1273.0,  3200.0] # Tested values of temperature [K]
        p_test = [1e0,    1e3,   1e5,     1e7,     1e8]    # Tested values of pressure [Pa]
        v_expt = [926.1211619878637, 0.007227415287350509, 0.43521435354778626, 17.02216182825162, 66.82013418840567]
        v_obs  = zero(p_test)
        test_pass = true
        for i in 1:5
            v_obs[i] = density.calc_rho_gas(t_test[i], p_test[i], aqua_H2O; phs_method=3)
            test_pass &= isapprox(v_expt[i], v_obs[i]; rtol=1e-3)
        end

        if !test_pass
            @error "Expected values = $(v_expt) kg m-3\n Modelled values = $(v_obs) kg m-3"
        end
        @test test_pass
    end

    # -------------
    # Test AQUA equation of state (phs_method=4)
    # This method uses a scaled density relative to some shifted-pressure evaluation.
    # -------------
    @testset "AQUA_EOS_PHSMETHOD4" begin
        aqua_H2O::species.Gas_t = species.load_gas("$RES_DIR/thermodynamics/", "H2O", true, true)
        t_test = [200.0,  300.0, 500.0,   1273.0,  3200.0] # Tested values of temperature [K]
        p_test = [1e0,    1e3,   1e5,     1e7,     1e8]    # Tested values of pressure [Pa]

        # index 1 (200 K, 1 Pa) is below the table floor once shifted by phs_dlogp, so
        # the correct value is the analytic ideal-gas density, not a table lookup
        v_expt = [density._rho_ideal(200.0, 1e0, aqua_H2O.mmw),
                    0.007227415287350509, 0.43521435354778626, 17.02216182825162, 66.82013418840567]
        v_obs  = zero(p_test)
        test_pass = true
        for i in 1:5
            v_obs[i] = density.calc_rho_gas(t_test[i], p_test[i], aqua_H2O; phs_method=4)
            test_pass &= isapprox(v_expt[i], v_obs[i]; rtol=1e-3)
        end

        if !test_pass
            @error "Expected values = $(v_expt) kg m-3\n Modelled values = $(v_obs) kg m-3"
        end
        @test test_pass

        # test that a broken fallback returns the ice density of ~926 kg/m3
        # at this point instead of the ~1e-5 ideal-gas value
        wrong_condensed_fallback = 926.1211619878637
        @test !isapprox(v_obs[1], wrong_condensed_fallback; rtol=1e-3)
    end


    # -------------
    # Test VdW equation of state
    # -------------
    @testset "vdw_EOS" begin
        vdw_CO2::species.Gas_t = species.load_gas("$RES_DIR/thermodynamics/", "CO2", true, true)
        t_test = [200.0,  300.0, 500.0,   1273.0,  3200.0] # Tested values of temperature [K]
        p_test = [1e0,    1e3,   1e5,     1e7,     1e8]    # Tested values of pressure [Pa]
        v_expt = [2.6465333190669653e-5, 0.017644297244863345, 1.0597595155879878, 41.19090754552491, 147.5212980424141]
        v_obs  = zero(p_test)
        test_pass = true
        for i in 1:5
            v_obs[i] = density.calc_rho_gas(t_test[i], p_test[i], vdw_CO2)
            test_pass &= isapprox(v_expt[i], v_obs[i]; rtol=1e-3)
        end
        if !test_pass
            @error "Expected values = $(v_expt) kg m-3\n Modelled values = $(v_obs) kg m-3"
        end
        @test test_pass
    end


    # -------------
    # Test that calc_rho_gas's vapour/condensed branch test uses the evaluation
    # pressure (total pressure, under Amagat's law) rather than a species' partial
    # pressure.
    # -------------
    @testset "calc_rho_gas_branch_test_uses_total_not_partial_pressure" begin
        aqua_H2O::species.Gas_t = species.load_gas("$RES_DIR/thermodynamics/", "H2O", true, true)
        tmp::Float64 = 540.0
        vmr_h2o::Float64 = 0.3
        phs_dlogp::Float64 = 0.3

        # column equals exactly at H2O's own saturation pressure (total pressure),
        # but H2O's partial pressure (vmr_h2o * prs) is only 30% of that
        prs_total::Float64 = 10.0 ^ aqua_H2O.sat_I(tmp)

        # branch test evaluated at the total (evaluation) pressure: correctly
        # identifies this point as needing phase-boundary handling
        @test species.is_vapour(aqua_H2O, tmp, prs_total; phs_εlogp=-phs_dlogp) == false

        # branch test evaluated at the partial pressure: misclassifies this point as
        # vapour, since the partial pressure is sub-saturated relative to phs_dlogp
        @test species.is_vapour(aqua_H2O, tmp, prs_total*vmr_h2o; phs_εlogp=-phs_dlogp) == true

        # calc_rho_gas (fixed) must use the total-pressure branch test, and so must
        # NOT reduce to the naive evaluation at this point
        rho_fixed::Float64 = density.calc_rho_gas(tmp, prs_total, aqua_H2O;
                                                    phs_method=4, phs_dlogp=phs_dlogp)
        naive_eval::Float64 = 10.0 ^ aqua_H2O.eos_I(tmp, log10(prs_total))
        @test isfinite(rho_fixed) && rho_fixed > 0.0
        @test isapprox(rho_fixed, 23.095477724934334; rtol=1e-3)

        # naive evaluation at p=psat inflates the density by roughly an order of
        # magnitude relative phase-boundary-aware evaluation
        @test naive_eval > 5.0 * rho_fixed
    end


    # -------------
    # Test that calc_rho_mix's ngas==1 fast path forwards phs_dlogp through to
    # calc_rho_gas, rather than always using PHS_DLOGP_DEFAULT regardless of the
    # caller's choice.
    # -------------
    @testset "calc_rho_mix_single_gas_forwards_phs_dlogp" begin
        aqua_H2O::species.Gas_t = species.load_gas("$RES_DIR/thermodynamics/", "H2O", true, true)
        tmp::Float64 = 270.0
        prs::Float64 = 1e5

        # calc_rho_mix with a single gas must delegate exactly to calc_rho_gas with
        # the same (non-default) phs_dlogp value
        rho_mix_dlogp01 = density.calc_rho_mix([aqua_H2O], [1.0], tmp, prs, aqua_H2O.mmw;
                                                    phs_method=4, phs_dlogp=0.1)
        rho_gas_dlogp01 = density.calc_rho_gas(tmp, prs, aqua_H2O; phs_method=4, phs_dlogp=0.1)
        @test isapprox(rho_mix_dlogp01, rho_gas_dlogp01; rtol=1e-10)

        # if phs_dlogp were silently ignored, the two values below would match instead
        rho_gas_default_dlogp = density.calc_rho_gas(tmp, prs, aqua_H2O;
                                                        phs_method=4, phs_dlogp=density.PHS_DLOGP_DEFAULT)
        @test !isapprox(rho_mix_dlogp01, rho_gas_default_dlogp; rtol=1e-6)
    end


    # -------------
    # Test calc_rho_mix against the H2O+H2 scenario to verify Amagat's law directly for
    # every phs_method, and checks the low-pressure limit where all
    # methods must agree and approach the ideal-gas mixture density.
    # -------------
    @testset "calc_rho_mix_amagat_law_holds_across_phs_methods" begin
        aqua_H2O::species.Gas_t = species.load_gas("$RES_DIR/thermodynamics/", "H2O", true, true)
        cms19_H2::species.Gas_t = species.load_gas("$RES_DIR/thermodynamics/", "H2",  true, true)
        tmp::Float64 = 270.0
        prs::Float64 = 1e5
        vmr::Array{Float64,1} = [0.7, 0.3]
        mmw::Float64 = vmr[1]*aqua_H2O.mmw + vmr[2]*cms19_H2.mmw
        mmr_h2o::Float64 = vmr[1]*aqua_H2O.mmw/mmw
        mmr_h2::Float64  = vmr[2]*cms19_H2.mmw/mmw

        for m in 1:4
            rho_mix = density.calc_rho_mix([aqua_H2O, cms19_H2], vmr, tmp, prs, mmw;
                                                phs_method=m, phs_dlogp=0.3)
            rho_h2o = density.calc_rho_gas(tmp, prs, aqua_H2O; phs_method=m, phs_dlogp=0.3)
            rho_h2  = density.calc_rho_gas(tmp, prs, cms19_H2; phs_method=m, phs_dlogp=0.3)

            # Amagat's law: 1/rho_mix = sum(mmr_i / rho_i); check this independently
            # of the internals of calc_rho_mix, using each component's own density
            rho_amagat = 1.0 / (mmr_h2o/rho_h2o + mmr_h2/rho_h2)
            @test isfinite(rho_mix) && rho_mix > 0.0
            @test isapprox(rho_mix, rho_amagat; rtol=1e-9)
        end

        # end member case at 1 Pa, where both components are far below their saturation
        # pressures at 270 K. All phs_method agree and the result must be close to the
        # ideal-gas mixture in this regime.
        prs_lo::Float64 = 1.0
        rho_mix_lo = [density.calc_rho_mix([aqua_H2O, cms19_H2], vmr, tmp, prs_lo, mmw;
                                                phs_method=m, phs_dlogp=0.3) for m in 1:4]
        for m in 2:4
            @test isapprox(rho_mix_lo[m], rho_mix_lo[1]; rtol=1e-9)
        end
        rho_ideal_h2o = density._rho_ideal(tmp, prs_lo, aqua_H2O.mmw)
        rho_ideal_h2  = density._rho_ideal(tmp, prs_lo, cms19_H2.mmw)
        rho_ideal_mix = 1.0 / (mmr_h2o/rho_ideal_h2o + mmr_h2/rho_ideal_h2)
        @test isapprox(rho_mix_lo[1], rho_ideal_mix; rtol=0.02)
    end


    # -------------
    # Test mixing ratios
    # -------------
    @testset "mixing_ratios" begin
        tmp_surf        = 200.0     # Surface temperature [kelvin]
        toa_heating     = 10000.00  # Instellation flux [W m-2]
        p_surf          = 1.0       # bar
        p_top           = 1e-8
        theta           = 65.0
        gravity         = 10.0
        nlev_centre     = 100
        radius          = 1.0e7    # metres
        mf_dict         = Dict([
                                ("H2O" , 0.5),
                                ("CO2" , 0.2),
                                ("N2"  , 0.1),
                                ("H2"  , 0.2)
                                ])
        spfile_name   ="$RES_DIR/spectral_files/Dayspring/48/Dayspring.sf"

        # Setup atmosphere
        atmos = atmosphere.Atmos_t()
        atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                                spfile_name,
                                toa_heating, 1.0, 0.0, theta,
                                tmp_surf,
                                gravity, radius,
                                nlev_centre, p_surf, p_top,
                                mf_dict, ""
                        )
        atmosphere.allocate!(atmos,"")

        dct_e::Dict{String, Float64} = mf_dict
        dct_o::Dict{String, Float64} = Dict()
        test_pass = true
        for k in keys(dct_e)
            dct_o[k] = atmos.gas_vmr[k][20]
            test_pass &= isapprox(dct_e[k], dct_o[k]; atol=1e-6)
        end

        if !test_pass
            @error "Expected values = $(dct_e)\n Modelled values = $(dct_o)"
        end
        atmosphere.deallocate!(atmos)
        @test test_pass
    end


    # -------------
    # Test saturation pressure lookup
    # At the boiling point of water (373.15 K) Psat should be ~1 atm = 101325 Pa
    # -------------
    @testset "Psat" begin
        gas_H2O::species.Gas_t = species.load_gas("$RES_DIR/thermodynamics/", "H2O", true, false)
        T_boil  = 373.15  # [K]
        Psat_e  = 101325.0  # expected ~1 atm [Pa]
        Psat_o  = species.get_Psat(gas_H2O, T_boil)
        @test isfinite(Psat_o)
        @test Psat_o > 0.0
        @test isapprox(Psat_o, Psat_e; rtol=0.05)
    end


    # -------------
    # Test latent heat lookup
    # At the boiling point of water (373.15 K) Lv ≈ 2.256e6 J/kg
    # -------------
    @testset "Lv" begin
        gas_H2O::species.Gas_t = species.load_gas("$RES_DIR/thermodynamics/", "H2O", true, false)
        T_boil  = 373.15   # [K]
        Lv_e    = 2.256e6  # expected latent heat [J/kg]
        Lv_o    = species.get_Lv(gas_H2O, T_boil)
        @test isfinite(Lv_o)
        @test Lv_o > 0.0
        @test isapprox(Lv_o, Lv_e; rtol=0.05)
    end


    # -------------
    # Test thermal conductivity (kinetic theory estimate)
    # Must be positive and in a physically reasonable range for water vapour
    # -------------
    @testset "Kc" begin
        gas_H2O::species.Gas_t = species.load_gas("$RES_DIR/thermodynamics/", "H2O", true, false)
        t_test  = [200.0, 500.0, 1000.0]  # [K]
        for t in t_test
            Kc = species.get_Kc(gas_H2O, t)
            @test isfinite(Kc) && Kc > 0.0
        end
        # conductivity must increase with temperature (∝ sqrt(T))
        @test species.get_Kc(gas_H2O, 500.0) > species.get_Kc(gas_H2O, 200.0)
    end


    # -------------
    # Test the H2O demixing-temperature fit (_Tdemix_H2O).
    # This is a pure analytic formula from Appendix A of the cited paper.
    # -------------
    @testset "Tdemix_H2O_fit_is_symmetric_about_peak_composition" begin
        d::Float64 = 0.4498  # peak molar fraction (Table A1 coefficient)
        for p in [1e7, 1e8, 5e8, 1e9]  # Pa, spanning ~0.1-10 kbar
            t_peak::Float64  = species._Tdemix_H2O(p, d)
            t_left::Float64  = species._Tdemix_H2O(p, d - 0.05)
            t_right::Float64 = species._Tdemix_H2O(p, d + 0.05)

            # exact symmetry of the Lorentzian term about x=d
            @test isapprox(t_left, t_right; rtol=1e-12)

            # the fit's Lorentzian term peaks at x=d for these (physically positive)
            # coefficients, so the demixing temperature there exceeds nearby values
            @test t_peak > t_left
        end

        # away from the peak the demixing temperature is substantially lower than peak
        t_peak_1e8::Float64 = species._Tdemix_H2O(1e8, d)
        t_far_1e8::Float64  = species._Tdemix_H2O(1e8, 0.05)
        @test !isapprox(t_peak_1e8, t_far_1e8; rtol=0.1)
        @test t_far_1e8 < t_peak_1e8

        # edge case: x at the boundaries of its physical domain (pure H2 / pure H2O)
        @test isfinite(species._Tdemix_H2O(1e8, 0.0))
        @test isfinite(species._Tdemix_H2O(1e8, 1.0))
    end


    # -------------
    # Test that the demixing temperature increases with pressure at fixed
    # composition, consistent with immiscibility being enhanced at higher pressure.
    # -------------
    @testset "Tdemix_H2O_fit_increases_with_pressure_at_fixed_composition" begin
        d::Float64 = 0.4498
        p_test::Array{Float64,1} = [1e7, 1e8, 5e8, 1e9]  # Pa
        t_test::Array{Float64,1} = [species._Tdemix_H2O(p, d) for p in p_test]
        @test all(isfinite.(t_test))
        # strictly increasing (not just non-decreasing) across this pressure range
        @test all(diff(t_test) .> 0.0)
    end


    # -------------
    # Test get_Tdemix, the gas-aware wrapper that dispatches to a species-specific
    # demixing fit (only implemented for H2O; other species have no known fit).
    # -------------
    @testset "get_Tdemix_dispatches_by_species" begin
        gas_H2O::species.Gas_t = species.load_gas("$RES_DIR/thermodynamics/", "H2O", true, false)
        gas_CO2::species.Gas_t = species.load_gas("$RES_DIR/thermodynamics/", "CO2", true, false)
        prs::Float64 = 1e8
        x::Float64   = 0.4

        # H2O uses the fitted demixing curve directly
        @test isapprox(species.get_Tdemix(gas_H2O, prs, x), species._Tdemix_H2O(prs, x); rtol=1e-12)

        # species without a demixing return a value far below any physical
        # temperature, so that "atmos.tmp < get_Tdemix(...)" never triggers
        tdemix_co2::Float64 = species.get_Tdemix(gas_CO2, prs, x)
        @test isapprox(tdemix_co2, -1.0 * AGNI.consts.BIGFLOAT; rtol=1e-12)
        @test tdemix_co2 < -1.0e5
    end


    # -------------
    # Test is_vapour
    # Ideal gas is always in the vapour phase (no real-gas condensation)
    # -------------
    @testset "is_vapour" begin
        gas_Ne::species.Gas_t = species.load_gas("$RES_DIR/thermodynamics/", "Ne", true, false)
        # ideal / no-sat stub is always vapour
        @test species.is_vapour(gas_Ne, 300.0, 1e5)
        @test species.is_vapour(gas_Ne, 100.0, 1e8)
    end


    # -------------
    # Test is_vapour()'s phs_εlogp offset against real AQUA saturation curve.
    # With phs_εlogp=-phs_dlogp, the switch away from "vapour" occurs phs_dlogp
    # below the true saturation pressure (not at psat itself).
    # -------------
    @testset "is_vapour_switch_boundary_is_offset_below_saturation_curve" begin
        gas_H2O::species.Gas_t = species.load_gas("$RES_DIR/thermodynamics/", "H2O", true, true)
        T_boil::Float64 = 373.15
        phs_dlogp::Float64 = 0.3
        psat::Float64 = 10.0 ^ gas_H2O.sat_I(T_boil)

        # default (near-zero) tolerance: switch essentially at psat itself
        @test species.is_vapour(gas_H2O, T_boil, psat*0.999)
        @test !species.is_vapour(gas_H2O, T_boil, psat*1.001)

        # with phs_dlogp=0.3, the switch is shifted to 0.3 dex below psat: still
        # vapour comfortably inside that band, no longer vapour once inside it
        @test species.is_vapour(gas_H2O, T_boil, psat*10.0^(-0.35); phs_εlogp=-phs_dlogp)
        @test !species.is_vapour(gas_H2O, T_boil, psat*10.0^(-0.25); phs_εlogp=-phs_dlogp)

        # pressures above psat are never classified as vapour, regardless of offset
        @test !species.is_vapour(gas_H2O, T_boil, psat*1.5; phs_εlogp=-phs_dlogp)
    end


    # -------------
    # Test _pretty_colour
    # Known gases return their lookup colour; unknown gases return a valid hex string
    # -------------
    @testset "pretty_colour" begin
        # known lookup entry
        @test style.pretty_colour("CO2") == lookup_colour["CO2"]
        # unknown molecule: must return a 7-character hex code starting with '#'
        col_sio = style.pretty_colour("SiO")
        @test length(col_sio) == 7
        @test col_sio[1] == '#'
    end


    # -------------
    # Test count_atoms with various formulas
    # -------------
    @testset "count_atoms_extended" begin
        # Simple molecules
        @test formulae.count_atoms("O2") == Dict("O" => 2)
        @test formulae.count_atoms("N2") == Dict("N" => 2)

        # Molecules with parentheses/brackets (should be skipped by parser)
        atoms_nh3 = formulae.count_atoms("NH3")
        @test atoms_nh3["N"] == 1
        @test atoms_nh3["H"] == 3

        # Two-letter element names
        atoms_ch4 = formulae.count_atoms("CH4")
        @test atoms_ch4["C"] == 1
        @test atoms_ch4["H"] == 4

        # Complex molecule
        atoms_h2so4 = formulae.count_atoms("H2SO4")
        @test atoms_h2so4["H"] == 2
        @test atoms_h2so4["S"] == 1
        @test atoms_h2so4["O"] == 4

        # Multi-digit stoichiometry
        atoms_glucose = formulae.count_atoms("C6H12O6")
        @test atoms_glucose == Dict("C" => 6, "H" => 12, "O" => 6)

        atoms_sucrose = formulae.count_atoms("C12H22O11")
        @test atoms_sucrose == Dict("C" => 12, "H" => 22, "O" => 11)

        atoms_c60 = formulae.count_atoms("C60")
        @test atoms_c60 == Dict("C" => 60)
    end


    # -------------
    # Test same_atoms function
    # -------------
    @testset "same_atoms" begin
        @test formulae.same_atoms(Dict("H"=>2, "O"=>1), Dict("O"=>1, "H"=>2)) == true
        @test formulae.same_atoms(Dict("C"=>1, "O"=>2), Dict("C"=>1, "O"=>2)) == true
        @test formulae.same_atoms(Dict("H"=>2, "O"=>1), Dict("H"=>2, "O"=>2)) == false
        @test formulae.same_atoms(Dict("H"=>2), Dict("H"=>2, "O"=>1)) == false
    end


    # -------------
    # Test _get_mmw for various molecules
    # -------------
    @testset "get_mmw" begin
        # Test known molecules
        @test isapprox(formulae.get_mmw("H2O"), lookup_mmw["H2O"]; rtol=1e-12)
        @test isapprox(formulae.get_mmw("CO2"), lookup_mmw["CO2"]; rtol=1e-12)
        @test isapprox(formulae.get_mmw("N2"), lookup_mmw["N2"]; rtol=1e-12)

        # Test that MMW is positive
        @test formulae.get_mmw("CH4") > 0.0
        @test formulae.get_mmw("O2") > 0.0
    end


    # -------------
    # Test _pretty_name generates subscripted text
    # -------------
    @testset "pretty_name" begin
        # Function replaces digits with unicode subscripts
        pn_h2o = style.pretty_name("H2O")
        @test pn_h2o != "H2O"  # should be different due to subscript
        @test !occursin("2", pn_h2o)  # regular '2' should be gone

        pn_co2 = style.pretty_name("CO2")
        @test pn_co2 != "CO2"
        @test !occursin("2", pn_co2)

        # Molecule without numbers stays same (except unicode conversion)
        pn_o = style.pretty_name("O")
        @test length(pn_o) >= 1
    end


    # -------------
    # Test error handling in load_gas for unsupported elements
    # -------------
    @testset "load_gas_unsupported_element" begin
        # Try to load a gas with an unsupported element (should error)
        # suppress error output for cleaner test logs
        with_logger(MinLevelLogger(current_logger(), Test.Logging.Error+1)) do
            gas = species.load_gas("$RES_DIR/thermodynamics/", "Xz2", true, false)
            @test gas.fail == true
        end
    end


    # -------------
    # Test load_gas's handling of malformed thermodynamic data files.
    # Each scenario writes a minimal synthetic NetCDF file to a scratch directory.
    # -------------
    @testset "load_gas_rejects_malformed_data_files" begin
        scratch_dir::String = mktempdir()

        # -- corrupted file: the companion .chk hash does not match the file --
        @testset "corrupt_checksum" begin
            formula = "H2O"
            fpath = joinpath(scratch_dir, "$formula.nc")
            write(fpath, "this is not a real netcdf file")
            write(fpath * ".chk", "0"^64)  # deliberately wrong hash

            species.ENABLE_CHECKSUM = true
            local gas::species.Gas_t
            with_logger(MinLevelLogger(current_logger(), Test.Logging.Error+1)) do
                gas = species.load_gas(scratch_dir, formula, true, false; check_integrity=true)
            end
            @test gas.fail == true
            @test gas.eos == species.EOS_IDEAL

            rm(fpath); rm(fpath * ".chk")
        end

        # -- missing creation date --
        @testset "missing_creation_date" begin
            formula = "H2O"
            fpath = joinpath(scratch_dir, "$formula.nc")
            NCDataset(fpath, "c") do ds
                defVar(ds, "mmw", 0.018, ())
            end

            local gas::species.Gas_t
            with_logger(MinLevelLogger(current_logger(), Test.Logging.Error+1)) do
                gas = species.load_gas(scratch_dir, formula, true, false; check_integrity=false)
            end
            @test gas.fail == true
            @test gas.eos == species.EOS_IDEAL

            rm(fpath)
        end

        # -- outdated creation date --
        @testset "outdated_creation_date" begin
            formula = "H2O"
            fpath = joinpath(scratch_dir, "$formula.nc")
            created_fixture::Int64 = 20200101  # predates MIN_DATA_VERSION
            NCDataset(fpath, "c") do ds
                defVar(ds, "created", created_fixture, ())
                defVar(ds, "mmw", 0.018, ())
            end

            # sanity check on the fixture itself: it must actually be outdated
            @test created_fixture < species.MIN_DATA_VERSION

            local gas::species.Gas_t
            with_logger(MinLevelLogger(current_logger(), Test.Logging.Error+1)) do
                gas = species.load_gas(scratch_dir, formula, true, false; check_integrity=false)
            end
            @test gas.fail == true
            @test gas.eos == species.EOS_IDEAL

            rm(fpath)
        end

        # -- EOS pressure axis not strictly ascending --
        @testset "eos_pressure_axis_not_ascending" begin
            formula = "CO2"
            fpath = joinpath(scratch_dir, "$formula.nc")
            NCDataset(fpath, "c") do ds
                defVar(ds, "created", species.MIN_DATA_VERSION, ())
                defVar(ds, "mmw", 0.044, ())
                defVar(ds, "JANAF", "CO2", ())
                defDim(ds, "vdw_T", 2)
                defDim(ds, "vdw_P", 3)
                defVar(ds, "vdw_T", [100.0, 200.0], ("vdw_T",))
                defVar(ds, "vdw_P", [10.0, 5.0, 20.0], ("vdw_P",))  # not ascending
                defVar(ds, "vdw_rho", zeros(2,3), ("vdw_T","vdw_P"))
            end

            local gas::species.Gas_t
            with_logger(MinLevelLogger(current_logger(), Test.Logging.Error+1)) do
                gas = species.load_gas(scratch_dir, formula, true, true; check_integrity=false)
            end
            @test gas.fail == true
            @test isapprox(gas.mmw, 0.044; rtol=1e-12)

            rm(fpath)
        end

        # -- EOS temperature axis not strictly ascending --
        @testset "eos_temperature_axis_not_ascending" begin
            formula = "CO2"
            fpath = joinpath(scratch_dir, "$formula.nc")
            NCDataset(fpath, "c") do ds
                defVar(ds, "created", species.MIN_DATA_VERSION, ())
                defVar(ds, "mmw", 0.044, ())
                defVar(ds, "JANAF", "CO2", ())
                defDim(ds, "vdw_T", 2)
                defDim(ds, "vdw_P", 3)
                defVar(ds, "vdw_T", [200.0, 100.0], ("vdw_T",))  # not ascending
                defVar(ds, "vdw_P", [10.0, 15.0, 20.0], ("vdw_P",))
                defVar(ds, "vdw_rho", zeros(2,3), ("vdw_T","vdw_P"))
            end

            local gas::species.Gas_t
            with_logger(MinLevelLogger(current_logger(), Test.Logging.Error+1)) do
                gas = species.load_gas(scratch_dir, formula, true, true; check_integrity=false)
            end
            @test gas.fail == true
            @test isapprox(gas.mmw, 0.044; rtol=1e-12)

            rm(fpath)
        end

        rm(scratch_dir; recursive=true, force=true)
    end

end
