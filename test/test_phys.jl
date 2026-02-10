using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
RES_DIR         = joinpath(ROOT_DIR,"res/")
OUT_DIR         = joinpath(ROOT_DIR,"out/")
TEST_DIR        = joinpath(ROOT_DIR,"test/")

@testset "phys" begin
    # atom counting
    @test AGNI.phys.count_atoms("H2O") == AGNI.consts._lookup_count_atoms["H2O"]
    @test AGNI.phys.count_atoms("CO2") == AGNI.consts._lookup_count_atoms["CO2"]

    # same_atoms
    @test AGNI.phys.same_atoms(Dict("H"=>2, "O"=>1), Dict("O"=>1, "H"=>2))

    # mean molecular weight
    @test isapprox(AGNI.phys._get_mmw("H2O"), AGNI.consts._lookup_mmw["H2O"]; rtol=1e-12)

    # pretty name replaces digits (subscript unicode); ensure result differs
    pn = AGNI.phys._pretty_name("H2O")
    @test pn != "H2O"
    @test !occursin("2", pn)

    # pretty colour for known gas
    @test AGNI.phys._pretty_color("H2O") == AGNI.consts._lookup_color["H2O"]

    # ideal density positive
    rho = AGNI.phys._rho_ideal(300.0, 1e5, AGNI.consts._lookup_mmw["CO2"])
    @test isfinite(rho) && rho > 0.0

    # planck positive
    p = AGNI.phys.evaluate_planck(500.0, 300.0)
    @test p > 0.0

    # gravity approximate
    g = AGNI.phys.grav_accel(5.972e24, 6.371e6)
    @test isapprox(g, 9.81; atol=1.0)

    # liquid density table and fallback
    @test AGNI.phys.liquid_rho("H2O") == AGNI.consts._lookup_liquid_rho["H2O"]
    @test AGNI.phys.liquid_rho("UNKNOWN") == AGNI.phys.BIGFLOAT

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
    ideal_H2O::phys.Gas_t = phys.load_gas("$RES_DIR/thermodynamics/", "H2O", true, false)
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
            v_obs[i] = phys.get_Cp(ideal_H2O, t_test[i]) * ideal_H2O.mmw # get value and convert units
            test_pass &= isapprox(v_expt[i], v_obs[i]; rtol=1e-3)
        end
        if !test_pass
            @info "Expected values = $(v_expt) J mol-1 K-1"
            @info "Modelled values = $(v_obs) J mol-1 K-1"
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
            v_obs[i] = phys.calc_rho_gas(t_test[i], p_test[i], ideal_H2O)
            test_pass &= isapprox(v_expt[i], v_obs[i]; rtol=1e-3, atol=1e-12)
        end

        if test_pass
            @test_logs (:warn, "Expected values = $(v_expt) kg m-3")
            @test_logs (:warn, "Modelled values = $(v_obs) kg m-3")
        end
        @test test_pass
    end


    # -------------
    # Test AQUA equation of state
    # -------------
    @testset "AQUA_EOS" begin
        aqua_H2O::phys.Gas_t = phys.load_gas("$RES_DIR/thermodynamics/", "H2O", true, true)
        t_test = [200.0,  300.0, 500.0,   1273.0,  3200.0] # Tested values of temperature [K]
        p_test = [1e0,    1e3,   1e5,     1e7,     1e8]    # Tested values of pressure [Pa]
        v_expt = [926.12116198786, 0.007222354920, 0.4333412952269, 17.038999553692, 66.87150907049]
        v_obs  = zero(p_test)
        test_pass = true
        for i in 1:5
            v_obs[i] = phys.calc_rho_gas(t_test[i], p_test[i], aqua_H2O)
            test_pass &= isapprox(v_expt[i], v_obs[i]; rtol=1e-3)
        end

        if !test_pass
            @test_logs (:warn, "Expected values = $(v_expt) kg m-3")
            @test_logs (:warn, "Modelled values = $(v_obs) kg m-3")
        end
        @test test_pass
    end


    # -------------
    # Test VdW equation of state
    # -------------
    @testset "vdw_EOS" begin
        vdw_CO2::phys.Gas_t = phys.load_gas("$RES_DIR/thermodynamics/", "CO2", true, true)
        t_test = [200.0,  300.0, 500.0,   1273.0,  3200.0] # Tested values of temperature [K]
        p_test = [1e0,    1e3,   1e5,     1e7,     1e8]    # Tested values of pressure [Pa]
        v_expt = [2.646533036586e-5, 0.01764355357724, 1.061227115599, 41.2385376189, 147.631823888]
        v_obs  = zero(p_test)
        test_pass = true
        for i in 1:5
            v_obs[i] = phys.calc_rho_gas(t_test[i], p_test[i], vdw_CO2)
            test_pass &= isapprox(v_expt[i], v_obs[i]; rtol=1e-3)
        end
        if test_pass
            @test_logs (:warn, "Expected values = $(v_expt) kg m-3")
            @test_logs (:warn, "Modelled values = $(v_obs) kg m-3")
        end
        @test test_pass
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
            @test_logs (:warn, "Expected values = $(dct_e)")
            @test_logs (:warn, "Modelled values = $(dct_o)")
        end
        atmosphere.deallocate!(atmos)
        @test test_pass
    end
end
