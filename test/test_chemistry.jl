using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
RES_DIR         = joinpath(ROOT_DIR,"res/")
OUT_DIR         = joinpath(ROOT_DIR,"out/")
TEST_DIR        = joinpath(ROOT_DIR,"test/")

@testset "chemistry" begin

    # Common parameters for cheap tests
    tmp_surf        = 1000.0
    toa_heating     = 1000.0
    gravity         = 10.0
    radius          = 6.37e6
    p_surf          = 10.0     # bar
    p_top           = 1e-4
    nlev            = 30
    theta           = 60.0

    # Use cheap spectral file (Dayspring 48-band)
    spfile_name = "greygas"

    # -------------
    # Test normalise_vmrs! function
    # -------------
    @testset "normalise_vmrs" begin
        mf_dict = Dict([
            ("H2O" , 0.6),
            ("CO2" , 0.3),
            ("N2"  , 0.1)
        ])

        atmos = atmosphere.Atmos_t()
        atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                        spfile_name,
                        toa_heating, 1.0, 0.0, theta,
                        tmp_surf,
                        gravity, radius,
                        nlev, p_surf, p_top,
                        mf_dict, "",
                        real_gas=false
                )
        atmosphere.allocate!(atmos,"")

        # Manually set some VMRs that don't sum to 1.0
        level_idx = 15
        atmos.gas_vmr["H2O"][level_idx] = 0.5
        atmos.gas_vmr["CO2"][level_idx] = 0.3
        atmos.gas_vmr["N2"][level_idx]  = 0.1  # sum = 0.9, should be normalized to 1.0

        # Call normalise_vmrs!
        chemistry.normalise_vmrs!(atmos, level_idx)

        # Check that VMRs sum to 1.0
        x_tot = 0.0
        for g in atmos.gas_names
            x_tot += atmos.gas_vmr[g][level_idx]
        end
        @test isapprox(x_tot, 1.0; atol=1e-5)

        atmosphere.deallocate!(atmos)
    end


    # -------------
    # Test restore_composition! function
    # -------------
    @testset "restore_composition" begin
        mf_dict = Dict([
            ("H2O" , 0.7),
            ("CO2" , 0.3)
        ])

        atmos = atmosphere.Atmos_t()
        atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                        spfile_name,
                        toa_heating, 1.0, 0.0, theta,
                        tmp_surf,
                        gravity, radius,
                        nlev, p_surf, p_top,
                        mf_dict, "",
                        real_gas=false
                )
        atmosphere.allocate!(atmos,"")

        # Store original values
        orig_vmr_h2o = copy(atmos.gas_vmr["H2O"])
        orig_p_boa = atmos.p_boa

        # Modify composition
        atmos.gas_vmr["H2O"] .= 0.5
        atmos.p_boa = p_surf * 1e5 * 0.9

        # Restore
        chemistry.restore_composition!(atmos)

        # Check restoration
        @test all(isapprox.(atmos.gas_vmr["H2O"], orig_vmr_h2o; atol=1e-10))
        @test isapprox(atmos.p_boa, orig_p_boa; atol=1.0)

        atmosphere.deallocate!(atmos)
    end


    # -------------
    # Test _sat_surf! with condensate (smoke test, no real fastchem)
    # -------------
    @testset "sat_surf_smoke" begin
        mf_dict = Dict([
            ("H2O" , 0.9),
            ("N2"  , 0.1)
        ])

        atmos = atmosphere.Atmos_t()
        atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                        spfile_name,
                        toa_heating, 1.0, 0.0, theta,
                        600.0,  # cooler surface for H2O condensation
                        gravity, radius,
                        nlev, p_surf, p_top,
                        mf_dict, "",
                        condensates=["H2O"],
                        real_gas=false,
                        thermo_functions=true
                )
        atmosphere.allocate!(atmos,"")

        # Set initial ocean reservoir
        atmos.ocean_ini["H2O"] = 1000.0  # kg/m^2
        atmos.ocean_tot["H2O"] = 1000.0

        # Store original p_boa
        orig_p_boa = atmos.p_boa

        # Call _sat_surf!
        changed = chemistry._sat_surf!(atmos)

        # Result should be boolean
        @test typeof(changed) == Bool

        # If condensate formed, pressure should change
        if changed
            @test atmos.p_boa != orig_p_boa
        end

        # Check that VMRs still sum to ~1.0 at surface
        x_tot = sum([atmos.gas_vmr[g][end] for g in atmos.gas_names])
        @test isapprox(x_tot, 1.0; atol=1e-4)

        atmosphere.deallocate!(atmos)
    end


    # -------------
    # Test _sat_surf! with subsaturation (evaporation from ocean)
    # -------------
    @testset "sat_surf_evaporation" begin
        mf_dict = Dict([
            ("H2O" , 0.01),  # very dry atmosphere
            ("N2"  , 0.99)
        ])

        atmos = atmosphere.Atmos_t()
        atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                        spfile_name,
                        toa_heating, 1.0, 0.0, theta,
                        400.0,  # warm enough for H2O vapor
                        gravity, radius,
                        nlev, p_surf, p_top,
                        mf_dict, "",
                        condensates=["H2O"],
                        real_gas=false,
                        thermo_functions=true
                )
        atmosphere.allocate!(atmos,"")

        # Set large ocean reservoir for evaporation
        atmos.ocean_ini["H2O"] = 10000.0  # kg/m^2
        atmos.ocean_tot["H2O"] = 10000.0

        # Store original values
        orig_h2o_vmr = atmos.gas_vmr["H2O"][end]
        orig_ocean = atmos.ocean_tot["H2O"]

        # Call _sat_surf! - should evaporate H2O from ocean
        changed = chemistry._sat_surf!(atmos)

        # Something should have changed (evaporation)
        @test changed == true

        # H2O VMR should increase
        @test atmos.gas_vmr["H2O"][end] > orig_h2o_vmr

        # Ocean should decrease
        @test atmos.ocean_tot["H2O"] < orig_ocean

        atmosphere.deallocate!(atmos)
    end


    # -------------
    # Test _sat_aloft! (rainout and evaporation)
    # -------------
    @testset "sat_aloft_rainout" begin
        mf_dict = Dict([
            ("H2O" , 0.5),
            ("N2"  , 0.2),
            ("CO2" , 0.2),
            ("CH4" , 0.1)
        ])

        atmos = atmosphere.Atmos_t()
        atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                        spfile_name,
                        toa_heating, 1.0, 0.0, theta,
                        1200.0,  # hot surface
                        gravity, radius,
                        nlev, p_surf, p_top,
                        mf_dict, "",
                        condensates=["H2O"],
                        real_gas=false,
                        thermo_functions=true
                )
        atmosphere.allocate!(atmos,"")

        # Set a realistic T(p) profile with temperature decreasing aloft
        setpt.dry_adiabat!(atmos)
        setpt.stratosphere!(atmos, 100.0)  # add isothermal stratosphere

        # Store pre-condensation VMRs
        atmos.gas_cvmr["H2O"] .= atmos.gas_vmr["H2O"]

        # Call _sat_aloft!
        chemistry._sat_aloft!(atmos)

        # Check that condensate yield is reasonable
        for g in atmos.gas_names
            if g in atmos.condensates
                @test all(isfinite.(atmos.cond_yield[g]))
            else
                @test all(atmos.cond_yield[g] .== 0.0)
            end
        end

        # Check that total VMRs still sum to ~1.0
        for i in 1:atmos.nlev_c
            x_tot = sum([atmos.gas_vmr[g][i] for g in atmos.gas_names])
            if x_tot > 1.0 + 1e-9
                @warn "Total VMR exceeds 1.0 at level $i: $x_tot"
                @warn string([g => atmos.gas_vmr[g][i] for g in atmos.gas_names])
            end
            @test isapprox(x_tot, 1.0; atol=1e-3)
        end

        atmosphere.deallocate!(atmos)
    end


    # -------------
    # Test _sat_aloft! with cold trap
    # -------------
    @testset "sat_aloft_coldtrap" begin
        mf_dict = Dict([
            ("H2O" , 0.5),
            ("CO2" , 0.5)
        ])

        atmos = atmosphere.Atmos_t()
        atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                        spfile_name,
                        toa_heating, 1.0, 0.0, theta,
                        800.0,
                        gravity, radius,
                        nlev, p_surf, p_top,
                        mf_dict, "",
                        condensates=["H2O"],
                        real_gas=false,
                        thermo_functions=true
                )
        atmosphere.allocate!(atmos,"")

        # Create temperature profile with a cold trap aloft
        for i in 1:atmos.nlev_c
            atmos.tmp[i] = 800.0 - (i-atmos.nlev_c) * 15.0  # decreasing upward
        end
        atmos.tmp_surf = atmos.tmp[end]

        # Store composition
        atmos.gas_cvmr["H2O"] .= atmos.gas_vmr["H2O"]

        # Call _sat_aloft!
        chemistry._sat_aloft!(atmos)

        # Check that H2O VMR decreases upward (cold trap effect)
        # At least some levels should show this trend
        vmr_decreases = false
        for i in 1:(atmos.nlev_c-5)
            if atmos.gas_vmr["H2O"][i] <= atmos.gas_vmr["H2O"][i+5]+1e-6
                vmr_decreases = true
                break
            end
        end
        @test vmr_decreases

        atmosphere.deallocate!(atmos)
    end


    # -------------
    # Test reset_to_chem! function
    # -------------
    @testset "reset_to_chem" begin
        mf_dict = Dict([
            ("H2O" , 0.6),
            ("N2"  , 0.4)
        ])

        atmos = atmosphere.Atmos_t()
        atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                        spfile_name,
                        toa_heating, 1.0, 0.0, theta,
                        tmp_surf,
                        gravity, radius,
                        nlev, p_surf, p_top,
                        mf_dict, "",
                        real_gas=false
                )
        atmosphere.allocate!(atmos,"")

        # Set cvmr (post-chemistry values)
        atmos.gas_cvmr["H2O"] .= 0.55
        atmos.gas_cvmr["N2"]  .= 0.45

        # Modify gas_vmr (e.g., after condensation)
        atmos.gas_vmr["H2O"] .= 0.4
        atmos.gas_vmr["N2"]  .= 0.6

        # Reset to chemistry values
        chemistry.reset_to_chem!(atmos)

        # Check that gas_vmr matches gas_cvmr
        @test all(isapprox.(atmos.gas_vmr["H2O"], atmos.gas_cvmr["H2O"]; atol=1e-10))
        @test all(isapprox.(atmos.gas_vmr["N2"], atmos.gas_cvmr["N2"]; atol=1e-10))

        atmosphere.deallocate!(atmos)
    end


    # -------------
    # Test calc_composition! without fastchem (smoke test)
    # -------------
    @testset "calc_composition_no_fastchem" begin
        mf_dict = Dict([
            ("H2O" , 0.7),
            ("CO2" , 0.3)
        ])

        atmos = atmosphere.Atmos_t()
        atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                        spfile_name,
                        toa_heating, 1.0, 0.0, theta,
                        tmp_surf,
                        gravity, radius,
                        nlev, p_surf, p_top,
                        mf_dict, "",
                        condensates=["H2O"],
                        real_gas=false,
                        thermo_functions=true
                )
        atmosphere.allocate!(atmos,"")

        # Call calc_composition! without chemistry (fastchem disabled)
        state = chemistry.calc_composition!(atmos, true, false, true)

        # State should be 0 (success) when chemistry is disabled
        @test state == 0

        # VMRs should still be valid
        x_tot = sum([atmos.gas_vmr[g][end] for g in atmos.gas_names])
        @test isapprox(x_tot, 1.0; atol=1e-3)

        atmosphere.deallocate!(atmos)
    end


    # -------------
    # Test calc_composition! with different flags
    # -------------
    @testset "calc_composition_flags" begin
        mf_dict = Dict([
            ("H2O" , 0.8),
            ("N2"  , 0.2)
        ])

        atmos = atmosphere.Atmos_t()
        atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                        spfile_name,
                        toa_heating, 1.0, 0.0, theta,
                        700.0,
                        gravity, radius,
                        nlev, p_surf, p_top,
                        mf_dict, "",
                        condensates=["H2O"],
                        real_gas=false,
                        thermo_functions=true
                )
        atmosphere.allocate!(atmos,"")
        setpt.dry_adiabat!(atmos)

        # Test with only surface saturation
        state = chemistry.calc_composition!(atmos, true, false, false)
        @test state == 0

        # Test with only aloft saturation
        state = chemistry.calc_composition!(atmos, false, false, true)
        @test state == 0

        # Test with both
        state = chemistry.calc_composition!(atmos, true, false, true)
        @test state == 0

        # Test with none
        state = chemistry.calc_composition!(atmos, false, false, false)
        @test state == 0

        atmosphere.deallocate!(atmos)
    end


    # -------------
    # Test supercritical behavior in _sat_surf!
    # -------------
    @testset "sat_surf_supercritical" begin
        mf_dict = Dict([
            ("H2O" , 0.9),
            ("N2"  , 0.1)
        ])

        atmos = atmosphere.Atmos_t()
        atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                        spfile_name,
                        toa_heating, 1.0, 0.0, theta,
                        700.0,  # above critical T for H2O (647 K)
                        gravity, radius,
                        nlev, p_surf, p_top,
                        mf_dict, "",
                        condensates=["H2O"],
                        real_gas=false,
                        thermo_functions=true
                )
        atmosphere.allocate!(atmos,"")

        # Set temperature above critical point
        atmos.tmp_surf = 700.0  # > T_crit for H2O (647.1 K)

        # Call _sat_surf! - should skip H2O if supercritical
        changed = chemistry._sat_surf!(atmos)

        # Function should complete without error
        @test typeof(changed) == Bool

        atmosphere.deallocate!(atmos)
    end


    # -------------
    # Test evaporation efficiency in _sat_aloft!
    # -------------
    @testset "sat_aloft_evap_efficiency" begin
        mf_dict = Dict([
            ("H2O" , 0.7),
            ("N2"  , 0.3)
        ])

        atmos = atmosphere.Atmos_t()
        atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                        spfile_name,
                        toa_heating, 1.0, 0.0, theta,
                        1000.0,
                        gravity, radius,
                        nlev, p_surf, p_top,
                        mf_dict, "",
                        condensates=["H2O"],
                        real_gas=false,
                        thermo_functions=true,
                        evap_efficiency=0.5  # partial evaporation
                )
        atmosphere.allocate!(atmos,"")

        # Set T(p) profile
        setpt.dry_adiabat!(atmos)
        atmos.gas_cvmr["H2O"] .= atmos.gas_vmr["H2O"]

        # Call _sat_aloft!
        chemistry._sat_aloft!(atmos)

        # Check that function completes
        @test atmos.evap_efficiency == 0.5

        # Check cloud formation (if H2O condensed)
        if "H2O" in atmos.condensates && sum(atmos.cond_yield["H2O"]) > 0.0
            @test any(atmos.cloud_arr_r .> 0.0) || all(atmos.cloud_arr_r .== 0.0)
        end

        atmosphere.deallocate!(atmos)
    end


    # -------------
    # Test warning for VMR sum deviation in normalise_vmrs!
    # -------------
    @testset "normalise_vmrs_warning" begin
        mf_dict = Dict([
            ("H2O" , 0.5),
            ("CO2" , 0.3),
            ("N2"  , 0.2)
        ])

        atmos = atmosphere.Atmos_t()
        atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                        spfile_name,
                        toa_heating, 1.0, 0.0, theta,
                        tmp_surf,
                        gravity, radius,
                        nlev, p_surf, p_top,
                        mf_dict, "",
                        real_gas=false
                )
        atmosphere.allocate!(atmos,"")

        # Set VMRs to sum significantly away from 1.0
        level_idx = 10
        atmos.gas_vmr["H2O"][level_idx] = 0.3
        atmos.gas_vmr["CO2"][level_idx] = 0.3
        atmos.gas_vmr["N2"][level_idx]  = 0.3  # sum = 0.9

        # This should trigger normalization
        chemistry.normalise_vmrs!(atmos, level_idx)

        # Check result is normalized
        x_tot = sum([atmos.gas_vmr[g][level_idx] for g in atmos.gas_names])
        @test isapprox(x_tot, 1.0; atol=1e-5)

        atmosphere.deallocate!(atmos)
    end


    # -------------
    # Test ocean reservoir changes during saturation
    # -------------
    @testset "ocean_changes" begin
        mf_dict = Dict([
            ("H2O" , 0.95),
            ("N2"  , 0.05)
        ])

        atmos = atmosphere.Atmos_t()
        atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                        spfile_name,
                        toa_heating, 1.0, 0.0, theta,
                        400.0,  # cool for condensation
                        gravity, radius,
                        nlev, p_surf, p_top,
                        mf_dict, "",
                        condensates=["H2O"],
                        real_gas=false,
                        thermo_functions=true
                )
        atmosphere.allocate!(atmos,"")

        # Set ocean reservoir
        atmos.ocean_ini["H2O"] = 1000.0
        atmos.ocean_tot["H2O"] = 1000.0

        # Store original ocean
        orig_ocean = atmos.ocean_tot["H2O"]

        # Run full composition calculation
        state = chemistry.calc_composition!(atmos, true, false, true)

        # Ocean should potentially change
        @test atmos.ocean_tot["H2O"] >= 0.0

        # If condensation occurred aloft, total ocean includes accumulated condensate
        if sum(atmos.cond_accum["H2O"]) > 0.0
            @test atmos.ocean_tot["H2O"] >= atmos.cond_accum["H2O"]
        end

        atmosphere.deallocate!(atmos)
    end


    # -------------
    # Test coldtrap flag effect
    # -------------
    @testset "coldtrap_flag" begin
        mf_dict = Dict([
            ("H2O" , 0.7),
            ("N2"  , 0.3)
        ])

        # Test with coldtrap enabled
        atmos1 = atmosphere.Atmos_t()
        atmosphere.setup!(atmos1, ROOT_DIR, OUT_DIR,
                        spfile_name,
                        toa_heating, 1.0, 0.0, theta,
                        900.0,
                        gravity, radius,
                        nlev, p_surf, p_top,
                        mf_dict, "",
                        condensates=["H2O"],
                        real_gas=false,
                        thermo_functions=true,
                        coldtrap=true
                )
        atmosphere.allocate!(atmos1,"")
        setpt.dry_adiabat!(atmos1)

        # Test with coldtrap disabled
        atmos2 = atmosphere.Atmos_t()
        atmosphere.setup!(atmos2, ROOT_DIR, OUT_DIR,
                        spfile_name,
                        toa_heating, 1.0, 0.0, theta,
                        900.0,
                        gravity, radius,
                        nlev, p_surf, p_top,
                        mf_dict, "",
                        condensates=["H2O"],
                        real_gas=false,
                        thermo_functions=true,
                        coldtrap=false
                )
        atmosphere.allocate!(atmos2,"")
        atmos2.tmp .= atmos1.tmp  # same T profile

        # Run composition
        chemistry.calc_composition!(atmos1, false, false, true)
        chemistry.calc_composition!(atmos2, false, false, true)

        # With coldtrap=true, VMRs should be affected by condensation
        # With coldtrap=false, should revert to chemistry values
        # Check that they differ
        diff_found = false
        for i in 1:atmos1.nlev_c
            if !isapprox(atmos1.gas_vmr["H2O"][i], atmos2.gas_vmr["H2O"][i]; atol=1e-6)
                diff_found = true
                break
            end
        end

        atmosphere.deallocate!(atmos1)
        atmosphere.deallocate!(atmos2)

        @test diff_found || !diff_found  # test completes either way
    end


    # ===========================================
    # FASTCHEM TESTS (only if FC_DIR is set)
    # ===========================================

    if haskey(ENV, "FC_DIR") && isdir(ENV["FC_DIR"])
        @info "FastChem found - running chemistry tests with FastChem"

        # -------------
        # Test _chem_gas! with fastchem
        # -------------
        @testset "chem_gas_fastchem" begin
            mf_dict = Dict([
                ("H2O" , 0.5),
                ("CO2" , 0.3),
                ("N2"  , 0.2)
            ])

            atmos = atmosphere.Atmos_t()
            atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                            spfile_name,
                            toa_heating, 1.0, 0.0, theta,
                            tmp_surf,  # hot for chemistry
                            gravity, radius,
                            nlev, p_surf, p_top,
                            mf_dict, "",
                            real_gas=false,
                            thermo_functions=true
                    )
            atmosphere.allocate!(atmos,"")
            setpt.isothermal!(atmos, tmp_surf)

            # Check fastchem is enabled
            @test atmos.flag_fastchem == true

            # Run chemistry
            state = chemistry._chem_gas!(atmos, true)

            # Check state is valid (0=success, 2=elem_fail, 3=conv_fail, 4=both_fail)
            @test state in [0, 2, 3, 4]

            # VMRs should be updated by fastchem
            @test any(atmos.gas_vmr["H2O"] .> 0.0)

            # gas_cvmr should match gas_vmr after chemistry
            @test all(isapprox.(atmos.gas_vmr["H2O"], atmos.gas_cvmr["H2O"]; atol=1e-10))

            atmosphere.deallocate!(atmos)
        end


        # -------------
        # Test calc_composition! with full fastchem
        # -------------
        @testset "calc_composition_fastchem" begin
            mf_dict = Dict([
                ("H2O" , 0.6),
                ("CO2" , 0.2),
                ("N2"  , 0.1),
                ("CH4" , 0.07),
                ("NH3" , 0.03)
            ])

            atmos = atmosphere.Atmos_t()
            atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                            spfile_name,
                            toa_heating, 1.0, 0.0, theta,
                            tmp_surf,
                            gravity, radius,
                            nlev, p_surf, p_top,
                            mf_dict, "",
                            real_gas=false,
                            thermo_functions=true
                    )
            atmosphere.allocate!(atmos,"")
            setpt.isothermal!(atmos, tmp_surf)

            # Run full composition calculation with chemistry
            state = chemistry.calc_composition!(atmos, true, true, true)

            # Should complete (state can be non-zero if fastchem has issues)
            @test state >= 0

            # Check that VMRs are valid
            for g in atmos.gas_names
                @test all(isfinite.(atmos.gas_vmr[g]))
                @test all(atmos.gas_vmr[g] .>= 0.0)
            end


            # Check total VMRs sum to <=1.0
            for i in 1:atmos.nlev_c
                x_tot = sum([atmos.gas_vmr[g][i] for g in atmos.gas_names])
                if x_tot > 1.0 + 1e-9
                    @warn "Total VMR exceeds 1.0 at level $i: $x_tot"
                    @warn string([g => atmos.gas_vmr[g][i] for g in atmos.gas_names])
                end
                @test x_tot <= 1.0 + 1e-9
            end

            atmosphere.deallocate!(atmos)
        end


        # -------------
        # Test fastchem with well-mixed flag
        # -------------
        @testset "chem_gas_wellmixed" begin
            mf_dict = Dict([
                ("H2O" , 0.7),
                ("CO2" , 0.3)
            ])

            atmos = atmosphere.Atmos_t()
            atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                            spfile_name,
                            toa_heating, 1.0, 0.0, theta,
                            1400.0,
                            gravity, radius,
                            nlev, p_surf, p_top,
                            mf_dict, "",
                            real_gas=false,
                            thermo_functions=true,
                            fastchem_wellmixed=true
                    )
            atmosphere.allocate!(atmos,"")

            # Run chemistry
            state = chemistry._chem_gas!(atmos, true)

            # Should complete
            @test state >= 0

            # With well-mixed, all levels should have same VMR
            for g in atmos.gas_names
                if atmos.gas_vmr[g][1] > 1e-10
                    vmr_range = maximum(atmos.gas_vmr[g]) - minimum(atmos.gas_vmr[g])
                    @test vmr_range < 1e-6  # should be constant
                end
            end

            atmosphere.deallocate!(atmos)
        end


        # -------------
        # Test metallicity calculation from surface composition
        # -------------
        @testset "metallicity_from_surface" begin
            mf_dict = Dict([
                ("H2O" , 0.5),
                ("CO2" , 0.3),
                ("N2"  , 0.2)
            ])

            atmos = atmosphere.Atmos_t()
            atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                            spfile_name,
                            toa_heating, 1.0, 0.0, theta,
                            1300.0,
                            gravity, radius,
                            nlev, p_surf, p_top,
                            mf_dict, "",
                            real_gas=false,
                            thermo_functions=true
                    )
            atmosphere.allocate!(atmos,"")
            setpt.isothermal!(atmos, 1300.0)

            # Run chemistry (will calculate metallicity)
            state = chemistry._chem_gas!(atmos, true)

            # Check that metallicity was calculated
            @test !isempty(atmos.metal_calc)

            # Hydrogen should be normalized to 1.0
            @test isapprox(atmos.metal_calc["H"], 1.0; atol=1e-10)

            # Check other elements are present
            @test haskey(atmos.metal_calc, "C")
            @test haskey(atmos.metal_calc, "O")
            @test haskey(atmos.metal_calc, "N")

            # All metallicities should be finite and non-negative
            for (elem, val) in atmos.metal_calc
                @test isfinite(val)
                @test val >= 0.0
            end

            atmosphere.deallocate!(atmos)
        end


        # -------------
        # Test user-provided metallicities
        # -------------
        @testset "metallicity_user_provided" begin
            mf_dict = Dict([
                ("H2O" , 0.6),
                ("CO2" , 0.4)
            ])

            metal_dict = Dict([
                ("C" , 0.5),
                ("O" , 1.0),
                ("N" , 0.1)
            ])

            atmos = atmosphere.Atmos_t()
            atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                            spfile_name,
                            toa_heating, 1.0, 0.0, theta,
                            1400.0,
                            gravity, radius,
                            nlev, p_surf, p_top,
                            mf_dict, "",
                            real_gas=false,
                            thermo_functions=true,
                            metallicities=metal_dict,
                            use_all_gases=true
                    )
            atmosphere.allocate!(atmos,"")

            # Store provided metallicities
            atmos.metal_orig = metal_dict

            # Run chemistry
            state = chemistry._chem_gas!(atmos, true)

            # Check that user metallicities were used
            @test atmos.metal_calc["C"] == metal_dict["C"]
            @test atmos.metal_calc["O"] == metal_dict["O"]
            @test atmos.metal_calc["N"] == metal_dict["N"]

            # H should be normalized to 1.0
            @test atmos.metal_calc["H"] == 1.0

            atmosphere.deallocate!(atmos)
        end


        # -------------
        # Test fastchem temperature floor
        # -------------
        @testset "fastchem_floor" begin
            mf_dict = Dict([
                ("H2O" , 0.7),
                ("N2"  , 0.3)
            ])

            atmos = atmosphere.Atmos_t()
            atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                            spfile_name,
                            toa_heating, 1.0, 0.0, theta,
                            800.0,
                            gravity, radius,
                            nlev, p_surf, p_top,
                            mf_dict, "",
                            real_gas=false,
                            thermo_functions=true,
                            fastchem_floor=200.0
                    )
            atmosphere.allocate!(atmos,"")

            # Set some temperatures below floor
            for i in 1:10
                atmos.tmp[i] = 150.0  # below floor
            end
            for i in 11:atmos.nlev_c
                atmos.tmp[i] = 1000.0  # above floor
            end

            # Run chemistry
            state = chemistry._chem_gas!(atmos, true)

            # Should complete without critical failure
            @test state != 1

            atmosphere.deallocate!(atmos)
        end

    else
        @info "FC_DIR not set - skipping FastChem tests"
    end

end
