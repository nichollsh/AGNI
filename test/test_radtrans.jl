
using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
RES_DIR         = joinpath(ROOT_DIR,"res/")
OUT_DIR         = joinpath(ROOT_DIR,"out/")
TEST_DIR        = joinpath(ROOT_DIR,"test/")


@testset "radtrans" begin

    # -------------
    # Test instellation
    # -------------
    @testset "instellation" begin
        tmp_surf        = 200.0    # Surface temperature [kelvin]
        toa_heating     = 1000.00     # Instellation flux [W m-2]
        p_surf          = 50.0    # bar
        theta           = 65.0
        mf_dict         = Dict([
                                ("H2O" , 0.8),
                                ("CO2" , 0.2),
                                ])
        spfile_name   = "$RES_DIR/spectral_files/Dayspring/48/Dayspring.sf"

        # Setup atmosphere
        atmos = atmosphere.Atmos_t()
        atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                                spfile_name,
                                toa_heating, 1.0, 0.0, theta,
                                tmp_surf,
                                gravity, radius,
                                nlev_centre, p_surf, p_top,
                                mf_dict,"",
                                flag_gcontinuum=false,
                                flag_rayleigh=false,
                                overlap_method="ro",
                                real_gas=false
                        )
        atmosphere.allocate!(atmos,"$RES_DIR/stellar_spectra/sun.txt")
        energy.radtrans!(atmos, false)

        # Test absorbed flux
        val_e = toa_heating * cosd(theta)
        val_o = atmos.flux_d_sw[1]
        test_check = isapprox(val_e, val_o; atol=2)
        if !test_check
            @error ("Expected value = $(val_e) W m-2\n Modelled value = $(val_o) W m-2")
        end
        @test test_check

        # Test no-scattering
        val_o = atmos.flux_u_sw[2]
        val_e = 0.0
        @info "Expected value = $(val_e) W m-2"
        @info "Modelled value = $(val_o) W m-2"
        test_check = isapprox(val_e,val_o; atol=1.0e-10)
        if !test_check
            @error ("Expected value = $(val_e) W m-2\n Modelled value = $(val_o) W m-2")
        end
        @test test_check
        atmosphere.deallocate!(atmos)
    end


    # -------------
    # Test radiative transfer effects
    # -------------
    @testset "radtrans" begin
        tmp_surf        = 1300.0    # Surface temperature [kelvin]
        toa_heating     = 1000.00    # Instellation flux [W m-2]
        p_surf          = 300.0    # bar
        theta           = 45.0
        mf_dict         = Dict([
                                ("H2O" , 1.0),
                                ("N2"  , 1.0e-9)
                                ])
        spfile_name   = "$RES_DIR/spectral_files/Oak/318/Oak.sf"

        # Setup atmosphere
        atmos = atmosphere.Atmos_t()
        atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                                spfile_name,
                                toa_heating, 1.0, 0.0, theta,
                                tmp_surf,
                                gravity, radius,
                                300, p_surf, p_top,
                                mf_dict,"",
                                flag_gcontinuum=true,
                                flag_rayleigh=false,
                                overlap_method="ee",
                                surface_material="greybody",
                                albedo_s=0.5,
                                real_gas=false
                        )
        atmosphere.allocate!(atmos,"$RES_DIR/stellar_spectra/sun.txt")
        setpt.dry_adiabat!(atmos)
        setpt.saturation!(atmos, "H2O")
        atmosphere.calc_layer_props!(atmos)
        energy.radtrans!(atmos, true; calc_cf=true)
        energy.radtrans!(atmos, false)

        # Test hydrostatic integrator
        @testset "hydrostatic" begin

            # check shape
            @test size(atmos.r) == (atmos.nlev_c,)

            # check that radius increases with height
            @test all(diff(atmos.r) .< 0.0)

            # check that gravity decreases with height
            @test all(diff(atmos.g) .> 0.0)

            # check that acceleration and gravity are the same, in this case
            @test all(isapprox.(atmos.a, atmos.g; atol=1e-5))

            # check surface value equals radius of planet
            @test isapprox(atmos.rl[end], atmos.rp; atol=1e-5)

            # check pressure decreases with height
            @test all(diff(atmos.pl) .> 0.0)
            @test all(diff(atmos.p) .> 0.0)

            # check all layers are bound
            @test all(atmos.layer_isbound)

            # check known value
            val_e = 1.0487713813847492e7 / AGNI.consts.R_earth
            val_o = atmos.r[1] / AGNI.consts.R_earth # height of topmost layer-centre

            test_check = isapprox(val_e, val_o; rtol=rtol)
            if !test_check
                @error ("Expected value = $(val_e) Rearth \n Modelled value = $(val_o) Rearth")
            end
            @test test_check
        end

        # Test contribution function
        @testset "contrib_function" begin
            # check all >0
            @test all(atmos.contfunc_band .>= 0.0)
        end


        # Test optical depth
        @testset "optical_depth" begin
            # check all >0
            @test all(atmos.tau_band .>= 0.0)

            # check shape
            @test size(atmos.tau_band) == (atmos.nlev_c, atmos.nbands)

            # check that tau increases with depth, over each band
            for b in 1:atmos.nbands
                 @test all(diff(atmos.tau_band[:,b]) .>= 0.0)
            end
        end

        # Test OLR is positive and summed correctly
        @testset "olr_summed" begin

            # bolometric value is positive
            @test all(atmos.flux_u_lw .>= 0.0)

            # sum over bands equals bolometric value
            val_e = sum(atmos.band_u_lw[1,:])
            val_o = atmos.flux_u_lw[1]
            @test isapprox(val_e, val_o; rtol=rtol)
        end



        # Test greenhouse
        @testset "greenhouse_limit" begin
            val_e = [270.0, 280.0]
            val_o = atmos.flux_u_lw[1]
            test_check = ( val_o > val_e[1]) && (val_o < val_e[2])
            if !test_check
                @error ("Expected range = $(val_e) W m-2 \n Modelled value = $(val_o) W m-2")
            end
            @test test_check
        end


        # Test surface albedo
        @testset "surf_albedo" begin

            # check that all upward fluxes are positive
            @test all(atmos.flux_u_sw .>= 0.0)

            # check known value
            val_e = 29.699515011666094  # known from previous tests
            val_o = atmos.flux_u_sw[end] # bottom level
            test_check = isapprox(val_e, val_o; rtol=1e-3)
            if !test_check
                @error ("Expected value = $(val_e) W m-2\n Modelled value = $(val_o) W m-2")
            end
            @test test_check
        end


        # Test write NetCDF
        @testset "write_ncdf" begin
            out_path::String = joinpath(OUT_DIR,"agni_atm.nc")
            rm(out_path, force=true)
            save.write_ncdf(atmos, out_path)
            @test isfile(out_path)
            rm(out_path, force=true)
        end

        # Test plots
        @testset "plot" begin
            function _assert_plot_output(path::String, plt)
                @test plt !== nothing
                @test isfile(path)
                @test filesize(path) > 0
            end

            out_path = joinpath(OUT_DIR, "agni_plot_tmp.png")
            rm(out_path, force=true)
            plt_pt = plotting.plot_pt(atmos, out_path; incl_magma=true, title="TP profile")
            _assert_plot_output(out_path, plt_pt)

            # no-save branch
            @test plotting.plot_pt(atmos, ""; incl_magma=false, title="") !== nothing

            out_path = joinpath(OUT_DIR, "agni_plot_rad.png")
            rm(out_path, force=true)
            plt_ra = plotting.plot_radius(atmos, out_path; title="Radius profile")
            _assert_plot_output(out_path, plt_ra)

            out_path = joinpath(OUT_DIR, "agni_plot_alb.png")
            rm(out_path, force=true)
            plt_alb = plotting.plot_albedo(atmos, out_path)
            _assert_plot_output(out_path, plt_alb)

            out_path = joinpath(OUT_DIR, "agni_plot_flux.png")
            rm(out_path, force=true)
            plt_fl = plotting.plot_fluxes(atmos, out_path;
                                            incl_eff=true, incl_mlt=true,
                                            incl_cdct=true, incl_latent=true,
                                            incl_deep=true, title="Flux budget")
            _assert_plot_output(out_path, plt_fl)

            out_path = joinpath(OUT_DIR, "agni_plot_flux_radonly.png")
            rm(out_path, force=true)
            plt_fl_rad = plotting.plot_fluxes(atmos, out_path;
                                                incl_eff=false, incl_mlt=false,
                                                incl_cdct=false, incl_latent=false,
                                                incl_deep=false, title="")
            _assert_plot_output(out_path, plt_fl_rad)

            out_path = joinpath(OUT_DIR, "agni_plot_vmr.png")
            rm(out_path, force=true)
            plt_vmr = plotting.plot_vmr(atmos, out_path)
            _assert_plot_output(out_path, plt_vmr)

            out_path = joinpath(OUT_DIR, "agni_plot_cloud.png")
            rm(out_path, force=true)
            plt_cld = plotting.plot_cloud(atmos, out_path; title="Cloud profile")
            _assert_plot_output(out_path, plt_cld)

            out_path = joinpath(OUT_DIR, "agni_plot_emission.png")
            rm(out_path, force=true)
            plt_emis = plotting.plot_emission(atmos, out_path)
            _assert_plot_output(out_path, plt_emis)

            out_path = joinpath(OUT_DIR, "agni_plot_contfunc1.png")
            rm(out_path, force=true)
            plt_cf1 = plotting.plot_contfunc1(atmos, out_path; cf_min=1e-7)
            _assert_plot_output(out_path, plt_cf1)

            out_path = joinpath(OUT_DIR, "agni_plot_contfunc2.png")
            rm(out_path, force=true)
            plt_cf2 = plotting.plot_contfunc2(atmos, out_path)
            _assert_plot_output(out_path, plt_cf2)

            out_path = joinpath(OUT_DIR, "agni_plot_tau.png")
            rm(out_path, force=true)
            plt_tau = plotting.plot_tau(atmos, out_path)
            _assert_plot_output(out_path, plt_tau)

            out_path = joinpath(OUT_DIR, "agni_plot_jacobian.png")
            rm(out_path, force=true)
            jac = reshape(collect(range(-5.0, 5.0, length=100)), 10, 10)
            pert = fill(true, 10)
            pert[2:2:end] .= false
            plt_jac = plotting.jacobian(jac, out_path; perturb=pert)
            _assert_plot_output(out_path, plt_jac)

            out_path = joinpath(OUT_DIR, "agni_plot_combined.png")
            rm(out_path, force=true)
            plt_combo = plotting.combined(plt_pt, plt_fl, plt_vmr, plt_ra, plt_cld, plt_jac, "integration diagnostics", out_path)
            _assert_plot_output(out_path, plt_combo)

            # Guard branches when radiance outputs are unavailable
            orig_is_out_lw = atmos.is_out_lw
            orig_is_out_sw = atmos.is_out_sw
            atmos.is_out_lw = false
            atmos.is_out_sw = false

            out_path = joinpath(OUT_DIR, "agni_plot_alb_skip.png")
            rm(out_path, force=true)
            @test plotting.plot_albedo(atmos, out_path) === nothing
            @test !isfile(out_path)

            out_path = joinpath(OUT_DIR, "agni_plot_emission_skip.png")
            rm(out_path, force=true)
            @test plotting.plot_emission(atmos, out_path) === nothing
            @test !isfile(out_path)

            out_path = joinpath(OUT_DIR, "agni_plot_contfunc1_skip.png")
            rm(out_path, force=true)
            @test plotting.plot_contfunc1(atmos, out_path) === nothing
            @test !isfile(out_path)

            out_path = joinpath(OUT_DIR, "agni_plot_contfunc2_skip.png")
            rm(out_path, force=true)
            @test plotting.plot_contfunc2(atmos, out_path) === nothing
            @test !isfile(out_path)

            out_path = joinpath(OUT_DIR, "agni_plot_tau_skip.png")
            rm(out_path, force=true)
            @test plotting.plot_tau(atmos, out_path) === nothing
            @test !isfile(out_path)

            atmos.is_out_lw = orig_is_out_lw
            atmos.is_out_sw = orig_is_out_sw
        end

        atmosphere.deallocate!(atmos)
    end



    # -------------
    # Test Rayleigh scattering
    # -------------
    @testset "rayleigh" begin
        tmp_surf        = 400.0    # Surface temperature [kelvin]
        toa_heating     = 1000.00    # Instellation flux [W m-2]
        p_surf          = 10.0    # bar
        theta           = 75.0
        mf_dict         = Dict([
                                ("H2O" , 0.6),
                                ("CO2" , 0.4),
                                ])
        spfile_name   = "$RES_DIR/spectral_files/Dayspring/48/Dayspring.sf"

        # Setup atmosphere
        atmos = atmosphere.Atmos_t()
        atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                                spfile_name,
                                toa_heating, 1.0, 0.0, theta,
                                tmp_surf,
                                gravity, radius,
                                nlev_centre, p_surf, p_top,
                                mf_dict, "",
                                flag_gcontinuum=true,
                                flag_rayleigh=true,
                                overlap_method="ro",
                                real_gas=false
                        )
        atmosphere.allocate!(atmos,"$RES_DIR/stellar_spectra/sun.txt")
        atmosphere.calc_layer_props!(atmos)
        energy.radtrans!(atmos, true)
        energy.radtrans!(atmos, false)

        # check less than instellation
        @test atmos.flux_u_sw[1] < toa_heating

        # check that bolometric equals sum of total over bands
        val_e = sum(atmos.band_u_sw[1,:])
        val_o = atmos.flux_u_sw[1]
        @test isapprox(val_e, val_o; rtol=rtol)

        # check known value
        val_e = 37.18288051811991
        val_o = atmos.flux_u_sw[20]
        test_check = isapprox(val_e, val_o; rtol=1e-3)
        if !test_check
            @error ("Expected value = $(val_e) W m-2\n Modelled value = $(val_o) W m-2")
        end
        atmosphere.deallocate!(atmos)
        @test test_check
    end



    # -------------
    # Test heating rate calculation
    # -------------
    @testset "total_fluxes" begin
        tmp_surf        = 2500.0    # Surface temperature [kelvin]
        toa_heating     = 1000.0    # Instellation flux [W m-2]
        p_surf          = 5.0     # bar
        theta           = 45.0
        mf_dict         = Dict([
                                ("H2O" , 1.0)
                                ])
        spfile_name   = "$RES_DIR/spectral_files/Oak/318/Oak.sf"

        # Setup atmosphere
        atmos = atmosphere.Atmos_t()
        atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                                spfile_name,
                                toa_heating, 1.0, 0.0, theta,
                                tmp_surf,
                                gravity, radius,
                                nlev_centre, p_surf, p_top,
                                mf_dict, "",
                                flag_gcontinuum=true,
                                flag_rayleigh=false,
                                overlap_method="ro",
                                thermo_functions=false,
                                real_gas=false
                        )
        atmosphere.allocate!(atmos,"$RES_DIR/stellar_spectra/sun.txt")
        setpt.isothermal!(atmos, 300.0)
        atmosphere.calc_layer_props!(atmos)
        atmos.flux_tot[:] .= 0.0
        energy.radtrans!(atmos, true)
        energy.radtrans!(atmos, false)
        atmos.flux_tot += atmos.flux_n
        energy.calc_hrates!(atmos)

        val_e = 6.144974916820797   # from previous tests
        val_o = atmos.heating_rate[atmos.nlev_c-10]
        test_check = isapprox(val_e, val_o; rtol=1e-3)
        if !test_check
            @error ("Expected value = $(val_e) K/day\n Modelled value = $(val_o) K/day")
        end
        @test test_check


        # -------------
        # Test flux calculation
        # -------------

        # evaluate
        energy.calc_fluxes!(atmos, radiative=true, convective=true, conductive=true, sens_heat=true, latent_heat=true, deep=true)

        # check that sum of flux components equals total flux
        val_e = atmos.flux_n + atmos.flux_cdct + atmos.flux_cdry + atmos.flux_l + atmos.flux_deep
        val_e[end] += atmos.flux_sens
        val_o = atmos.flux_tot
        for i in 1:atmos.nlev_c
            test_check = isapprox(val_e[i], val_o[i]; rtol=rtol)
            if !test_check
                @error ("At level $i: Expected value = $(val_e[i]) W m-2\n Modelled value = $(val_o[i]) W m-2")
            end
            @test test_check
        end

        # check against known value
        val_e = 8235.347576033042  # from previous tests
        val_o = atmos.flux_tot[atmos.nlev_c-10]
        test_check = isapprox(val_e, val_o; rtol=1e-3)
        if !test_check
            @error ("Expected value = $(val_e) W m-2\n Modelled value = $(val_o) W m-2")
        end
        @test test_check

        # -------------
        # Transparent atmosphere solver (sol_type = 4)
        # -------------
        @info " "
        @info "Testing transparent solver (sol_type=4)"
        atmosphere.make_transparent!(atmos)
        atmos.target_olr = 5000.0
        solver.solve_transparent!(atmos; sol_type=4)
        val_e = atmos.target_olr
        val_o = atmos.flux_u_lw[1]
        test_check = isapprox(val_e, val_o; rtol=1e-3)
        if !test_check
            @error ("Expected value = $(val_e) W m-2\n Modelled value = $(val_o) W m-2")
        end
        @test test_check

        # -------------
        # Transparent atmosphere solver (sol_type = 3)
        # -------------
        atmosphere.make_transparent!(atmos)
        atmos.flux_int = 1200.0
        solver.solve_transparent!(atmos; sol_type=3)
        val_e = atmos.flux_int
        val_o = atmos.flux_tot[1]
        test_check = isapprox(val_e, val_o; rtol=1e-2, atol=0.5)
        if !test_check
            @error ("Expected value = $(val_e) W m-2\n Modelled value = $(val_o) W m-2")
        end
        @test test_check


        atmosphere.deallocate!(atmos)
    end
end
