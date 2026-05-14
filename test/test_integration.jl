
using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
RES_DIR         = joinpath(ROOT_DIR,"res/")
OUT_DIR         = joinpath(ROOT_DIR,"out/")
TEST_DIR        = joinpath(ROOT_DIR,"test/")


@testset "integration" begin

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
        energy.radtrans!(atmos, true)
        energy.radtrans!(atmos, false)

        # -------------
        # Test hydrostatic integrator
        # -------------
        @testset "hydrostatic" begin
            val_e = 1.0487713813847492e7   # known from previous tests
            val_o = atmos.r[1] # height of topmost layer-centre

            test_check = isapprox(val_e, val_o; rtol=rtol)
            if !test_check
                @error ("Expected value = $(val_e) m \n Modelled value = $(val_o) m")
            end
            @test test_check
        end


        # -------------
        # Test greenhouse
        # -------------
        @testset "greenhouse" begin
            val_e = [270.0, 280.0]
            val_o = atmos.flux_u_lw[1]
            test_check = ( val_o > val_e[1]) && (val_o < val_e[2])
            if !test_check
                @error ("Expected range = $(val_e) W m-2 \n Modelled value = $(val_o) W m-2")
            end
            @test test_check
        end


        # -------------
        # Test surface albedo
        # -------------
        @testset "surf_albedo" begin
            val_e = 29.699515011666094  # known from previous tests
            val_o = atmos.flux_u_sw[end] # bottom level
            test_check = isapprox(val_e, val_o; rtol=1e-3)
            if !test_check
                @error ("Expected value = $(val_e) W m-2\n Modelled value = $(val_o) W m-2")
            end
            @test test_check
        end


        # -------------
        # Test write NetCDF
        # -------------
        @testset "write_ncdf" begin
            out_path::String = joinpath(OUT_DIR,"agni_atm.nc")
            rm(out_path, force=true)
            save.write_ncdf(atmos, out_path)
            @test isfile(out_path)
            rm(out_path, force=true)
        end

        # -------------
        # Test plots
        # -------------
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
    @testset "total_fluxes_heating" begin
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
        energy.calc_fluxes!(atmos, radiative=true, convective=true, conductive=true, sens_heat=true, latent_heat=true, deep=true)
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
    end

    # -------------
    # Run example TOML config file (writes to output folder)
    # -------------
    @testset "solver" begin
        cfg = AGNI.open_config(joinpath(TEST_DIR, "test.toml"))

        # check return code is fine
        succ = AGNI.run_from_config(cfg)
        @test succ

        # check io directory exists
        @test isdir(cfg["files"]["io_dir"])

        # -------------
        # Compare result from NetCDF
        # -------------

        # Variable to store expected values
        arr_e::Array{Float64,1}     = zeros(Float64, cfg["execution"]["num_levels"]+1)
        arr_o_tmp::Array{Float64,1} = zero(arr_e)
        arr_o_prs::Array{Float64,1} = zero(arr_e)
        arr_o_flN::Array{Float64,1} = zero(arr_e)
        arr_o_flC::Array{Float64,1} = zero(arr_e)
        arr_o_Kzz::Array{Float64,1} = zero(arr_e)

        atol = Float64(cfg["execution"]["converge_atol"])
        rtol = 1e-3

        # Read result
        @debug "ALL DEBUG SUPPRESSED"
        with_logger(MinLevelLogger(current_logger(), Logging.Info-200)) do

            # open netcdf
            ds = Dataset(joinpath(OUT_DIR,"atm.nc"),"r")

            # read profiles from netCDF file
            arr_o_tmp[:] .= ds["tmpl"][:]
            arr_o_prs[:] .= ds["pl"][:]
            arr_o_flN[:] .= ds["fl_N"][:]
            arr_o_flC[:] .= ds["fl_cnvct"][:]
            arr_o_Kzz[:] .= ds["Kzz"][:]

            # close netcdf
            close(ds)
        end
        @debug "ALL DEBUG RESTORED"

        # pressure profile
        @debug ("Pressure...")
        arr_e[:] .= [1.0, 1.6136902021082093, 2.603996068380034, 4.20204294187316, 6.780795524138674, 10.94210329980178, 17.65716488534604, 28.493193972492023, 45.979187940179074, 74.19616507995893, 119.72962462353301, 193.20652215708895, 311.7754717882973, 503.10902408243965, 811.8621027540557, 1310.093920677189, 2114.085723638309, 3411.4794186519816, 5505.070912572516, 8883.478993529165, 14335.183012492122, 23132.544372686574, 37328.7602040377, 60237.05459810247, 97203.94480881537, 156857.05334425246, 253118.69011318486, 408455.1502061103, 659120.073888237, 1.063615605246289e6, 1.7163460809953287e6, 2.7696508543289844e6, 4.469358446891312e6, 7.2121599354580715e6, 1.1638191823886061e7, 1.878043611646084e7, 3.0305805752451997e7, 4.890418180972638e7, 7.891619902847396e7, 1.2734629715987003e8, 2.054974720016427e8, 3.316092571270568e8, 5.3511460915431327e8, 6.462552062782464e8, 6.475677072692152e8]
        test_check = all(abs.(arr_e[:] .- arr_o_prs[:]) .< abs.(arr_e[:].*rtol) .+ 0.5)
        if !test_check
            @error ("Fail 'prs' \n Observed: $(arr_o_prs) \n Expected: $(arr_e)")
        end
        @test test_check

        # temperature profile
        @debug ("Temperature...")
        arr_e[:] .= [160.7127101166621, 172.06733891418713, 178.0617740057226, 179.28281901827108, 179.94763785340575, 179.33949405875046, 179.58738754884365, 181.72978534249938, 183.18320879768305, 182.8016447868322, 183.5353032987346, 186.34839140215888, 190.45553633955188, 195.4643979476455, 201.19521105696884, 207.55331188730077, 214.3042160886999, 221.64134265131094, 231.40105198506797, 245.97587232511935, 263.6284307590409, 288.41210045310174, 321.2183508951017, 356.36918399885417, 393.97230277533515, 431.28194847349585, 462.3693741151927, 488.43279699163713, 515.6021969359344, 543.9069896074404, 570.5354397592575, 593.1878110245217, 609.3083111546545, 618.3878186476225, 622.8567390061546, 625.272523893357, 626.7809360984941, 627.6414677920743, 628.003740829369, 628.110790216126, 628.1357162810501, 628.1410023864355, 628.1420940850873, 628.1421621847325, 628.1421620706775]
        test_check = all(abs.(arr_e[:] .- arr_o_tmp[:]) .< abs.(arr_e[:].*rtol) .+ 0.5)
        if !test_check
            @error ("Fail 'tmp' \n Observed: $(arr_o_tmp) \n Expected: $(arr_e)")
        end
        @test test_check

        # radiative flux
        @debug ("Radiative flux...")
        arr_e[:] .= [0.0417602436712059, 0.040815536211653125, 0.04172261471171623, 0.041653288013094425, 0.041788727242419554, 0.04180431898845427, 0.04168757446302607, 0.04158119980883157, 0.041770728666449486, 0.041801355801226237, 0.04164637467539478, 0.04129035838781192, 0.04002280019369664, 0.03661858175132693, 0.026763352582236166, 0.0009650966914591663, -0.062215841421448204, -0.2433562358924064, -1.2839451770144592, -8.985092341025648, -35.549738097339144, -38.37989351938455, -39.45366404204091, -26.250891294472154, -9.841526399065742, -0.0033711782976411087, -0.0020071803637051744, -0.001607917687522331, -0.0015186961662365661, -0.0013755143447937712, -0.0011838951126890151, -0.0009284323220732915, -0.0005668893964649868, -0.0002817085090942584, -0.00014323282234496304, -8.900588993698832e-5, -5.6428294185950634e-5, -2.6557126526161046e-5, -8.335344284363111e-6, -1.9906749580949484e-6, -4.2501324017039753e-7, -8.790732784095407e-8, 2.364686138706958e-9, -2.141496224811098e-8, 2.5779554562175857e-6]
        test_check = all(abs.(arr_e[:] .- arr_o_flN[:]) .< abs.(arr_e[:].*rtol) .+ atol)
        if !test_check
            @error ("Fail 'flN' \n Observed: $(arr_o_flN) \n Expected: $(arr_e)")
        end
        @test test_check

        # convective flux
        arr_e[:] .= [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 38.54271299376948, 39.61889257514243, 26.246268683682946, 10.214637499173634, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        test_check = all(abs.(arr_e[:] .- arr_o_flC[:]) .< abs.(arr_e[:].*rtol) .+ atol)
        if !test_check
            @error ("Fail 'flC' \n Observed: $(arr_o_flC) \n Expected: $(arr_e)")
        end
        @test test_check

        # eddy diffusion coefficients
        @debug ("Kzz...")
        arr_e[:] .= [989159.0065087004, 816841.9642355904, 674543.5164072266, 557034.2556443322, 459995.7666391871, 379861.9262314778, 313687.8499003042, 259041.6685117144, 213915.1581633123, 176649.93881077305, 145876.52950720975, 120464.02056138197, 99478.51308798684, 82148.79862119754, 67838.0175318585, 56020.25470723378, 46261.21239450335, 38202.249943233546, 31547.203914151887, 26051.504198834664, 21513.186172339287, 17765.468579216737, 17022.526315658855, 13990.776253039552, 9560.569990791864, 9560.569990791864, 9560.569990791864, 9560.569990791864, 9560.569990791864, 9560.569990791864, 9560.569990791864, 9560.569990791864, 9560.569990791864, 9560.569990791864, 9560.569990791864, 9560.569990791864, 9560.569990791864, 9560.569990791864, 9560.569990791864, 9560.569990791864, 9560.569990791864, 9560.569990791864, 9560.569990791864, 9560.569990791864, 9560.569990791864]
        test_check = all(abs.(arr_e[:] .- arr_o_Kzz[:]) .< abs.(arr_e[:].*rtol) .+ 1e4)
        if !test_check
            @error ("Fail 'Kzz' \n Observed: $(arr_o_Kzz) \n Expected: $(arr_e)")
        end
        @test test_check
    end
end
