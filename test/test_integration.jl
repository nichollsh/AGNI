
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
        test_pass = isapprox(val_e, val_o; atol=2)
        if !test_pass
            @warn ("Expected value = $(val_e) W m-2\n Modelled value = $(val_o) W m-2")
        end
        @test test_pass

        # Test no-scattering
        val_o = atmos.flux_u_sw[2]
        val_e = 0.0
        @info "Expected value = $(val_e) W m-2"
        @info "Modelled value = $(val_o) W m-2"
        test_pass = isapprox(val_e,val_o; atol=1.0e-10)
        if !test_pass
            @warn ("Expected value = $(val_e) W m-2\n Modelled value = $(val_o) W m-2")
        end
        @test test_pass
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

            test_pass = isapprox(val_e, val_o; rtol=rtol)
            if !test_pass
                @warn ("Expected value = $(val_e) m \n Modelled value = $(val_o) m")
            end
            @test test_pass
        end


        # -------------
        # Test greenhouse
        # -------------
        @testset "greenhouse" begin
            val_e = [270.0, 280.0]
            val_o = atmos.flux_u_lw[1]
            test_pass = ( val_o > val_e[1]) && (val_o < val_e[2])
            if !test_pass
                @warn ("Expected range = $(val_e) W m-2 \n Modelled value = $(val_o) W m-2")
            end
            @test test_pass
        end


        # -------------
        # Test surface albedo
        # -------------
        @testset "surf_albedo" begin
            val_e = 29.699515011666094  # known from previous tests
            val_o = atmos.flux_u_sw[end] # bottom level
            test_pass = isapprox(val_e, val_o; rtol=1e-3)
            if !test_pass
                @warn ("Expected value = $(val_e) W m-2\n Modelled value = $(val_o) W m-2")
            end
            @test test_pass
        end


        # -------------
        # Test write NetCDF
        # -------------
        @testset "write_ncdf" begin
            out_path::String = joinpath(OUT_DIR,"agni_atm.nc")
            rm(out_path, force=true)
            save.write_ncdf(atmos, out_path)
            if isfile(out_path)
                rm(out_path, force=true)
                @test true
            else
                @test false
            end
        end

        # -------------
        # Test plots
        # -------------
        @testset "plot" begin
            out_path = joinpath(OUT_DIR,"agni_plot_tmp.png")
            rm(out_path, force=true)
            plotting.plot_pt(atmos, out_path)
            @test isfile(out_path)


            out_path = joinpath(OUT_DIR,"agni_plot_rad.png")
            rm(out_path, force=true)
            plotting.plot_radius(atmos, out_path)
            @test isfile(out_path)


            out_path = joinpath(OUT_DIR,"agni_plot_alb.png")
            rm(out_path, force=true)
            plotting.plot_albedo(atmos, out_path)
            @test isfile(out_path)
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
        energy.radtrans!(atmos, true)
        energy.radtrans!(atmos, false)

        val_e = 37.18288051811991
        val_o = atmos.flux_u_sw[20]
        test_pass = isapprox(val_e, val_o; rtol=1e-3)
        if !test_pass
            @warn ("Expected value = $(val_e) W m-2\n Modelled value = $(val_o) W m-2")
        end
        atmosphere.deallocate!(atmos)
        @test test_pass
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
        test_pass = isapprox(val_e, val_o; rtol=1e-3)
        if !test_pass
            @warn ("Expected value = $(val_e) K/day\n Modelled value = $(val_o) K/day")
        end
        @test test_pass


        # -------------
        # Test flux calculation
        # -------------
        energy.calc_fluxes!(atmos, radiative=true, convective=true, conductive=true, sens_heat=true, latent_heat=true, deep=true)
        val_e = 8235.347576033042  # from previous tests
        val_o = atmos.flux_tot[atmos.nlev_c-10]
        test_pass = isapprox(val_e, val_o; rtol=1e-3)
        if !test_pass
            @warn ("Expected value = $(val_e) W m-2\n Modelled value = $(val_o) W m-2")
        end
        @test test_pass

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
        test_pass = isapprox(val_e, val_o; rtol=1e-3)
        if !test_pass
            @warn ("Expected value = $(val_e) W m-2\n Modelled value = $(val_o) W m-2")
        end
        @test test_pass

        # -------------
        # Transparent atmosphere solver (sol_type = 3)
        # -------------
        atmosphere.make_transparent!(atmos)
        atmos.flux_int = 1200.0
        solver.solve_transparent!(atmos; sol_type=3)
        val_e = atmos.flux_int
        val_o = atmos.flux_tot[1]
        test_pass = isapprox(val_e, val_o; rtol=1e-2, atol=0.5)
        if !test_pass
            @warn ("Expected value = $(val_e) W m-2\n Modelled value = $(val_o) W m-2")
        end
        @test test_pass
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
        arr_e[:] .= [1.0, 1.61372830078497, 2.60411902875434, 4.20234057531355, 6.78143591592047, 10.9433950574805, 17.6596663109266, 28.4979033083612, 45.9878730817361, 74.2119322849048, 119.757895384089, 193.256705023749, 311.863814213277, 503.263462986711, 812.130492972705, 1310.5579604405, 2114.88447058187, 3412.8489230686, 5507.4108934593, 8887.46482282671, 14341.9535068263, 23144.0162625079, 37348.1540366365, 60269.7731509967, 97259.0386156535, 156949.663121218, 253274.113177377, 408715.604290549, 659555.937616089, 1064344.08248185, 1717562.16767397, 2771678.67833306, 4472736.32390834, 7217781.18783982, 11647537.7716905, 18795961.3366388, 30331574.7493941, 48946920.5804721, 78987030.9769821, 127463607.282535, 205691630.391968, 331930405.19812, 535645488.759229, 646902370.004265, 648216250.228125 ]
        test_pass = all(abs.(arr_e[:] .- arr_o_prs[:]) .< abs.(arr_e[:].*rtol) .+ 0.5)
        if !test_pass
            @warn ("Fail 'prs' \n Observed: $(arr_o_prs) \n Expected: $(arr_e)")
        end
        @test test_pass

        # temperature profile
        @debug ("Temperature...")
        arr_e[:] .= [160.6951945486835, 172.1746786971098, 178.14914386251837, 179.3174970360867, 179.99367397418297, 179.36630061651562, 179.61486432349, 181.76182287005383, 183.21555860820072, 182.8310793945806, 183.5654069142934, 186.37892037403157, 190.48580094893484, 195.49508157697107, 201.22720345759492, 207.58851577389171, 214.35040002077855, 221.74153252547325, 231.6691046616695, 246.41073491399578, 263.4670022579926, 287.093684861275, 319.6845577217464, 355.27365904095996, 392.97312349421725, 429.96569214612794, 461.0531314431462, 487.49045050571283, 514.8861402680147, 543.5031606896113, 570.52469404162, 593.695649254799, 610.4590899795674, 620.0545016788653, 624.7560661682035, 627.2267737068485, 628.7455649318175, 629.6238504081975, 630.0052214218207, 630.1224213875696, 630.1507066331749, 630.1569595627943, 630.158312274035, 630.158413879461, 630.1584141552416]
        test_pass = all(abs.(arr_e[:] .- arr_o_tmp[:]) .< abs.(arr_e[:].*rtol) .+ 0.5)
        if !test_pass
            @warn ("Fail 'tmp' \n Observed: $(arr_o_tmp) \n Expected: $(arr_e)")
        end
        @test test_pass

        # radiative flux
        @debug ("Radiative flux...")
        arr_e[:] .= [1.1870906746480614e-7, -0.0001695253985758427, -4.829445401810517e-6, -1.9480345088140893e-5, 5.5380432968377136e-6, 7.924799490410805e-6, -1.28889058146342e-5, -3.178004538995083e-5, 1.6810952274681767e-6, 6.602606219985319e-6, -2.1631864171922643e-5, -0.00033096340735028207, -0.001622381046161081, -0.0051919088898557675, -0.015542161698306245, -0.042683713302892556, -0.10943129500662963, -0.30541037706086627, -1.4637689353679093, -9.831207438208537, -37.91071184320069, -40.75303719751662, -37.62829881131762, -23.043609644640355, -7.454170712616104, -0.00023558832512549088, -0.00017526234823606046, -0.00016422760828049832, -0.0001785546421331219, -0.0001712531209037138, -0.00015205559718367567, -0.0001204060679818042, -7.441142472330853e-5, -3.639908325325791e-5, -1.776827483812582e-5, -1.068927611447279e-5, -6.814122859744032e-6, -3.3173632762451e-6, -1.0873069782313394e-6, -2.681373675842043e-7, -5.9402425212115286e-8, -1.3076818613219839e-8, -6.566552094707885e-9, -1.264015736614347e-8, -6.251862942009202e-6]
        test_pass = all(abs.(arr_e[:] .- arr_o_flN[:]) .< abs.(arr_e[:].*rtol) .+ atol)
        if !test_pass
            @warn ("Fail 'flN' \n Observed: $(arr_o_flN) \n Expected: $(arr_e)")
        end
        @test test_pass

        # convective flux
        arr_e[:] .= [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.145716045044864, 21.454737904012404, 17.345291408293214, 7.450548718999271, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        test_pass = all(abs.(arr_e[:] .- arr_o_flC[:]) .< abs.(arr_e[:].*rtol) .+ atol)
        if !test_pass
            @warn ("Fail 'flC' \n Observed: $(arr_o_flC) \n Expected: $(arr_e)")
        end
        @test test_pass

        # eddy diffusion coefficients
        @debug ("Kzz...")
        arr_e[:] .= [427059.004585827, 352660.69825513515, 291223.38308828394, 240489.1139755759, 198593.3043131531, 163996.19868873604, 135426.28376808832, 111833.55761827467, 92350.94002119426, 76262.4055286653, 62976.668084633005, 52005.45007134976, 42945.537123193695, 35463.95918637879, 29285.753198646, 24183.857642759118, 19970.767578287156, 16491.643457278595, 13618.620458920112, 11246.10919976332, 9286.915111152905, 7669.033863158218, 13813.003476924308, 12130.702827887766, 8567.057003230395, 8567.057003230395, 8567.057003230395, 8567.057003230395, 8567.057003230395, 8567.057003230395, 8567.057003230395, 8567.057003230395, 8567.057003230395, 8567.057003230395, 8567.057003230395, 8567.057003230395, 8567.057003230395, 8567.057003230395, 8567.057003230395, 8567.057003230395, 8567.057003230395, 8567.057003230395, 8567.057003230395, 8567.057003230395, 8567.057003230395]
        test_pass = all(abs.(arr_e[:] .- arr_o_Kzz[:]) .< abs.(arr_e[:].*rtol) .+ 1e4)
        if !test_pass
            @warn ("Fail 'Kzz' \n Observed: $(arr_o_Kzz) \n Expected: $(arr_e)")
        end
        @test test_pass
    end
end
