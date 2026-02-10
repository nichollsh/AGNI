
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
            @test_logs (:warn, "Expected value = $(val_e) W m-2")
            @test_logs (:warn, "Modelled value = $(val_o) W m-2")
        end
        @test test_pass

        # Test no-scattering
        val_o = atmos.flux_u_sw[2]
        val_e = 0.0
        @info "Expected value = $(val_e) W m-2"
        @info "Modelled value = $(val_o) W m-2"
        test_pass = isapprox(val_e,val_o; atol=1.0e-10)
        if !test_pass
            @test_logs (:warn, "Expected value = $(val_e) W m-2")
            @test_logs (:warn, "Modelled value = $(val_o) W m-2")
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
            val_e = 1.0438945347722374e7   # known from previous tests
            val_o = atmos.r[1] # height of topmost layer-centre

            test_pass = isapprox(val_e, val_o; rtol=rtol)
            if !test_pass
                @test_logs (:warn,"Expected value = $(val_e) m")
                @test_logs (:warn,"Modelled value = $(val_o) m")
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
                @test_logs (:warn,"Expected range = $(val_e) W m-2")
                @test_logs (:warn,"Modelled value = $(val_o) W m-2")
            end
            @test test_pass
        end


        # -------------
        # Test surface albedo
        # -------------
        @testset "surf_albedo" begin
            val_e = 30.10484135704629  # known from previous tests
            val_o = atmos.flux_u_sw[end] # bottom level
            test_pass = isapprox(val_e, val_o; rtol=1e-3)
            if !test_pass
                @test_logs (:warn,"Expected value = $(val_e) W m-2")
                @test_logs (:warn,"Modelled value = $(val_o) W m-2")
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

        val_e = 37.125717835788286
        val_o = atmos.flux_u_sw[20]
        test_pass = isapprox(val_e, val_o; rtol=1e-3)
        if !test_pass
            @test_logs (:warn, "Expected value = $(val_e) W m-2")
            @test_logs (:warn, "Modelled value = $(val_o) W m-2")
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

        val_e = 6.315779421803655   # from previous tests
        val_o = atmos.heating_rate[atmos.nlev_c-10]
        test_pass = isapprox(val_e, val_o; rtol=1e-3)
        if !test_pass
            @test_logs (:warn, "Expected value = $(val_e) K/day")
            @test_logs (:warn, "Modelled value = $(val_o) K/day")
        end
        @test test_pass


        # -------------
        # Test flux calculation
        # -------------
        energy.calc_fluxes!(atmos, radiative=true, convective=true, conductive=true, sens_heat=true, latent_heat=true, advective=true)
        val_e = 8518.663363443498  # from previous tests
        val_o = atmos.flux_tot[atmos.nlev_c-10]
        test_pass = isapprox(val_e, val_o; rtol=1e-3)
        if !test_pass
            @test_logs (:warn, "Expected value = $(val_e) W m-2")
            @test_logs (:warn, "Modelled value = $(val_o) W m-2")
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
            @test_logs (:warn, "Expected value = $(val_e) W m-2")
            @test_logs (:warn, "Modelled value = $(val_o) W m-2")
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
        test_pass = isapprox(val_e, val_o; rtol=1e-3)
        if test_pass
            @test_logs (:warn, "Expected value = $(val_e) W m-2")
            @test_logs (:warn, "Modelled value = $(val_o) W m-2")
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
        @test_logs (:debug,"Pressure...")
        arr_e[:] .= [1.0, 1.61372830078497, 2.60411902875434, 4.20234057531355, 6.78143591592047, 10.9433950574805, 17.6596663109266, 28.4979033083612, 45.9878730817361, 74.2119322849048, 119.757895384089, 193.256705023749, 311.863814213277, 503.263462986711, 812.130492972705, 1310.5579604405, 2114.88447058187, 3412.8489230686, 5507.4108934593, 8887.46482282671, 14341.9535068263, 23144.0162625079, 37348.1540366365, 60269.7731509967, 97259.0386156535, 156949.663121218, 253274.113177377, 408715.604290549, 659555.937616089, 1064344.08248185, 1717562.16767397, 2771678.67833306, 4472736.32390834, 7217781.18783982, 11647537.7716905, 18795961.3366388, 30331574.7493941, 48946920.5804721, 78987030.9769821, 127463607.282535, 205691630.391968, 331930405.19812, 535645488.759229, 646902370.004265, 648216250.228125 ]
        test_pass = all(abs.(arr_e[:] .- arr_o_prs[:]) .< abs.(arr_e[:].*rtol) .+ 0.5)
        if !test_pass
            @test_logs (:warn, "Fail \n Observed: $(arr_o_prs)")
        end
        @test test_pass

        # temperature profile
        @test_logs (:debug,"Temperature...")
        arr_e[:] .= [156.2474171038719, 170.6236741177102, 177.82871758820932, 178.94959028534362, 179.77037596942523, 179.12404725519193, 179.34793569336753, 181.4834298466808, 182.98942045091283, 182.62303290117092, 183.33435351296674, 186.1488535646183, 190.2748414388916, 195.3117522536504, 201.08240309534239, 207.4856845107921, 214.28519267049717, 221.70968236521145, 231.6842325926396, 246.46951674232272, 263.3838117380063, 287.63649940134104, 320.37606242084445, 354.9305386393268, 392.2046754217648, 429.3311256102809, 460.63198921942706, 487.1091977164687, 514.5095828943975, 543.2112043967588, 570.3990955083275, 593.8081054178406, 610.8348521171442, 620.6303433540376, 625.4327055968067, 627.9447710671932, 629.4879169996773, 630.3857447766541, 630.7797347716225, 630.9020479893142, 630.9317199234893, 630.9382840929698, 630.9397032745547, 630.9398105301121, 630.9398108391669]
        test_pass = all(abs.(arr_e[:] .- arr_o_tmp[:]) .< abs.(arr_e[:].*rtol) .+ 0.5)
        if !test_pass
            @test_logs (:warn, "Fail \n Observed: $(arr_o_tmp)")
        end
        @test test_pass

        # radiative flux
        @test_logs (:debug,"Radiative flux...")
        arr_e[:] .= [-4.9872131683059706e-5, -0.00021685517390324094, -5.446596077263166e-5, -6.956502420507604e-5, -4.464714481855481e-5, -4.1937900391531e-5, -6.189397629441373e-5, -8.169547623992912e-5, -4.90136637267824e-5, -4.301526104200093e-5, -7.095810627788524e-5, -0.00037362132775342616, -0.0016387527005576885, -0.005163568319403566, -0.015416830214064703, -0.04250362273199926, -0.10932251298316942, -0.30647951633102366, -1.4784710919412305, -9.927186737124345, -38.40847329339533, -39.38202174156268, -40.03873661212435, -24.729291535171853, -8.090326572761398, -0.0002255708594987027, -0.0001643575964322963, -0.00015236288642483942, -0.00017095811467982003, -0.00016738894379386693, -0.0001507433066763042, -0.00012077610864480448, -7.536869504320975e-5, -3.700287694208271e-5, -1.7974538663700912e-5, -1.075028051067406e-5, -6.867279150668537e-6, -3.3999808010110044e-6, -1.1308877208085488e-6, -2.810934676129688e-7, -6.233177590498777e-8, -1.3682293434019401e-8, -4.055073132242046e-8, 2.369779394782707e-8, -2.6217860067846884e-5]
        test_pass = all(abs.(arr_e[:] .- arr_o_flN[:]) .< abs.(arr_e[:].*rtol) .+ atol)
        if !test_pass
            @test_logs (:warn, "Fail \n Observed: $(arr_o_flN)")
        end
        @test test_pass

        # convective flux
        arr_e[:] .= [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1758088369679132, 21.882272290150052, 19.522283264947543, 8.769662805655303, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        test_pass = all(abs.(arr_e[:] .- arr_o_flC[:]) .< abs.(arr_e[:].*rtol) .+ atol)
        if !test_pass
            @test_logs (:warn, "Fail \n Observed: $(arr_o_flC)")
        end
        @test test_pass

        # eddy diffusion coefficients
        @test_logs (:debug,"Kzz...")
        arr_e[:] .= [293515.330261667, 242381.028037358, 200155.006214063, 165285.323017843, 136490.40572033, 112711.948729345, 93076.0174630662, 76860.9284503441, 63470.7251477871, 52413.2746248124, 43282.1800995569, 35741.8445533194, 29515.1364634408, 24373.2043307382, 20127.0656526946, 16620.6612101877, 13725.1193905086, 11334.0197421416, 9359.48168174701, 7728.93460077965, 6382.45066279981, 5270.54226322015, 13528.8039526634, 12294.5144847605, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378, 8732.90905977378]
        test_pass = all(abs.(arr_e[:] .- arr_o_Kzz[:]) .< abs.(arr_e[:].*rtol) .+ 1e4)
        if !test_pass
            @test_logs (:warn, "Fail \n Observed: $(arr_o_Kzz)")
        end
        @test test_pass
    end
end
