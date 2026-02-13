
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
        energy.calc_fluxes!(atmos, radiative=true, convective=true, conductive=true, sens_heat=true, latent_heat=true, advective=true)
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
        if test_pass
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
            @warn ("Fail \n Observed: $(arr_o_prs)")
        end
        @test test_pass

        # temperature profile
        @debug ("Temperature...")
        arr_e[:] .= [160.7088664780815, 172.19032126178496, 178.16561136706684, 179.33373012757735, 180.00977169552914, 179.38249915236452, 179.63139981968027, 181.7784723456627, 183.23201439666911, 182.84793377402212, 183.58309955137906, 186.39754255672958, 190.50608275317828, 195.51792694295384, 201.2539426805421, 207.6206162995812, 214.39033710208366, 221.79737139979238, 231.75042336773424, 246.49731386931796, 263.55300687852633, 287.98847780008725, 320.73277446340387, 355.3021779910596, 392.6169738903187, 429.6576402027081, 460.8349624615243, 487.31898647282645, 514.7451164371932, 543.3888106490932, 570.4324135285711, 593.6178077887256, 610.3907977993852, 619.9928759088517, 624.69827311216, 627.1707368815037, 628.6896155575437, 629.5684702400786, 629.9508386423341, 630.0682907450105, 630.0966361922129, 630.1029022380585, 630.1042558931739, 630.1043594621661, 630.1043597908382]
        test_pass = all(abs.(arr_e[:] .- arr_o_tmp[:]) .< abs.(arr_e[:].*rtol) .+ 0.5)
        if !test_pass
            @warn ("Fail \n Observed: $(arr_o_tmp)")
        end
        @test test_pass

        # radiative flux
        @debug ("Radiative flux...")
        arr_e[:] .= [-4.9872131683059706e-5, -0.00021685517390324094, -5.446596077263166e-5, -6.956502420507604e-5, -4.464714481855481e-5, -4.1937900391531e-5, -6.189397629441373e-5, -8.169547623992912e-5, -4.90136637267824e-5, -4.301526104200093e-5, -7.095810627788524e-5, -0.00037362132775342616, -0.0016387527005576885, -0.005163568319403566, -0.015416830214064703, -0.04250362273199926, -0.10932251298316942, -0.30647951633102366, -1.4784710919412305, -9.927186737124345, -38.40847329339533, -39.38202174156268, -40.03873661212435, -24.729291535171853, -8.090326572761398, -0.0002255708594987027, -0.0001643575964322963, -0.00015236288642483942, -0.00017095811467982003, -0.00016738894379386693, -0.0001507433066763042, -0.00012077610864480448, -7.536869504320975e-5, -3.700287694208271e-5, -1.7974538663700912e-5, -1.075028051067406e-5, -6.867279150668537e-6, -3.3999808010110044e-6, -1.1308877208085488e-6, -2.810934676129688e-7, -6.233177590498777e-8, -1.3682293434019401e-8, -4.055073132242046e-8, 2.369779394782707e-8, -2.6217860067846884e-5]
        test_pass = all(abs.(arr_e[:] .- arr_o_flN[:]) .< abs.(arr_e[:].*rtol) .+ atol)
        if !test_pass
            @warn ("Fail \n Observed: $(arr_o_flN)")
        end
        @test test_pass

        # convective flux
        arr_e[:] .= [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5717014255377733, 22.0839382511119, 18.39491183960335, 7.661612110746325, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        test_pass = all(abs.(arr_e[:] .- arr_o_flC[:]) .< abs.(arr_e[:].*rtol) .+ atol)
        if !test_pass
            @warn ("Fail \n Observed: $(arr_o_flC)")
        end
        @test test_pass

        # eddy diffusion coefficients
        @debug ("Kzz...")
        arr_e[:] .= [333668.65612205677, 275540.00776454015, 227537.991614982, 187898.4400422248, 155164.52228356962, 128133.20307543786, 105811.02876317696, 87377.61594338887, 72155.50077524007, 59585.24086419178, 49204.85466386036, 40632.8427539614, 33554.16699321721, 27708.672253775378, 22881.52521331677, 18895.318811832065, 15603.552196467479, 12885.246529707978, 10640.498781356933, 8786.810097503372, 7256.053806881907, 5991.971632950747, 13734.470378102425, 12192.377627146934, 8529.677815084575, 8529.677815084575, 8529.677815084575, 8529.677815084575, 8529.677815084575, 8529.677815084575, 8529.677815084575, 8529.677815084575, 8529.677815084575, 8529.677815084575, 8529.677815084575, 8529.677815084575, 8529.677815084575, 8529.677815084575, 8529.677815084575, 8529.677815084575, 8529.677815084575, 8529.677815084575, 8529.677815084575, 8529.677815084575, 8529.677815084575]
        test_pass = all(abs.(arr_e[:] .- arr_o_Kzz[:]) .< abs.(arr_e[:].*rtol) .+ 1e4)
        if !test_pass
            @warn ("Fail \n Observed: $(arr_o_Kzz)")
        end
        @test test_pass
    end
end
