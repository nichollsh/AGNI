using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
OUT_DIR  = joinpath(ROOT_DIR,"out/")

@testset "rfm" begin

    tmp_surf    = 900.0
    toa_heating = 500.0
    gravity     = 10.0
    radius      = 6.37e6
    p_surf      = 10.0
    p_top       = 1e-4
    nlev        = 30
    theta       = 45.0
    mf_dict = Dict([
        ("H2O" , 0.7),
        ("CO2" , 0.3)
    ])

    # this will be used as a small partition file
    parfile = joinpath(RES_DIR, "parfiles", "h2o-co2_4000-5000.par")

    # setup agni atmosphere object
    atmos = atmosphere.Atmos_t()
    atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                    "greygas",
                    toa_heating, 1.0, 0.0, theta,
                    tmp_surf,
                    gravity, radius,
                    nlev, p_surf, p_top,
                    mf_dict, "",
                    rfm_parfile=parfile,
                    real_gas=false
            )
    atmosphere.allocate!(atmos, "")
    setpt.isothermal!(atmos, tmp_surf)

    @testset "write_profile_and_driver" begin

        # write atmosphere profile for RFM unit test
        rfm.write_profile(atmos)
        profile_path = joinpath(atmos.rfm_work, "profile.atm")
        @test isfile(profile_path)
        profile_str = read(profile_path, String)
        @test occursin("*HGT [km]", profile_str)
        @test occursin("*PRE [mb]", profile_str)
        @test occursin("*TEM [K]", profile_str)

        # write RFM config 'driver' file
        rfm.write_driver(atmos, 5000.0, 5050.0, 2.5) # 2.5 clipped to 1.0 in wrapper
        driver_path = joinpath(atmos.rfm_work, "rfm.drv")
        @test isfile(driver_path)
        driver_str = read(driver_path, String)
        @test occursin("*SPC\n  5000.0000 5050.0000 1.0000", driver_str) # spectral range
        @test occursin("*HIT\n  $(atmos.rfm_parfile)", driver_str) # HITRAN partition file
        @test occursin("*OUT\n  radfil=$(joinpath(atmos.rfm_work, "fluxes.asc"))", driver_str)
    end

    @testset "read_smoke_test_fluxes" begin
        # write fake RFM output for testing flux reading function
        fluxpath = joinpath(atmos.rfm_work, "fluxes.asc")
        open(fluxpath, "w") do f
            write(f, "! synthetic RFM output\n")
            write(f, "4 4000.0 0.5\n")
            write(f, "1.0 2.0 3.0 4.0\n")
        end

        # check read ok
        @test rfm.read_fluxes(atmos)
        @test atmos.rfm_npts == 4
        @test length(atmos.rfm_wn) == 4
        @test isapprox(atmos.rfm_wn[end], 4001.5; atol=1e-10)
        @test all(atmos.rfm_fl .> 0.0)
    end

    # test that read_fluxes fails on irregular grid
    @testset "read_smoke_test_fluxes_irreg" begin
        fluxpath = joinpath(atmos.rfm_work, "fluxes.asc")
        open(fluxpath, "w") do f
            write(f, "! synthetic RFM output\n")
            write(f, "2 4000.0 0.0\n")
            write(f, "1.0 2.0\n")
        end
        @test !rfm.read_fluxes(atmos)
    end

    # run rfm and check that it completes and produces expected output files
    @testset "run_rfm" begin
        # remove old file
        rm(joinpath(atmos.rfm_work, "fluxes.asc"), force=true)

        # run rfm
        success = rfm.run_rfm(atmos, 4200.0, 4220.0; nures=1.0, keeplog=true)

        # returns true on success
        @test success isa Bool
        @test success

        # input files exist
        @test isfile(joinpath(atmos.rfm_work, "profile.atm"))
        @test isfile(joinpath(atmos.rfm_work, "rfm.drv"))

        @testset "run_rfm_output" begin
            # output file exists
            fluxpath = joinpath(atmos.rfm_work, "fluxes.asc")
            @test isfile(fluxpath)

            # output file has expected format
            flux_str = read(fluxpath, String)
            @test occursin("RFM flux observing downwards from TOA", flux_str) # header line
        end

        @testset "run_rfm_fluxes_real" begin
            # check that fluxes read ok
            @test rfm.read_fluxes(atmos)
            @test length(atmos.rfm_wn) == atmos.rfm_npts
            @test length(atmos.rfm_fl) == atmos.rfm_npts
            @test all(atmos.rfm_fl .> 0.0)

            # check that fluxes are as expected
            fl_exp = Float64[3368.0288502914664, 3365.044808509454, 3362.062745930814, 3359.082913882957, 3356.105249534032, 3353.1293444769926, 3350.155858446296, 3347.184383034898, 3344.2148239950193, 3341.2474954859235, 3338.282114764274, 3335.318964573408, 3332.3577621699874, 3329.3987902973504, 3326.441766212159, 3323.486689914413, 3320.533844147451, 3317.583197495346, 3314.6345300466137, 3311.687810385326, 3308.743289838896]
            fl_obs = atmos.rfm_fl
            test_check = all(abs.(fl_exp[:] .- fl_obs[:]) .< abs.(fl_exp[:].*rtol) .+ 0.5)
            if !test_check
                @error ("Fail 'rfm_fl' \n Observed: $(fl_obs) \n Expected: $(fl_exp)")
            end
            @test test_check

        end

    end

    atmosphere.deallocate!(atmos)
end
