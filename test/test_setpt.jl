using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
RES_DIR = joinpath(ROOT_DIR,"res/")
OUT_DIR = joinpath(ROOT_DIR,"out/")

@testset "setpt" begin

    p_surf = 1.0       # bar
    p_top = 1e-8
    theta = 65.0
    gravity = 10.0
    nlev_centre = 50
    radius = 1.0e7
    tmp_surf = 300.0
    toa_heating = 1000.0
    mf_dict = Dict("H2O" => 1.0)
    spfile = "greygas"

    atmos = AGNI.atmosphere.Atmos_t()
    AGNI.atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                          spfile,
                          toa_heating, 1.0, 0.0, theta,
                          tmp_surf,
                          gravity, radius,
                          nlev_centre, p_surf, p_top,
                          mf_dict, ""
                  )
    AGNI.atmosphere.allocate!(atmos, "")  # Empty stellar spectrum string

    # Test _parse_tmp_str with different inputs
    @testset "parse_tmp_str" begin
        # Test with numeric string
        result = AGNI.setpt._parse_tmp_str(atmos, "500.0")
        @test isapprox(result, 500.0; atol=1e-10)

        # Test with numeric value
        result = AGNI.setpt._parse_tmp_str(atmos, 600.0)
        @test isapprox(result, 600.0; atol=1e-10)

        # Test with "tsurf" keyword
        result = AGNI.setpt._parse_tmp_str(atmos, "tsurf")
        @test isapprox(result, atmos.tmp_surf; atol=1e-10)

        # Test with "teq" keyword
        result = AGNI.setpt._parse_tmp_str(atmos, "teq")
        expected_teq = AGNI.phys.calc_Teq(atmos.instellation, atmos.albedo_b)
        @test isapprox(result, expected_teq; rtol=1e-3)

        # Test with invalid input (should error)
        @test_throws ErrorException AGNI.setpt._parse_tmp_str(atmos, Dict())
    end

    @testset "isothermal!" begin
        # Test with numeric value
        AGNI.setpt.isothermal!(atmos, 400.0)
        @test all(atmos.tmp .≈ 400.0)
        @test all(atmos.tmpl .≈ 400.0)

        # Test with string
        AGNI.setpt.isothermal!(atmos, "500.0")
        @test all(atmos.tmp .≈ 500.0)
        @test all(atmos.tmpl .≈ 500.0)

        # Test with keyword
        AGNI.setpt.isothermal!(atmos, "tsurf")
        @test all(atmos.tmp .≈ atmos.tmp_surf)
        @test all(atmos.tmpl .≈ atmos.tmp_surf)
    end

    @testset "add!" begin
        # Set initial isothermal profile
        AGNI.setpt.isothermal!(atmos, 300.0)

        # Add 50 K
        AGNI.setpt.add!(atmos, 50.0)
        @test all(atmos.tmp .≈ 350.0)
        @test all(atmos.tmpl .≈ 350.0)

        # Test with negative value
        AGNI.setpt.add!(atmos, -25.0)
        @test all(atmos.tmp .≈ 325.0)
        @test all(atmos.tmpl .≈ 325.0)
    end

    @testset "dry_adiabat!" begin
        AGNI.setpt.dry_adiabat!(atmos)

        # Check surface temperature is preserved
        @test atmos.tmpl[end] ≈ atmos.tmp_surf

        # Check temperatures are positive and decreasing with altitude
        @test all(atmos.tmp .> 0.0)
        @test all(atmos.tmpl .> 0.0)

        # Temperature should generally decrease with altitude (lower pressure)
        # Check that temperatures are monotonic (allowing for small violations)
        decreasing_count = sum(atmos.tmp[1:end-1] .< atmos.tmp[2:end])
        @test decreasing_count > atmos.nlev_c * 0.8  # At least 80% decreasing
    end

    @testset "stratosphere!" begin
        # Set initial profile
        AGNI.setpt.dry_adiabat!(atmos)
        initial_tmp = copy(atmos.tmp)

        # Apply stratosphere at 200K
        AGNI.setpt.stratosphere!(atmos, 200.0)

        # Check that some upper levels are capped at 200K
        @test any(atmos.tmp .≈ 200.0)

        # Check that lower levels are not affected
        @test any(atmos.tmp .> 200.0)

        # Test with string input
        AGNI.setpt.dry_adiabat!(atmos)
        AGNI.setpt.stratosphere!(atmos, "250.0")
        @test any(atmos.tmp .≈ 250.0)
    end

    @testset "loglinear!" begin
        # Set log-linear profile from 300K at surface to 150K at top
        AGNI.setpt.loglinear!(atmos, 150.0)

        # Check surface temperature
        @test atmos.tmpl[end] ≈ atmos.tmp_surf

        # Temperatures should be reasonable and within bounds
        @test all(atmos.tmp .> 0.0)
        @test all(atmos.tmpl .> 0.0)

        # Top should be cooler than or equal to bottom
        @test atmos.tmpl[1] <= atmos.tmpl[end]

        # Test with keyword
        AGNI.setpt.loglinear!(atmos, "tsurf")
        @test atmos.tmpl[end] ≈ atmos.tmp_surf
    end

    @testset "fromarrays!" begin
        # Create test pressure and temperature arrays
        test_pl = [1e2, 1e3, 1e4, 1e5]  # Pa
        test_tmpl = [150.0, 200.0, 250.0, 300.0]  # K

        AGNI.setpt.fromarrays!(atmos, test_pl, test_tmpl)

        # Check that temperatures are set and reasonable
        @test all(atmos.tmp .> 0.0)
        @test all(atmos.tmpl .> 0.0)
        @test all(atmos.tmp .>= atmos.tmp_floor)
        @test all(atmos.tmp .<= atmos.tmp_ceiling)

        # Test with reversed arrays (should auto-flip)
        test_pl_rev = reverse(test_pl)
        test_tmpl_rev = reverse(test_tmpl)
        AGNI.setpt.fromarrays!(atmos, test_pl_rev, test_tmpl_rev)
        @test all(atmos.tmp .> 0.0)

        # Test extrapolation option
        AGNI.setpt.fromarrays!(atmos, test_pl, test_tmpl; extrap=true)
        @test all(atmos.tmp .> 0.0)
    end

    @testset "fromcsv!" begin
        # Create a temporary CSV file
        tmpfile = tempname() * ".csv"
        open(tmpfile, "w") do io
            println(io, "# Pressure [Pa], Temperature [K]")
            println(io, "1e2, 150.0")
            println(io, "1e3, 200.0")
            println(io, "1e4, 250.0")
            println(io, "1e5, 300.0")
        end

        try
            AGNI.setpt.fromcsv!(atmos, tmpfile)

            # Check that temperatures are set
            @test all(atmos.tmp .> 0.0)
            @test all(atmos.tmpl .> 0.0)
        finally
            rm(tmpfile, force=true)
        end

        # Test with non-existent file (should log error but not throw)
        AGNI.setpt.fromcsv!(atmos, "/nonexistent/file.csv")
    end

    @testset "request!" begin
        # Test single request
        result = AGNI.setpt.request!(atmos, ["iso", 400.0])
        @test result == true
        @test all(atmos.tmp .≈ 400.0)

        # Test dry adiabat request
        result = AGNI.setpt.request!(atmos, ["dry"])
        @test result == true
        @test atmos.tmpl[end] ≈ atmos.tmp_surf

        # Test stratosphere request
        result = AGNI.setpt.request!(atmos, ["dry", "str", 200.0])
        @test result == true

        # Test loglinear request
        result = AGNI.setpt.request!(atmos, ["loglin", 150.0])
        @test result == true

        # Test add request
        result = AGNI.setpt.request!(atmos, ["iso", 300.0, "add", 50.0])
        @test result == true
        @test all(atmos.tmp .≈ 350.0)

        # Test analytic request
        result = AGNI.setpt.request!(atmos, ["ana"])
        @test result == true
        @test all(atmos.tmp .> 0.0)

        # Test invalid request
        result = AGNI.setpt.request!(atmos, ["invalid_command"])
        @test result == false

        # Test saturation request (requires gas in atmosphere)
        result = AGNI.setpt.request!(atmos, ["sat", "H2O"])
        @test result == true
    end

    @testset "saturation!" begin
        # Set initial profile
        AGNI.setpt.isothermal!(atmos, 300.0)

        # Apply saturation for water
        AGNI.setpt.saturation!(atmos, "H2O")

        # Check that temperatures are still reasonable
        @test all(atmos.tmp .> 0.0)
        @test all(atmos.tmpl .> 0.0)

        # Test with non-existent gas (should return without error)
        AGNI.setpt.saturation!(atmos, "NonExistentGas")
        @test all(atmos.tmp .> 0.0)

        # Test with custom dTdew
        AGNI.setpt.isothermal!(atmos, 300.0)
        AGNI.setpt.saturation!(atmos, "H2O"; dTdew=0.1)
        @test all(atmos.tmp .> 0.0)
    end

    @testset "analytic!" begin
        AGNI.setpt.analytic!(atmos)

        # Check that temperatures are set and reasonable
        @test all(atmos.tmp .> 0.0)
        @test all(atmos.tmpl .> 0.0)

        # Temperatures should be in a reasonable range
        @test all(atmos.tmp .< 2000.0)  # Not extremely hot
        @test all(atmos.tmp .> 50.0)    # Not extremely cold
    end

    atmosphere.deallocate!(atmos)  # Clean up

end
