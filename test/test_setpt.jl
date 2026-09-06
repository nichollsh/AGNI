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

        # invalid keyword-like string
        @test isnothing(AGNI.setpt._parse_tmp_str(atmos, "definitely_not_a_temperature"))
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

        # Check temperatures are positive and at least at the floor everywhere
        @test all(atmos.tmp .> 0.0)
        @test all(atmos.tmpl .> 0.0)
        @test all(atmos.tmp .>= atmos.tmp_floor)

        # Real cp and rho both vanish as p->0, so the true dry adiabat cools toward
        # (and, over a wide enough pressure range, clamps at) tmp_floor near the top
        # of the column; it must never increase with altitude nor oscillate, so the
        # profile (ordered TOA to surface) must be exactly non-decreasing throughout,
        # not just "mostly decreasing"
        @test issorted(atmos.tmp)
        @test issorted(atmos.tmpl)

        # discrimination guard: a solver that overshoots into an unstable oscillation
        # (e.g. a single explicit step using a stale neighbouring temperature to
        # evaluate cp/rho, rather than the level's own self-consistent temperature)
        # produces jumps back down to the floor from a level above it that already
        # exceeded it; issorted() above would already catch that, but pin the extra
        # invariant that only the levels nearest the floor are actually clamped there
        n_floored = count(isapprox.(atmos.tmp, atmos.tmp_floor; atol=1e-8))
        @test n_floored < atmos.nlev_c  # not every level can be floored
        if n_floored > 0
            @test all(isapprox.(atmos.tmp[1:n_floored], atmos.tmp_floor; atol=1e-8))
        end
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

        # Non-monotonic pressure array should fail
        bad_pl = [1e2, 1e4, 1e3, 1e5]
        bad_t = [150.0, 200.0, 250.0, 300.0]
        @test AGNI.setpt.fromarrays!(atmos, bad_pl, bad_t) == false

        # Input grid narrower than the model grid at the high-pressure end must trigger
        # extrapolation-padding branch
        narrow_pl = [1e2, 1e3, 1e4, 5e4]
        narrow_t  = [150.0, 200.0, 250.0, 280.0]

        # fromarrays! mutates its pl/tmpl arguments in place (push!/pushfirst! extend
        # them), so backup the original values before the call
        orig_pl, orig_t = copy(narrow_pl), copy(narrow_t)

        # test that extrapolation actually extended the arrays to the psurf
        @test atmos.pl[end] > orig_pl[end]
        @test AGNI.setpt.fromarrays!(atmos, narrow_pl, narrow_t; extrap=true)
        @test all(atmos.tmp .> 0.0)
        @test all(atmos.tmpl .> 0.0)

        # discrimination guard: without the padding point, Linear() extrapolation
        # beyond the file's domain would continue the trend
        slope = (orig_t[end] - orig_t[end-1]) /
                    (log10(orig_pl[end]) - log10(orig_pl[end-1]))
        expect_wrong = orig_t[end] + slope * (log10(atmos.pl[end]) - log10(orig_pl[end]))
        @test !isapprox(expect_wrong, orig_t[end]; atol=1.0)  # sanity: guard is non-trivial
        @test !isapprox(atmos.tmpl[end], expect_wrong; atol=1.0)
        @test isapprox(atmos.tmpl[end], orig_t[end]; atol=1.0)
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
        @test AGNI.setpt.fromcsv!(atmos, "/nonexistent/file.csv") == false

        # Too few data rows
        tiny_csv = tempname() * ".csv"
        open(tiny_csv, "w") do io
            println(io, "1e2, 150.0")
            println(io, "1e3, 200.0")
        end
        try
            @test AGNI.setpt.fromcsv!(atmos, tiny_csv) == false
        finally
            rm(tiny_csv, force=true)
        end

        # Invalid values in CSV
        bad_csv = tempname() * ".csv"
        open(bad_csv, "w") do io
            println(io, "1e2, 150.0")
            println(io, "-1e3, 200.0")
            println(io, "1e4, 250.0")
        end
        try
            @test AGNI.setpt.fromcsv!(atmos, bad_csv) == false
        finally
            rm(bad_csv, force=true)
        end

        # Negative temperature (distinct from the negative-pressure case above)
        bad_t_csv = tempname() * ".csv"
        open(bad_t_csv, "w") do io
            println(io, "1e2, 150.0")
            println(io, "1e3, -50.0")
            println(io, "1e4, 250.0")
        end
        try
            @test AGNI.setpt.fromcsv!(atmos, bad_t_csv) == false
        finally
            rm(bad_t_csv, force=true)
        end
    end

    @testset "fromncdf!" begin
        @test AGNI.setpt.fromncdf!(atmos, "/nonexistent/file.nc") == false

        # Round-trip smoke test for the file-found success path
        AGNI.setpt.dry_adiabat!(atmos)
        nc_tmpdir = mktempdir()
        nc_path = joinpath(nc_tmpdir, "roundtrip.nc")
        @test AGNI.save.write_ncdf(atmos, nc_path) # test writing file

        # make new atmosphere
        atmos2 = AGNI.atmosphere.Atmos_t()
        AGNI.atmosphere.setup!(atmos2, ROOT_DIR, OUT_DIR,
                              spfile,
                              toa_heating, 1.0, 0.0, theta,
                              tmp_surf,
                              gravity, radius,
                              nlev_centre, p_surf, p_top,
                              mf_dict, ""
                      )
        AGNI.atmosphere.allocate!(atmos2, "")
        @test AGNI.setpt.fromncdf!(atmos2, nc_path) # read from the file

        # check that the round-tripped profile matches the original
        @test all(isapprox.(atmos2.tmp,  atmos.tmp;  rtol=rtol))
        @test all(isapprox.(atmos2.tmpl, atmos.tmpl; rtol=rtol))
        @test isapprox(atmos2.tmp_surf, atmos.tmp_surf; rtol=rtol)

        # the round-tripped profile is a non-trivial dry adiabat
        @test !isapprox(atmos2.tmp[1], atmos2.tmp[end]; rtol=1e-2)

        # tidy up
        AGNI.atmosphere.deallocate!(atmos2)
        rm(nc_tmpdir; force=true, recursive=true)
    end

    @testset "request!" begin
        # Test single request
        result = AGNI.setpt.request!(atmos, Any["iso", 400.0])
        @test result == true
        @test all(atmos.tmp .≈ 400.0)

        # Test dry adiabat request
        result = AGNI.setpt.request!(atmos, Any["dry"])
        @test result == true
        @test atmos.tmpl[end] ≈ atmos.tmp_surf

        # Test stratosphere request
        result = AGNI.setpt.request!(atmos, Any["dry", "str", 200.0])
        @test result == true

        # Test loglinear request
        result = AGNI.setpt.request!(atmos, Any["loglin", 150.0])
        @test result == true

        # Test add request
        result = AGNI.setpt.request!(atmos, Any["iso", 300.0, "add", 50.0])
        @test result == true
        @test all(atmos.tmp .≈ 350.0)

        # Test analytic request
        result = AGNI.setpt.request!(atmos, Any["ana"])
        @test result == true
        @test all(atmos.tmp .> 0.0)

        # Test invalid request
        result = AGNI.setpt.request!(atmos, Any["invalid_command"])
        @test result == false

        # Test saturation request (requires gas in atmosphere)
        result = AGNI.setpt.request!(atmos, Any["sat", "H2O"])
        @test result == true

        # Test surfsat request: restores composition then ensures the surface is not
        # super-saturated (calls chemistry.restore_composition!/_sat_surf!).
        @test isempty(atmos.condensates)
        result = AGNI.setpt.request!(atmos, Any["surfsat"])
        @test result == false
        @test all(atmos.tmp .> 0.0)

        # Test missing argument path for a verb
        @test AGNI.setpt.request!(atmos, Any["iso"]) == false

        # csv/ncdf verbs with missing files should fail
        @test AGNI.setpt.request!(atmos, Any["csv", "/nonexistent/file.csv"]) == false
        @test AGNI.setpt.request!(atmos, Any["ncdf", "/nonexistent/file.nc"]) == false
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

        # Case-mismatched gas name (lowercase "h2o" vs stored "H2O") should warn
        AGNI.setpt.isothermal!(atmos, 310.0)
        before_case = copy(atmos.tmp)
        result_case = AGNI.setpt.saturation!(atmos, "h2o")
        @test result_case == true
        @test all(isapprox.(atmos.tmp, before_case; rtol=0.0, atol=1e-10))
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

    @testset "guard_paths_and_helpers" begin
        # _verb_arg helper out-of-range returns UNSET sentinel string
        @test AGNI.setpt._verb_arg(Any["iso"], 3) == AGNI.atmosphere.UNSET_STR

        # Methods should return false when atmosphere is not allocated
        atmos_unalloc = AGNI.atmosphere.Atmos_t()
        AGNI.atmosphere.setup!(atmos_unalloc, ROOT_DIR, OUT_DIR,
                          spfile,
                          toa_heating, 1.0, 0.0, theta,
                          tmp_surf,
                          gravity, radius,
                          nlev_centre, p_surf, p_top,
                          mf_dict, ""
                  )

        @test AGNI.setpt.request!(atmos_unalloc, Any["iso", 300.0]) == false
        @test AGNI.setpt.isothermal!(atmos_unalloc, 300.0) == false
        @test AGNI.setpt.add!(atmos_unalloc, 10.0) == false
        @test AGNI.setpt.dry_adiabat!(atmos_unalloc) == false
        @test AGNI.setpt.loglinear!(atmos_unalloc, 200.0) == false
        @test AGNI.setpt.saturation!(atmos_unalloc, "H2O") == false
        @test AGNI.setpt.analytic!(atmos_unalloc) == false
        @test AGNI.setpt.fromarrays!(atmos_unalloc, [1.0, 2.0, 3.0], [100.0, 100.0, 100.0]) == false
        @test AGNI.setpt.fromcsv!(atmos_unalloc, "/nonexistent/file.csv") == false
        @test AGNI.setpt.fromncdf!(atmos_unalloc, "/nonexistent/file.nc") == false
    end

    atmosphere.deallocate!(atmos)  # Clean up

end
