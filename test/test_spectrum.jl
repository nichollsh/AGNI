"""
Tests for `src/energy/spectrum.jl`.

Invariants exercised:
- Spectral file parsing and block insertion succeed with valid inputs.
- Aerosol headers and aerosol averaging produce non-empty outputs.
- Guard paths return gracefully on invalid inputs.
"""
const _SPECTRUM_TESTS_DOC = nothing
using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
RES_DIR = joinpath(ROOT_DIR,"res/")
RAD_DIR = AGNI.atmosphere.RAD_DIR

temp_sf = tempname() * ".sf"

# Use 16 to ensure the non-default precision is preserved rather than defaulted.
module DummySocratesPrecision
    const SOCRATES_REAL_BYTES = 4
end

# Missing precision simulates older SOCRATES Julia builds.
module DummySocratesNoPrecision
end


@testset "spectrum" begin

    # Check SOCRATES metainfo
    @testset "socrates_meta" begin

        # Check socrates was found
        @test isfile(joinpath(RAD_DIR,"version"))

        # Check SOCRATES precision getter
        @testset "socrates_precision" begin

            # test defined as single
            precision = AGNI.spectrum.get_socrates_precision(DummySocratesPrecision)
            @test !isempty(precision)
            @test precision == "single"

            # test fallback to double
            fallback = AGNI.spectrum.get_socrates_precision(DummySocratesNoPrecision)
            @test !isempty(fallback)
            @test fallback == "double"
        end

        # Test get_socrates_version
        @testset "socrates_version" begin
            version = AGNI.spectrum.get_socrates_version(RAD_DIR)
            @test version isa String
            @test startswith(version, "2") # this millenium
        end
    end

    # Test count_gases with a real spectral file
    @testset "count_gases_valid" begin
        # Use one of the test spectral files
        spfile = joinpath(RES_DIR, "spectral_files", "Dayspring", "48", "Dayspring.sf")
        @test isfile(spfile)
        num_gases = AGNI.spectrum.count_gases(spfile)
        @test num_gases > 0
        @test num_gases isa Int
    end

    # Test count_gases with non-existent file
    @testset "count_gases_missing" begin
        rm(temp_sf, force=true)  # Ensure the file does not exist
        num_gases = AGNI.spectrum.count_gases(temp_sf)
        @test num_gases == -1
    end

    # Test count_gases with invalid file (no gas count line)
    @testset "count_gases_invalid" begin
        # Create a temporary file with invalid content
        write(temp_sf, "This is not a valid spectral file\nNo gas information here\n")
        num_gases = AGNI.spectrum.count_gases(temp_sf)
        @test num_gases == -1
        rm(temp_sf, force=true)
    end

    # Test count_gases with malformed gas count line (wrong format)
    @testset "count_gases_malformed_line" begin
        write(temp_sf, "Total number of gaseous absorbers not_a_number\n")
        num_gases = AGNI.spectrum.count_gases(temp_sf)
        @test num_gases == -1
        rm(temp_sf, force=true)
    end

    # Test count_gases with wrong split count
    @testset "count_gases_wrong_split" begin
        write(temp_sf, "Total number of gaseous absorbers = 5 = extra\n")
        num_gases = AGNI.spectrum.count_gases(temp_sf)
        @test num_gases == -1
        rm(temp_sf, force=true)
    end

    # Test count_gases with unparseable number (to hit catch block)
    @testset "count_gases_unparseable" begin
        write(temp_sf, "Total number of gaseous absorbers = not_a_number\n")
        num_gases = AGNI.spectrum.count_gases(temp_sf)
        @test num_gases == -1
        rm(temp_sf, force=true)
    end

    # Test insert_aerosol_header
    @testset "insert_aerosol_header" begin
        # Create a minimal spectral file
        sf_content = """
Line 1
Line 2
Line 3
*END
Some more content
*BLOCK: TYPE =    1
"""
        write(temp_sf, sf_content)

        # Try inserting aerosol header with a valid aerosol name
        success = AGNI.spectrum.insert_aerosol_header(temp_sf, ["dust"])
        @test success

        # Check that file was modified
        content = read(temp_sf, String)
        @test contains(content, "Total number of aerosols")
        @test contains(content, "Dust-like Aerosol")

        rm(temp_sf, force=true)
    end

    # Test insert_aerosol_header with existing aerosol data
    @testset "insert_aerosol_header_existing" begin
        write(temp_sf, "This file already has aerosols in it\n*END\n")
        success = AGNI.spectrum.insert_aerosol_header(temp_sf, ["dust"])
        @test !success
        rm(temp_sf, force=true)
    end

    # Test insert_aerosol_header with missing *END marker
    @testset "insert_aerosol_header_no_end" begin
        write(temp_sf, "Line 1\nLine 2\nNo END marker here\n")
        success = AGNI.spectrum.insert_aerosol_header(temp_sf, ["sulphuric_acid"])
        @test !success
        rm(temp_sf, force=true)
    end

    # Test blackbody_star
    @testset "blackbody_star" begin
        Teff = 5800.0  # Sun-like star
        S0 = 1361.0     # Solar constant

        wl, fl = AGNI.spectrum.blackbody_star(Teff, S0)

        @test length(wl) == length(fl)
        @test length(wl) > 0
        @test all(fl .>= AGNI.spectrum.SMALLFLOAT)
        @test all(fl .<= AGNI.spectrum.BIGFLOAT)
        @test wl[1] == 1.0
        @test wl[end] == 100e3
    end

    # Test load_from_file with valid file
    @testset "load_from_file_valid" begin
        # Create a test stellar spectrum file
        test_data = """
# Header line 1
# Header line 2
100.0 1.5e10
200.0 2.3e10
300.0 1.8e10
"""
        write(temp_sf, test_data)
        wl, fl = AGNI.spectrum.load_from_file(temp_sf)

        @test length(wl) == 3
        @test length(fl) == 3
        @test wl[1] == 100.0
        @test fl[1] == 1.5e10
        rm(temp_sf, force=true)
    end

    # Test write_to_socrates_format
    @testset "write_to_socrates_format" begin
        # Create test wavelength and flux arrays
        wl = collect(range(100.0, 1000.0, length=1000))
        fl = ones(Float64, 1000) .* 1e10

        temp_star = tempname() * ".txt"
        success = AGNI.spectrum.write_to_socrates_format(wl, fl, temp_star, 500)
        @test success
        @test isfile(temp_star)

        # Check file content
        content = read(temp_star, String)
        @test contains(content, "Star spectrum at TOA")
        @test contains(content, "*BEGIN_DATA")
        @test contains(content, "*END")

        rm(temp_star, force=true)
    end

    # Test write_to_socrates_format with mismatched arrays
    @testset "write_to_socrates_format_mismatch" begin
        wl = collect(range(100.0, 1000.0, length=1000))
        fl = ones(Float64, 999)  # Different length

        temp_star = tempname() * ".txt"
        success = AGNI.spectrum.write_to_socrates_format(wl, fl, temp_star)
        @test !success
        rm(temp_star, force=true)
    end

    # Test write_to_socrates_format with short spectrum
    @testset "write_to_socrates_format_short" begin
        wl = collect(range(100.0, 1000.0, length=100))
        fl = ones(Float64, 100) .* 1e10

        temp_star = tempname() * ".txt"
        success = AGNI.spectrum.write_to_socrates_format(wl, fl, temp_star)
        @test success
        rm(temp_star, force=true)
    end

    # Test write_to_socrates_format with invalid wavelength (too small)
    @testset "write_to_socrates_format_wl_too_small" begin
        wl = [1e-50, 100.0, 200.0]
        fl = ones(Float64, 3) .* 1e10

        temp_star = tempname() * ".txt"
        success = AGNI.spectrum.write_to_socrates_format(wl, fl, temp_star)
        @test !success
        rm(temp_star, force=true)
    end

    # Test write_to_socrates_format with descending wavelength array
    @testset "write_to_socrates_format_descending" begin
        wl = [1000.0, 500.0, 100.0]  # Descending order
        fl = ones(Float64, 3) .* 1e10

        temp_star = tempname() * ".txt"
        success = AGNI.spectrum.write_to_socrates_format(wl, fl, temp_star)
        @test !success
        rm(temp_star, force=true)
    end

    # Test write_to_socrates_format with duplicate wavelengths
    @testset "write_to_socrates_format_duplicates" begin
        wl = [100.0, 200.0, 200.0, 300.0]  # Has duplicate
        fl = ones(Float64, 4) .* 1e10

        temp_star = tempname() * ".txt"
        success = AGNI.spectrum.write_to_socrates_format(wl, fl, temp_star)
        @test success  # Should succeed after removing duplicates
        rm(temp_star, force=true)
    end

    # Clean up any remaining temp files
    rm(temp_sf, force=true)


    @testset "aerosol_guard_invalid_input" begin
        tmpdir = mktempdir()
        orig = joinpath(tmpdir, "base.sf")
        star = joinpath(tmpdir, "star.dat")
        outp = joinpath(tmpdir, "out.sf")

        # insert_blocks: missing original file
        @test !AGNI.spectrum.insert_blocks( RAD_DIR,
            orig, star, outp, false, false; aerosol_avg_files=Dict{String,String}()
        )

        # write minimal original + _k to pass cp step, but missing star file should fail
        write(orig, "Total number of gaseous absorbers = 1\n*END\n")
        write(orig * "_k", "dummy k-table\n")
        @test !AGNI.spectrum.insert_blocks( RAD_DIR,
            orig, star, outp, false, false; aerosol_avg_files=Dict{String,String}()
        )

        # write star file, but malformed spectral file => gas count parse failure
        write(star, "header\nheader\n1.0 1.0\n")
        write(orig, "No gas count line here\n*END\n")
        @test !AGNI.spectrum.insert_blocks( RAD_DIR,
            orig, star, outp, false, false; aerosol_avg_files=Dict{String,String}()
        )

        # valid gas-count line, but missing prep binaries / execution path should still be handled
        write(orig, "Total number of gaseous absorbers = 1\n*END\n")
        @test !AGNI.spectrum.insert_blocks( RAD_DIR,
            orig, star, outp, true, false; aerosol_avg_files=Dict{String,String}()
        )

        # generate_aerosol_avg_files guard paths
        missing_orig = joinpath(tmpdir, "missing.sf")
        @test isempty(AGNI.spectrum.generate_aerosol_avg_files( RAD_DIR,
            missing_orig, ["dust"], tmpdir, 2, star, tmpdir
        ))

        @test isempty(AGNI.spectrum.generate_aerosol_avg_files(RAD_DIR,
            orig, ["dust"], tmpdir, 0, star, tmpdir
        ))

        @test isempty(AGNI.spectrum.generate_aerosol_avg_files(RAD_DIR,
            orig, ["dust"], tmpdir, 2, "", tmpdir
        ))

        @test isempty(AGNI.spectrum.generate_aerosol_avg_files(RAD_DIR,
            orig, ["dust"], tmpdir, 2, joinpath(tmpdir, "missing_star.dat"), tmpdir
        ))

        # no species => early empty return without external tool
        @test isempty(AGNI.spectrum.generate_aerosol_avg_files(RAD_DIR,
            orig, String[], tmpdir, 2, star, tmpdir
        ))

        # species provided but missing .mon file
        @test isempty(AGNI.spectrum.generate_aerosol_avg_files(RAD_DIR,
            orig, ["dust"], tmpdir, 2, star, tmpdir
        ))

        rm(tmpdir; force=true, recursive=true)
    end

    # Check successful data insersion
    @testset "aerosol_insertion_success" begin

        # Setup paths
        tmpdir = mktempdir()
        spfile = joinpath(RES_DIR, "spectral_files", "Dayspring", "48", "Dayspring.sf")

        # A real spectral file
        @test isfile(spfile)
        @test isfile(spfile * "_k")

        # Some SOCRATES scattering data is pre-computed
        scattering_dir = joinpath(RES_DIR, "scattering")
        @test isdir(scattering_dir)

        # Choose the 1,2 aerosol species with a matching .mon file in the repo data.
        available_species = [
            s for s in AGNI.spectrum.input_head_pcf.aerosol_suffix
            if isfile(joinpath(scattering_dir, s * ".mon"))
        ]
        @test !isempty(available_species)
        species = [available_species[1], available_species[2]] # first two

        # Get stellar spectrum
        star_file = joinpath(tmpdir, "star.dat")
        wl = collect(range(100.0, 1000.0, length=1000))
        fl = ones(Float64, 1000) .* 1e10
        @test AGNI.spectrum.write_to_socrates_format(wl, fl, star_file, 500)
        @test isfile(star_file)

        # Generate aerosol average-properties file(s)
        avg_files = AGNI.spectrum.generate_aerosol_avg_files(
            RAD_DIR, spfile, species, tmpdir, 1, star_file, scattering_dir
        )

        # should be 2
        @test length(avg_files) == length(species)

        # check output from input
        @test haskey(avg_files, species[1])
        avg_file = avg_files[species[1]]
        @test isfile(avg_file)
        @test filesize(avg_file) > 0
        @test length(readlines(avg_file)) > 1

        # insert data into spectral file
        outp_file = joinpath(tmpdir, "runtime.sf")
        success = AGNI.spectrum.insert_blocks(
            RAD_DIR, spfile, star_file, outp_file, false, true;
            aerosol_avg_files=avg_files
        )

        # check the spectral file
        @test success
        @test isfile(outp_file)
        @test isfile(outp_file * "_k")
        @test filesize(outp_file) > 0
        @test filesize(outp_file * "_k") > 0

        # check the spectral file contains aerosol data
        out_lines = readlines(outp_file)
        @test any(contains.(out_lines, "Total number of aerosols"))
        @test any(contains.(out_lines, "List of indexing numbers of aerosols"))

        # remove this folder
        rm(tmpdir; force=true, recursive=true)
    end

end
