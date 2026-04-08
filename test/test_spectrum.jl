using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
RES_DIR = joinpath(ROOT_DIR,"res/")

temp_sf = tempname() * ".sf"

@testset "spectrum" begin

    # Test get_socrates_version
    @testset "socrates_version" begin
        version = AGNI.spectrum.get_socrates_version()
        @test !isempty(version)
        @test version isa String
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
        @test all(fl .>= AGNI.spectrum.FLOAT_SML)
        @test all(fl .<= AGNI.spectrum.FLOAT_BIG)
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

end
