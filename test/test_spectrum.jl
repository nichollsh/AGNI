using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
RES_DIR = joinpath(ROOT_DIR,"res/")

temp_sf = "/tmp/test_spectral_file_xyz123.sf"

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
        if isfile(spfile)
            num_gases = AGNI.spectrum.count_gases(spfile)
            @test num_gases > 0
            @test num_gases isa Int
        end
    end

    # Test count_gases with non-existent file
    @testset "count_gases_missing" begin
        rm(temp_sf, force=true)  # Ensure the file does not exist
        num_gases = AGNI.spectrum.count_gases(temp_sf; quiet=true)
        @test num_gases == -1
    end

    # Test count_gases with invalid file (no gas count line)
    @testset "count_gases_invalid" begin
        # Create a temporary file with invalid content
        write(temp_sf, "This is not a valid spectral file\nNo gas information here\n")
        num_gases = AGNI.spectrum.count_gases(temp_sf; quiet=true)
        @test num_gases == -1
        rm(temp_sf, force=true)
    end

end
