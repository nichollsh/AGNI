using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
RES_DIR = joinpath(ROOT_DIR,"res/")

test_file = joinpath(RES_DIR, "config", "default.toml")
temp_file = tempname() * "_blake.txt"

@testset "blake" begin

    # Test hash_file with a file that exists
    #     Use one of the config files that should be present
    @test isfile(test_file) 
    hash_result = AGNI.phys.blake.hash_file(test_file)
    @test length(hash_result) == 128  # BLAKE2b produces 512-bit = 64 byte = 128 hex chars
    @test occursin(r"^[0-9a-f]{128}$", hash_result)  # hex string

    # Test hash_file with non-existent file
    hash_result = AGNI.phys.blake.hash_file(temp_file; quiet=true)
    @test startswith(hash_result, "FILE_NOT_FOUND:")

    # Test valid_file with missing .chk file
    #     Create a temporary file without a .chk
    write(temp_file, "test content")
    validity = AGNI.phys.blake.valid_file(temp_file; quiet=true)
    @test validity == false  # should fail because .chk is missing
    rm(temp_file, force=true)

    # Test valid_file with matching hash
    #     Create a temp file and its .chk file with correct hash
    correct_hash = AGNI.phys.blake.hash_file(temp_file; quiet=true)
    write(temp_file * ".chk", correct_hash)
    validity2 = AGNI.phys.blake.valid_file(temp_file; quiet=true)
    @test validity2 == true
    rm(temp_file, force=true)
    rm(temp_file * ".chk", force=true)

    # Test valid_file with mismatched hash
    write(temp_file, "another test content")
    wrong_hash = "0" ^ 128  # all zeros
    write(temp_file * ".chk", wrong_hash)
    validity = AGNI.phys.blake.valid_file(temp_file; quiet=true)
    @test validity == false
    rm(temp_file, force=true)
    rm(temp_file * ".chk", force=true)

end
