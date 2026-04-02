using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)),"../"))
RES_DIR = joinpath(ROOT_DIR,"res/")
OUT_DIR = joinpath(ROOT_DIR,"out/")

@testset "blake" begin

    # Test hash_file with a file that exists
    # Use one of the config files that should be present
    test_file = joinpath(RES_DIR, "config", "default.toml")
    if isfile(test_file)
        hash_result = AGNI.phys.blake.hash_file(test_file)
        @test length(hash_result) == 128  # BLAKE2b produces 512-bit = 64 byte = 128 hex chars
        @test occursin(r"^[0-9a-f]{128}$", hash_result)  # hex string
    end

    # Test hash_file with non-existent file
    nonexistent_file = joinpath(OUT_DIR, "this_file_does_not_exist_xyz123.tmp")
    hash_result = AGNI.phys.blake.hash_file(nonexistent_file; quiet=true)
    @test startswith(hash_result, "FILE_NOT_FOUND:")

    # Test valid_file with missing .chk file
    # Create a temporary file without a .chk
    temp_file = joinpath(OUT_DIR, "test_blake_temp.txt")
    write(temp_file, "test content")
    validity = AGNI.phys.blake.valid_file(temp_file; quiet=true)
    @test validity == false  # should fail because .chk is missing
    rm(temp_file, force=true)

    # Test valid_file with matching hash
    # Create a temp file and its .chk file with correct hash
    temp_file2 = joinpath(OUT_DIR, "test_blake_temp2.txt")
    write(temp_file2, "test content for validation")
    correct_hash = AGNI.phys.blake.hash_file(temp_file2; quiet=true)
    write(temp_file2 * ".chk", correct_hash)
    validity2 = AGNI.phys.blake.valid_file(temp_file2; quiet=true)
    @test validity2 == true
    rm(temp_file2, force=true)
    rm(temp_file2 * ".chk", force=true)

    # Test valid_file with mismatched hash
    temp_file3 = joinpath(OUT_DIR, "test_blake_temp3.txt")
    write(temp_file3, "another test content")
    wrong_hash = "0" ^ 128  # all zeros
    write(temp_file3 * ".chk", wrong_hash)
    validity3 = AGNI.phys.blake.valid_file(temp_file3; quiet=true)
    @test validity3 == false
    rm(temp_file3, force=true)
    rm(temp_file3 * ".chk", force=true)

end
