using Test
using AGNI

@testset "blake" begin

    # Paths
    test_file = joinpath(AGNI.paths.RES_DIR, "config", "default.toml")
    temp_file = tempname() * "_blake.txt"

    # Test core Blake2b implementation against known hash
    hash_txt = "I was angry with my friend: I told my wrath, my wrath did end. I was angry with my foe: I told it not, my wrath did grow."
    hash_exp = "e4fcad40e634f623934d4b36255cb6c6341ec9e98d5abdf738ded677e0021f50f82c4bb09c23479d0253a9a400d3fd07d1f6b565469311c02a3f5819d199bf75"
    hash_obs = AGNI.blake.hash_string(hash_txt)
    @test hash_obs == hash_exp

    # Test hash_file with a file that exists
    #     Use one of the config files that should be present
    @test isfile(test_file)
    hash_result = AGNI.blake.hash_file(test_file)
    @test length(hash_result) == 128  # BLAKE2b produces 512-bit = 64 byte = 128 hex chars
    @test occursin(r"^[0-9a-f]{128}$", hash_result)  # hex string

    # Test hash_file with non-existent file
    hash_result = AGNI.blake.hash_file(temp_file)
    @test startswith(hash_result, "FILE_NOT_FOUND:")

    # Test valid_file with missing .chk file
    #     Create a temporary file without a .chk
    write(temp_file, "test content")
    validity = AGNI.blake.valid_file(temp_file)
    @test validity == false  # should fail because .chk is missing
    rm(temp_file, force=true)

    # Test valid_file with matching hash
    #     Create a temp file and its .chk file with correct hash
    write(temp_file, "test content")
    correct_hash = AGNI.blake.hash_file(temp_file)
    write(temp_file * ".chk", correct_hash)
    validity2 = AGNI.blake.valid_file(temp_file)
    @test validity2 == true
    rm(temp_file, force=true)
    rm(temp_file * ".chk", force=true)

    # Test valid_file with mismatched hash
    write(temp_file, "another test content")
    wrong_hash = "0" ^ 128  # all zeros
    write(temp_file * ".chk", wrong_hash)
    validity = AGNI.blake.valid_file(temp_file)
    @test validity == false
    rm(temp_file, force=true)
    rm(temp_file * ".chk", force=true)

end
