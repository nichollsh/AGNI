using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)), "../"))

@testset "paths" begin
    @testset "get_dir" begin
        @test paths.get_dir("thermodynamics") == joinpath(paths.RES_DIR, "thermodynamics")
        @test paths.get_dir("scattering") == joinpath(paths.RES_DIR, "scattering")
        withenv("AGNI_REFRACTIVE_DIR" => nothing) do
            @test paths.get_dir("refractive") == joinpath(paths.RES_DIR, "refractive")
        end
        @test paths.get_dir("config") == joinpath(paths.RES_DIR, "config")
        @test paths.get_dir("stellar_spectra") == joinpath(paths.RES_DIR, "stellar_spectra")
        @test paths.get_dir("spectral_files") == joinpath(paths.RES_DIR, "spectral_files")
        @test paths.get_dir("blobs") == joinpath(paths.RES_DIR, "blobs")
        @test paths.get_dir("out") == joinpath(paths.ROOT_DIR, "out")
        @test isnothing(paths.get_dir("does_not_exist"))
    end

    @testset "refractive_dir_override" begin
        default = joinpath(paths.RES_DIR, "refractive")
        tmp = mktempdir()
        material = first(AGNI.density.list_condensate_rho())
        write(joinpath(tmp, material * ".txt"), "# test\n1.0 1.5 0.01\n2.0 1.4 0.02\n")
        withenv("AGNI_REFRACTIVE_DIR" => tmp) do
            @test paths.get_dir("refractive") == tmp
            @test AGNI.aerosol_optics.nk_path(material) == joinpath(tmp, material * ".txt")
            @test AGNI.aerosol_optics.list_materials() == [material]
        end
        withenv("AGNI_REFRACTIVE_DIR" => "") do
            @test paths.get_dir("refractive") == default
        end
        missing_dir = joinpath(tmp, "absent")
        withenv("AGNI_REFRACTIVE_DIR" => missing_dir) do
            @test_logs (:warn, r"not a directory") paths.get_dir("refractive")
            @test AGNI.aerosol_optics.list_materials() == String[]
        end
        rm(tmp; force=true, recursive=true)
    end

    @testset "is_safe_dir" begin
        # check explicit unsafe paths
        @test paths.is_safe_dir("") == false
        @test paths.is_safe_dir("/") == false
        @test paths.is_safe_dir(homedir()) == false
        @test paths.is_safe_dir(paths.ROOT_DIR) == false
        @test paths.is_safe_dir(paths.RES_DIR) == false
        @test paths.is_safe_dir(pwd()) == false

        # make a safe temp dir
        tmp_safe = mktempdir()
        tmp_git = mktempdir()
        @test paths.is_safe_dir(tmp_safe) == true

        # make it unsafe
        mkdir(joinpath(tmp_git, ".git"))
        @test paths.is_safe_dir(tmp_git) == false

        # tidy up
        rm(tmp_safe; force=true, recursive=true)
        rm(tmp_git; force=true, recursive=true)
    end

    @testset "constpaths" begin
        @test normpath(paths.ROOT_DIR) == normpath(ROOT_DIR)
        @test normpath(paths.RES_DIR) == normpath(joinpath(ROOT_DIR, "res"))
        @test normpath(paths.FWL_DATA) == normpath(joinpath(get(ENV, "FWL_DATA", paths.RES_DIR)))
    end

    @testset "get_avail_space" begin
        # Existing directory: must report a real, physically-plausible disk
        # quantity, cross-checked against an independent diskstat() call
        tmp_existing = mktempdir()
        avail = paths.get_avail_space(tmp_existing)
        @test avail isa Int64
        stat_direct = diskstat(tmp_existing)
        @test avail == Int64(stat_direct.available)

        # Available space cannot exceed total space, and cannot be negative.
        @test 0 <= avail <= stat_direct.total
        rm(tmp_existing; force=true, recursive=true)

        # Nonexistent path: get_avail_space() must fall back to querying "/"
        nonexistent = joinpath(tmp_existing, "does_not_exist", "deeper")
        @test !ispath(nonexistent)
        @test_throws Base.IOError diskstat(nonexistent)
        @test paths.get_avail_space(nonexistent) == paths.get_avail_space("/")
        @test paths.get_avail_space(nonexistent) > 0

        # Edge case: empty string path also triggers the same fallback
        @test !ispath("")
        @test paths.get_avail_space("") == paths.get_avail_space("/")
    end
end
