using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)), "../"))

@testset "paths" begin
    @testset "get_dir" begin
        @test paths.get_dir("thermodynamics") == joinpath(paths.RES_DIR, "thermodynamics")
        @test paths.get_dir("scattering") == joinpath(paths.RES_DIR, "scattering")
        @test paths.get_dir("config") == joinpath(paths.RES_DIR, "config")
        @test paths.get_dir("stellar_spectra") == joinpath(paths.RES_DIR, "stellar_spectra")
        @test paths.get_dir("spectral_files") == joinpath(paths.RES_DIR, "spectral_files")
        @test paths.get_dir("blobs") == joinpath(paths.RES_DIR, "blobs")
        @test paths.get_dir("out") == joinpath(paths.ROOT_DIR, "out")
        @test isnothing(paths.get_dir("does_not_exist"))
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
end
