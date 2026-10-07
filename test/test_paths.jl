# Tests for src/interface/paths.jl. A res folder resolves to AGNI_DIR_<name>, else to the
# res root: AGNI_DIR_res, else the configured res, else the res/ of AGNI. Blank variables
# count as unset; "out" stays in the AGNI root; unsafe directories are refused.
using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)), "../"))

@testset "paths" begin
    unset = [("AGNI_DIR_" * n) => nothing for n in ("res", paths.RES_NAMES...)]

    # With no override, every res folder is the res/ of AGNI
    @testset "get_dir" begin
        withenv(unset...) do
            for name in ("thermodynamics", "scattering", "refractive", "config",
                         "stellar_spectra", "spectral_files", "blobs")
                @test paths.get_dir(name) == joinpath(paths.RES_DIR, name)
            end
        end
        @test paths.get_dir("out") == joinpath(paths.ROOT_DIR, "out")
        @test isnothing(paths.get_dir("does_not_exist"))
    end

    # Precedence AGNI_DIR_<name> > AGNI_DIR_res > configured res > res/; blanks count as
    # unset and a trailing separator is dropped
    @testset "get_dir_overrides" begin
        mktempdir() do tmp
            env, cfg, own = joinpath.(tmp, ("env", "cfg", "own"))
            nk, res = "AGNI_DIR_refractive", "AGNI_DIR_res"
            # (variables, configured res, expected res root, refractive folder if it moved)
            for (vars, cfg_res, root, moved) in (
                    ([], cfg, cfg, nothing),
                    ([res => env], nothing, env, nothing),
                    ([res => env], cfg, env, nothing),
                    ([nk => own * "/"], cfg, cfg, own),
                    ([nk => own, res => env], nothing, env, own),
                    ([res => " ", nk => ""], cfg, cfg, nothing))
                withenv(unset..., vars...) do
                    for name in paths.RES_NAMES
                        want = (name == "refractive" && !isnothing(moved)) ? moved :
                               joinpath(root, name)
                        @test paths.get_dir(name; res=cfg_res) == want
                    end
                end
            end
        end
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
