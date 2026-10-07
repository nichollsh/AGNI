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

    # Precedence AGNI_DIR_<name> > AGNI_DIR_res > configured res > res/; blank is unset
    @testset "get_dir_overrides" begin
        mktempdir() do tmp
            env_res, cfg_res = joinpath(tmp, "env"), joinpath(tmp, "cfg")
            own, env_nk = joinpath(tmp, "own"), joinpath(tmp, "env", "refractive")
            # (variables, configured res, root of the other folders, refractive folder)
            cases = [
                ([], cfg_res, cfg_res, joinpath(cfg_res, "refractive")),
                (["AGNI_DIR_res" => env_res], nothing, env_res, env_nk),
                (["AGNI_DIR_res" => env_res], cfg_res, env_res, env_nk),
                (["AGNI_DIR_refractive" => own], cfg_res, cfg_res, own),
                (["AGNI_DIR_refractive" => own, "AGNI_DIR_res" => env_res], nothing,
                    env_res, own),
                # blanks fall through, and a trailing separator is dropped
                (["AGNI_DIR_res" => "  ", "AGNI_DIR_refractive" => ""], cfg_res, cfg_res,
                    joinpath(cfg_res, "refractive")),
                (["AGNI_DIR_refractive" => own * "/"], nothing, paths.RES_DIR, own),
            ]
            for (vars, res, root, refractive) in cases
                withenv(unset..., vars...) do
                    @test paths.get_dir("refractive"; res=res) == refractive
                    for name in filter(!=("refractive"), paths.RES_NAMES)
                        @test paths.get_dir(name; res=res) == joinpath(root, name)
                    end
                end
            end
        end
    end

    # A refractive override reaches nk_path, list_materials and read_nk; n and k differ so a
    # swapped column would fail
    @testset "refractive_override_reaches_the_readers" begin
        mktempdir() do tmp
            material = first(AGNI.density.list_condensate_rho())
            write(joinpath(tmp, material * ".txt"), "# test\n1.0 1.5 0.01\n2.0 1.3 0.07\n")
            withenv(unset..., "AGNI_DIR_refractive" => tmp) do
                @test AGNI.aerosol_optics.list_materials() == [material]
                nk = AGNI.aerosol_optics.nk_path(material)
                λ, n, k = AGNI.aerosol_optics.read_nk(nk)
                @test λ ≈ [1.0e-6, 2.0e-6]
                @test n ≈ [1.5, 1.3]
                @test k ≈ [0.01, 0.07]
            end
            @test AGNI.aerosol_optics.list_materials(tmp) == [material]
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
