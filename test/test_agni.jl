using Test
using AGNI
using Logging

const AGNI_CORE_ROOT = abspath(joinpath(dirname(abspath(@__FILE__)), "../"))
const AGNI_CORE_TEST_CFG = joinpath(AGNI_CORE_ROOT, "test", "test.toml")
const AGNI_CORE_OUT_DIR = joinpath(AGNI_CORE_ROOT, "out/")
const AGNI_CORE_SF = joinpath(AGNI_CORE_ROOT, "res", "spectral_files", "Dayspring", "16", "Dayspring.sf")
_base_cfg() = AGNI.open_config(AGNI_CORE_TEST_CFG)

@testset "agni_core" begin
    tmpdir = mktempdir()
    logpath = joinpath(tmpdir, "agni.log")

    # file logger writes formatted records
    logger = AGNI.make_logger(logpath; to_term=false)
    with_logger(logger) do
        @info "hello logger"
    end
    @test isfile(logpath)
    logtxt = read(logpath, String)
    @test occursin("INFO", logtxt)
    @test occursin("hello logger", logtxt)

    # existing log file is replaced, not appended to
    open(logpath, "w") do io
        write(io, "stale")
    end
    logger = AGNI.make_logger(logpath; to_term=false)
    with_logger(logger) do
        @warn "fresh logger"
    end
    logtxt = read(logpath, String)
    @test !occursin("stale", logtxt)
    @test occursin("WARN", logtxt)

    # setup logging branches
    old_logger = current_logger()
    try
        @test AGNI.setup_logging(logpath, 0) === nothing
        @test AGNI.setup_logging(logpath, 2) === nothing
    finally
        global_logger(old_logger)
    end

    # open_config success case
    cfg_ok = joinpath(tmpdir, "cfg_ok.toml")
    out_ok = joinpath(tmpdir, "output_ok")
    write(cfg_ok, """
title = "unit-test"

[plots]
[planet]
[execution]
[physics]

[files]
output_dir = "$out_ok"
""")
    cfg = AGNI.open_config(cfg_ok)
    @test cfg["title"] == "unit-test"
    @test cfg["files"]["output_dir"] == out_ok

    # missing required header
    cfg_missing = joinpath(tmpdir, "cfg_missing.toml")
    write(cfg_missing, """
title = "unit-test"

[plots]
[planet]
[execution]

[files]
output_dir = "$out_ok"
""")
    @test_throws ErrorException AGNI.open_config(cfg_missing)

    # unsafe output directory
    cfg_unsafe = joinpath(tmpdir, "cfg_unsafe.toml")
    write(cfg_unsafe, """
title = "unit-test"

[plots]
[planet]
[execution]
[physics]

[files]
output_dir = "/"
""")
    @test_throws ErrorException AGNI.open_config(cfg_unsafe)
    rm(tmpdir; force=true, recursive=true)
end

@testset "test_config_parser" begin

    # check the config has expected values
    cfg = _base_cfg()
    @test haskey(cfg, "planet")
    @test haskey(cfg["planet"], "tmp_surf")
    @test haskey(cfg["composition"], "vmr_dict")
    @test haskey(cfg["physics"], "chemistry")
    @test cfg["physics"]["chemistry"] == false

    # hide error messages for these since we're testing that they throw
    with_logger(MinLevelLogger(current_logger(), Test.Logging.Error+1)) do
        # missing required key
        cfg = _base_cfg()
        delete!(cfg["planet"], "tmp_surf")
        @test AGNI.run_from_config(cfg) == false

        # overspecified gravity + mass
        cfg = _base_cfg()
        cfg["planet"]["mass"] = 5.972e24
        @test AGNI.run_from_config(cfg) == false

        # greybody without albedo_s
        cfg = _base_cfg()
        delete!(cfg["planet"], "albedo_s")
        @test AGNI.run_from_config(cfg) == false

        # p_surf with no composition source
        cfg = _base_cfg()
        delete!(cfg["composition"], "vmr_dict")
        @test AGNI.run_from_config(cfg) == false

        # composition overspecified (vmr_dict + vmr_file)
        cfg = _base_cfg()
        cfg["composition"]["vmr_file"] = "/tmp/vmr.csv"
        @test AGNI.run_from_config(cfg) == false

        # metallicities require chemistry
        cfg = _base_cfg()
        delete!(cfg["composition"], "vmr_dict")
        cfg["composition"]["metallicities"] = Dict("C" => 1.0)
        cfg["physics"]["chemistry"] = false
        @test AGNI.run_from_config(cfg) == false

        # transparent mode incompatible with chemistry
        cfg = _base_cfg()
        cfg["composition"]["transparent"] = true
        cfg["physics"]["chemistry"] = true
        @test AGNI.run_from_config(cfg) == false

        # partial pressures cannot be combined with VMR definition
        cfg = _base_cfg()
        delete!(cfg["composition"], "p_surf")
        cfg["composition"]["p_dict"] = Dict("N2" => 1.0e5)
        cfg["composition"]["vmr_dict"] = Dict("N2" => 1.0)
        @test AGNI.run_from_config(cfg) == false

        # greygas requires explicit opacities
        cfg = _base_cfg()
        cfg["files"]["input_sf"] = "greygas"
        @test AGNI.run_from_config(cfg) == false

        # invalid plot extension
        cfg = _base_cfg()
        cfg["plots"]["extension"] = "invalidext"
        @test AGNI.run_from_config(cfg) == false

        # latent heat requires rainout
        cfg = _base_cfg()
        cfg["physics"]["latent_heat"] = true
        cfg["physics"]["rainout"] = false
        @test AGNI.run_from_config(cfg) == false

        # sensible heat needs roughness + wind_speed
        cfg = _base_cfg()
        delete!(cfg["planet"], "roughness")
        @test AGNI.run_from_config(cfg) == false

        # sol_type=2 requires conductive-skin parameters
        cfg = _base_cfg()
        cfg["execution"]["solution_type"] = 2
        delete!(cfg["planet"], "skin_k")
        @test AGNI.run_from_config(cfg) == false

        # sol_type=3 requires flux_int
        cfg = _base_cfg()
        cfg["execution"]["solution_type"] = 3
        delete!(cfg["planet"], "flux_int")
        @test AGNI.run_from_config(cfg) == false

        # sol_type=4 requires target_olr
        cfg = _base_cfg()
        cfg["execution"]["solution_type"] = 4
        @test AGNI.run_from_config(cfg) == false
    end
end

# Greygas RT rejects scattering flags while accepting boundary opacity inputs.
@testset "greygas_config" begin
    mkpath(AGNI_CORE_OUT_DIR)
    κ_lw = 1.0e-4
    κ_sw = 0.0

    atmos_ok = atmosphere.Atmos_t()
    ok = atmosphere.setup!(atmos_ok, AGNI_CORE_ROOT, AGNI_CORE_OUT_DIR,
                            "  Greygas ",
                            1200.0, 1.0, 0.0, 0.0,
                            350.0,
                            10.0, 1.0e7,
                            30, 10.0, 1e-6,
                            Dict("N2" => 1.0), "";
                            real_gas=false,
                            thermo_functions=false,
                            flag_rayleigh=false,
                            flag_cloud=false,
                            flag_aerosol=false,
                            κ_grey_lw=κ_lw,
                            κ_grey_sw=κ_sw)
    @test ok
    @test atmos_ok.rt_scheme == atmosphere.RT_GREYGAS
    @test atmos_ok.SOCRATES_VERSION == "0000"
    @test isapprox(atmos_ok.κ_grey_lw, κ_lw; rtol=1e-12, atol=0.0)
    @test isapprox(atmos_ok.κ_grey_sw, κ_sw; rtol=0.0, atol=1e-12)
    @test atmos_ok.κ_grey_lw > atmos_ok.κ_grey_sw

    atmos_bad = atmosphere.Atmos_t()
    bad = with_logger(MinLevelLogger(current_logger(), Test.Logging.Error+1)) do
        atmosphere.setup!(atmos_bad, AGNI_CORE_ROOT, AGNI_CORE_OUT_DIR,
                        "greygas",
                        1200.0, 1.0, 0.0, 0.0,
                        350.0,
                        10.0, 1.0e7,
                        30, 10.0, 1e-6,
                        Dict("N2" => 1.0), "";
                        real_gas=false,
                        thermo_functions=false,
                        flag_rayleigh=true,
                        flag_cloud=false,
                        flag_aerosol=false,
                        κ_grey_lw=κ_lw,
                        κ_grey_sw=κ_sw)
    end
    @test !bad
    @test atmos_bad.rt_scheme == atmosphere.RT_GREYGAS
end

