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

    # missing output_dir key entirely (distinct from the present-but-empty-string
    # case exercised implicitly by requiring the key for cfg_ok/cfg_unsafe above)
    cfg_no_key = joinpath(tmpdir, "cfg_no_key.toml")
    write(cfg_no_key, """
title = "unit-test"

[plots]
[planet]
[execution]
[physics]

[files]
""")
    @test_throws ErrorException AGNI.open_config(cfg_no_key)

    # make_logger: file logger must route all four log levels (existing tests above
    # only exercise INFO and WARN)
    logpath_all = joinpath(tmpdir, "agni_alllevels.log")
    logger_all = AGNI.make_logger(logpath_all; to_term=false)
    with_logger(logger_all) do
        @debug "dbg line"
        @info "info line"
        @warn "warn line"
        @error "err line"
    end
    logtxt_all = read(logpath_all, String)
    @test occursin("DEBUG", logtxt_all) && occursin("dbg line", logtxt_all)
    @test occursin("ERROR", logtxt_all) && occursin("err line", logtxt_all)

    # make_logger: terminal logger routes ERROR to stderr and other levels to
    # stdout. Levels are logged in this order (ERROR last) deliberately.
    out_path = joinpath(tmpdir, "term_stdout.txt")
    err_path = joinpath(tmpdir, "term_stderr.txt")
    open(out_path, "w") do out_io
        open(err_path, "w") do err_io
            redirect_stdout(out_io) do
                redirect_stderr(err_io) do
                    # constructed *inside* the redirect scope: make_logger captures
                    # the current stdout/stderr streams by value at call time
                    term_logger = AGNI.make_logger(""; to_term=true)
                    with_logger(term_logger) do
                        @debug "term debug"
                        @info "term info"
                        @warn "term warn"
                        @error "term error"
                    end
                end
            end
        end
    end
    out_txt = read(out_path, String)
    err_txt = read(err_path, String)
    @test occursin("term debug", out_txt)
    @test occursin("term info", out_txt)
    @test occursin("term warn", out_txt)
    @test occursin("term error", err_txt)

    # discrimination guard: ERROR must not also appear on stdout, and the levels
    # logged before it must not appear on stderr
    @test !occursin("term error", out_txt)
    @test !occursin("term debug", err_txt)
    @test !occursin("term info", err_txt)
    @test !occursin("term warn", err_txt)

    # make_logger: to_term=false & empty outpath falls back to a NullLogger, with a
    # diagnostic warning printed directly to stderr (not via the logging system)
    null_err_path = joinpath(tmpdir, "null_logger_stderr.txt")
    null_logger = open(null_err_path, "w") do io
        redirect_stderr(io) do
            AGNI.make_logger(""; to_term=false)
        end
    end
    @test null_logger isa Logging.NullLogger
    @test occursin("NullLogger", read(null_err_path, String))

    rm(tmpdir; force=true, recursive=true)
end

@testset "agni_config_parser" begin

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

# These cases use Test.collect_test_logs to inspect *which* validation branch fired,
# rather than only checking the boolean return value.
@testset "agni_from_config_extra" begin
    # gravity computed from `mass` alone (single-key success path), distinct from the
    # already-tested "mass AND gravity both given" overspecified-error branch above
    cfg = _base_cfg()
    delete!(cfg["planet"], "gravity")
    cfg["planet"]["mass"] = 5.972e24  # kg, Earth-like
    delete!(cfg["planet"], "albedo_s")
    logs, ok = Test.collect_test_logs() do
        AGNI.run_from_config(cfg)
    end
    @test ok == false
    errs = [l.message for l in logs if l.level == Logging.Error]
    @test any(occursin("albedo_s", m) for m in errs)

    # discrimination guard: if the mass-only branch were broken (e.g. treating the
    # default `gravity=0.0` as "provided"), the function would never reach albedo_s branch
    @test !any(occursin("provide `planet.mass` OR `planet.gravity`", m) for m in errs)

    # metallicities-driven composition + chemistry=true uses the dummy-VMR-dict path
    cfg = _base_cfg()
    delete!(cfg["composition"], "vmr_dict")
    cfg["composition"]["metallicities"] = Dict("C" => 1.0)
    cfg["physics"]["chemistry"] = true
    cfg["plots"]["extension"] = "invalidext"
    logs, ok = Test.collect_test_logs() do
        AGNI.run_from_config(cfg)
    end
    @test ok == false
    errs = [l.message for l in logs if l.level == Logging.Error]
    @test any(occursin("Plot extension", m) for m in errs)

    # discrimination guard: metallicities-without-chemistry must NOT be the failure reason
    @test !any(occursin("must enable FastChem", m) for m in errs)

    # neither `p_surf` nor `p_dict` provided (non-transparent)
    cfg = _base_cfg()
    delete!(cfg["composition"], "p_surf")
    @test !haskey(cfg["composition"], "p_dict")
    logs, ok = Test.collect_test_logs() do
        AGNI.run_from_config(cfg)
    end
    @test ok == false
    errs = [l.message for l in logs if l.level == Logging.Error]
    @test any(occursin("Must provide either", m) for m in errs)

    # RFM requested with only one of the two required wavenumber bounds
    cfg = _base_cfg()
    cfg["files"]["rfm_parfile"] = "dummy.par"
    cfg["execution"]["rfm_wn_min"] = 500.0
    logs, ok = Test.collect_test_logs() do
        AGNI.run_from_config(cfg)
    end
    @test ok == false
    errs = [l.message for l in logs if l.level == Logging.Error]
    @test any(occursin("RFM calculation enabled", m) for m in errs)

    # RFM requested with both wavenumber bounds present
    cfg = _base_cfg()
    cfg["files"]["rfm_parfile"] = "dummy.par"
    cfg["execution"]["rfm_wn_min"] = 500.0
    cfg["execution"]["rfm_wn_max"] = 1500.0
    cfg["plots"]["extension"] = "invalidext"
    logs, ok = Test.collect_test_logs() do
        AGNI.run_from_config(cfg)
    end
    @test ok == false
    errs = [l.message for l in logs if l.level == Logging.Error]
    @test any(occursin("Plot extension", m) for m in errs)
    @test !any(occursin("RFM calculation enabled", m) for m in errs)

    # grey opacities (`grey_lw`/`grey_sw`) success-assignment path
    cfg = _base_cfg()
    cfg["execution"]["grey_start"] = true
    cfg["physics"]["grey_lw"] = 1.0e-4
    cfg["physics"]["grey_sw"] = 0.0
    cfg["plots"]["extension"] = "invalidext"
    logs, ok = Test.collect_test_logs() do
        AGNI.run_from_config(cfg)
    end
    @test ok == false
    errs = [l.message for l in logs if l.level == Logging.Error]
    @test any(occursin("Plot extension", m) for m in errs)
    @test !any(occursin("Grey-gas calculation enabled", m) for m in errs)
end

# Greygas RT rejects scattering flags while accepting boundary opacity inputs.
@testset "agni_greygas_config" begin
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

# Minimal from-scratch greygas configuration dict, used to reach the parts of
# run_from_config that lie *after* atmosphere allocation
function _minimal_greygas_cfg(; out_dir::String=mktempdir())
    return Dict(
        "title" => "minimal greygas smoke config",
        "files" => Dict(
            "output_dir" => out_dir,
            "input_sf"   => "greygas",
            "input_star" => "",
        ),
        "planet" => Dict(
            "radius"          => 1.0e7,
            "surface_material"=> "greybody",
            "albedo_s"        => 0.2,
            "instellation"    => 1200.0,
            "s0_fact"         => 1.0,
            "albedo_b"        => 0.0,
            "zenith_angle"    => 0.0,
            "tmp_surf"        => 350.0,
            "gravity"         => 10.0,
        ),
        "composition" => Dict(
            "p_surf"      => 10.0,
            "p_top"       => 1.0e-6,
            "vmr_dict"    => Dict("H2O" => 1.0),
            "condensates" => String[],
        ),
        "execution" => Dict(
            "num_levels"     => 30,
            "solver"         => "newton",
            "solution_type"  => 1,
            "converge_atol"  => 1.0,
            "converge_rtol"  => 1.0e-2,
            "grey_start"     => false,
            "initial_state"  => Any["iso", "300"],
        ),
        "physics" => Dict(
            "continua"        => true,
            "rayleigh"        => false,
            "cloud"           => false,
            "rainout"         => false,
            "oceans"          => false,
            "overlap_method"  => "ee",
            "thermo_funct"    => false,
            "convection_crit" => "l",
            "convection"      => true,
            "conduction"      => false,
            "sensible_heat"   => false,
            "latent_heat"     => false,
            "chemistry"       => false,
            "real_gas"        => false,
            "grey_lw"         => 1.0e-4,
            "grey_sw"         => 0.0,
        ),
        "plots" => Dict(),
    )
end

@testset "agni_from_config_grey" begin
    # Invalid solver name
    cfg = _minimal_greygas_cfg()
    cfg["execution"]["solver"] = "notarealsolver"
    logs, ok = Test.collect_test_logs() do
        AGNI.run_from_config(cfg)
    end
    @test ok == false
    warns = [l.message for l in logs if l.level == Logging.Warn]
    @test any(occursin("Invalid solver", m) for m in warns)

    # Multi-column ("globe") simulation requested but `lons`/`lats` omitted
    cfg2 = _minimal_greygas_cfg()
    cfg2["planet"]["globe"] = Dict{String,Any}()
    logs2, ok2 = Test.collect_test_logs() do
        AGNI.run_from_config(cfg2)
    end
    @test ok2 == false
    errs2 = [l.message for l in logs2 if l.level == Logging.Error]
    @test any(occursin("globe.lons", m) for m in errs2)
end

@testset "agni_main" begin
    # nonexistent config path raises before any file I/O
    @test_throws ErrorException AGNI.main(cfg_path="/nonexistent/path/to/config.toml")

    # a config that parses successfully via open_config, but is missing required key
    tmpdir = mktempdir()
    out_dir = joinpath(tmpdir, "main_out")
    cfg_path = joinpath(tmpdir, "cfg_incomplete.toml")
    write(cfg_path, """
title = "incomplete config"

[plots]
[planet]
[execution]
    clean_output = false
    verbosity = 0
[physics]

[files]
output_dir = "$out_dir"
""")
    # main() calls setup_logging(), which mutates the process-global logger
    old_global_logger = global_logger()
    logs, ok = try
        Test.collect_test_logs() do
            AGNI.main(cfg_path=cfg_path)
        end
    finally
        global_logger(old_global_logger)
    end
    @test ok == false
    errs = [l.message for l in logs if l.level == Logging.Error]
    @test any(occursin("missing required key", m) for m in errs)
    rm(tmpdir; force=true, recursive=true)
end

