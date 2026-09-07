# Tests for src/state/atmosphere.jl
# Covers:
#   - atmosphere.setup!() composition loading when mf_source==1 (VMR set from a CSV file),
#     including the pressure-domain extension guards and malformed-column handling.
#   - atmosphere.allocate!() rejecting a missing SOCRATES spectral file.
#   - atmosphere.allocate!() surface reflectance/emissivity loading: the greybody
#     path (Kirchhoff's law, albedo_s clamping) and the spectral-file path (r/e/w
#     header conventions, interpolation onto the model's spectral bands, and the
#     missing-file/bad-header/malformed-column error/gap paths).

using Test
using AGNI
using Logging

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)), "../"))
OUT_DIR  = joinpath(ROOT_DIR, "out/")

# Common cheap-to-run atmosphere parameters (greygas RT; no SOCRATES spectral read needed
# for setup! itself, since only allocate! touches the spectral file).
const _TMP_SURF  = 500.0
const _GRAVITY   = 10.0
const _RADIUS    = 6.37e6
const _NLEV      = 30
const _P_SURF    = 100.0  # bar -> p_boa = 1e7 Pa
const _P_TOP     = 1e-6   # bar -> p_toa = 0.1 Pa
const _THETA     = 60.0

# Cheap fixture used by the tests below that only require atmosphere.setup!()
function _setup_only(; condensates::Array{String,1}=String[], gravity::Float64=_GRAVITY)
    atmos = atmosphere.Atmos_t()
    ok = atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                            "greygas",
                            1000.0, 1.0, 0.0, _THETA,
                            _TMP_SURF,
                            gravity, _RADIUS,
                            _NLEV, _P_SURF, _P_TOP,
                            Dict("H2O" => 1.0), "";
                            real_gas=false,
                            thermo_functions=false,
                            flag_rayleigh=false,
                            flag_cloud=false,
                            condensates=condensates)
    ok || error("Failed to setup test atmosphere")
    return atmos
end

function _setup_with_vmr_file(mf_path::String)
    atmos = atmosphere.Atmos_t()
    ok = atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                            "greygas",
                            1000.0, 1.0, 0.0, _THETA,
                            _TMP_SURF,
                            _GRAVITY, _RADIUS,
                            _NLEV, _P_SURF, _P_TOP,
                            Dict{String,Float64}(), mf_path;
                            real_gas=false,
                            thermo_functions=false,
                            flag_rayleigh=false,
                            flag_cloud=false)
    return atmos, ok
end

# Cheap fixture for surface-reflectance/albedo tests: greygas RT gives
# nbands=1, bands_cen=[1e-6] m = 1000 nm (atmosphere.jl:1613-1617), and the
# surface-properties block (atmosphere.jl:2216-2315) runs unconditionally for
# both RT schemes, so allocate!(atmos, ""; check_safe_gas=false) reaches it
# without needing a real SOCRATES spectral file.
function _setup_with_surface(surface_material::String, albedo_s::Float64=0.0)
    atmos = atmosphere.Atmos_t()
    ok = atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                            "greygas",
                            1000.0, 1.0, 0.0, _THETA,
                            _TMP_SURF,
                            _GRAVITY, _RADIUS,
                            _NLEV, _P_SURF, _P_TOP,
                            Dict("H2O" => 1.0), "";
                            real_gas=false,
                            thermo_functions=false,
                            flag_rayleigh=false,
                            flag_cloud=false,
                            surface_material=surface_material,
                            albedo_s=albedo_s)
    ok || error("Failed to setup test atmosphere")
    return atmos
end

@testset "atmosphere" begin

    # -----------------------------------------------------------------
    # mf_source == 1 : composition read from a VMR CSV file
    # -----------------------------------------------------------------
    @testset "vmr_file_composition" begin

        tmpdir = mktempdir()

        # -------------------------------------------------------------
        # Happy path + pressure-domain extension edge case.
        #
        # File spans only 3 decades in pressure (1e2 to 1e5 Pa), which is
        # narrower than the model's full domain (p_toa=0.1 Pa, p_boa=1e7 Pa),
        # so both the low-pressure and high-pressure extension branches
        # (atmosphere.jl:1121-1130) must trigger.
        #
        # CO2 is defined as an exactly log10(p)-linear profile, and N2 is
        # its complement (N2 = 1 - CO2), so the per-level VMR normalisation
        # (sum to 1) is a no-op and does not obscure the interpolated values.
        # -------------------------------------------------------------
        vmr_path = joinpath(tmpdir, "vmr_profile.csv")
        p1, p2 = 1.0e2, 1.0e5   # Pa
        co2_1, co2_2 = 0.3, 0.9
        write(vmr_path, """
        #p,CO2,N2
        #[Pa],[-],[-]
        $(p1),$(co2_1),$(1-co2_1)
        $(p2),$(co2_2),$(1-co2_2)
        """)

        atmos, ok = _setup_with_vmr_file(vmr_path)
        @test ok
        @test "CO2" in atmos.gas_names
        @test "N2"  in atmos.gas_names

        # Sanity-check assumption behind the extension checks below: the model's
        # pressure grid must actually extend beyond the file's [p1, p2] range.
        @test atmos.p[1]   < p1
        @test atmos.p[end] > p2

        # Boundary extension: levels outside the file's pressure range take the
        # nearest edge value verbatim (not extrapolated or zeroed).
        @test isapprox(atmos.gas_vmr["CO2"][1],   co2_1;     atol=1e-8)
        @test isapprox(atmos.gas_vmr["N2"][1],    1-co2_1;   atol=1e-8)
        @test isapprox(atmos.gas_vmr["CO2"][end], co2_2;     atol=1e-8)
        @test isapprox(atmos.gas_vmr["N2"][end],  1-co2_2;   atol=1e-8)

        # Interior interpolation is linear in log10(pressure), not linear in
        # pressure. Pick a level strictly inside (p1, p2) and compare against
        # an independently-computed log-linear value.
        idx = findfirst(p -> (p > p1) && (p < p2), atmos.p)
        @test !isnothing(idx)
        logp = log10(atmos.p[idx])
        slope = (co2_2 - co2_1) / (log10(p2) - log10(p1))
        expect_correct = co2_1 + slope * (logp - log10(p1))
        @test isapprox(atmos.gas_vmr["CO2"][idx], expect_correct; atol=1e-6)

        # Discrimination guard: interpolating linearly in pressure (the wrong
        # basis, given pressure spans decades) gives a substantially different
        # value, so this pin cannot be satisfied by an incorrect implementation.
        frac_linear = (atmos.p[idx] - p1) / (p2 - p1)
        expect_wrong = co2_1 + (co2_2 - co2_1) * frac_linear
        @test abs(expect_correct - expect_wrong) > 1e-3

        # Composition must remain physically normalised at every level.
        for i in 1:atmos.nlev_c
            @test isapprox(atmos.gas_vmr["CO2"][i] + atmos.gas_vmr["N2"][i], 1.0; atol=1e-8)
        end

        atmosphere.deallocate!(atmos)

        # -------------------------------------------------------------
        # Duplicate gas column: a repeated header name must be skipped
        # (warned, not merged or overwritten), keeping the first occurrence.
        # -------------------------------------------------------------
        dup_path = joinpath(tmpdir, "vmr_dup.csv")
        write(dup_path, """
        #p,CO2,CO2,N2
        #[Pa],[-],[-],[-]
        $(p1),0.3,0.99,0.7
        $(p2),0.9,0.01,0.1
        """)

        atmos_dup, ok_dup = with_logger(MinLevelLogger(current_logger(), Logging.Error+1)) do
            _setup_with_vmr_file(dup_path)
        end
        @test ok_dup
        @test count(==("CO2"), atmos_dup.gas_names) == 1
        @test isapprox(atmos_dup.gas_vmr["CO2"][1],   0.3; atol=1e-8)
        @test isapprox(atmos_dup.gas_vmr["CO2"][end], 0.9; atol=1e-8)
        # discrimination guard: the discarded duplicate column's values must
        # NOT have been used instead
        @test !isapprox(atmos_dup.gas_vmr["CO2"][1], 0.99; atol=1e-8)
        atmosphere.deallocate!(atmos_dup)

        # -------------------------------------------------------------
        # Invalid gas name (non-alphanumeric column header): skipped with a
        # warning, while the remaining valid columns still load correctly.
        # -------------------------------------------------------------
        bad_path = joinpath(tmpdir, "vmr_badname.csv")
        write(bad_path, """
        #p,H2O,N2*,CO2
        #[Pa],[-],[-],[-]
        $(p1),0.7,0.5,0.3
        $(p2),0.4,0.5,0.6
        """)

        atmos_bad, ok_bad = with_logger(MinLevelLogger(current_logger(), Logging.Error+1)) do
            _setup_with_vmr_file(bad_path)
        end
        @test ok_bad
        @test "H2O" in atmos_bad.gas_names
        @test "CO2" in atmos_bad.gas_names
        @test !("N2*" in atmos_bad.gas_names)
        @test atmos_bad.gas_num == 2
        atmosphere.deallocate!(atmos_bad)

        # -------------------------------------------------------------
        # Missing file: setup! must fail gracefully (return false), not throw.
        # -------------------------------------------------------------
        missing_path = joinpath(tmpdir, "does_not_exist.csv")
        @test !isfile(missing_path)
        atmos_missing = atmosphere.Atmos_t()
        ok_missing = with_logger(MinLevelLogger(current_logger(), Logging.Error+1)) do
            (_, ok) = _setup_with_vmr_file(missing_path)
            ok
        end
        @test !ok_missing
    end

    # -----------------------------------------------------------------
    # allocate!() must reject a spectral file that does not exist on disk,
    # rather than erroring deep inside the SOCRATES wrapper.
    # -----------------------------------------------------------------
    @testset "spectral_file_not_found" begin
        atmos = atmosphere.Atmos_t()
        bogus_sf = joinpath(ROOT_DIR, "res", "spectral_files", "Dayspring", "16", "NoSuchFile.sf")
        @test !isfile(bogus_sf)

        ok = atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                                bogus_sf,
                                1000.0, 1.0, 0.0, _THETA,
                                _TMP_SURF,
                                _GRAVITY, _RADIUS,
                                _NLEV, _P_SURF, _P_TOP,
                                Dict("N2" => 1.0), "";
                                real_gas=false,
                                thermo_functions=false,
                                flag_rayleigh=false,
                                flag_cloud=false)
        # setup! itself only records the path; it does not check for existence
        @test ok
        @test atmos.rt_scheme == atmosphere.RT_SOCRATES
        @test atmos.spectral_file == abspath(bogus_sf)

        alloc_ok = with_logger(MinLevelLogger(current_logger(), Logging.Error+1)) do
            atmosphere.allocate!(atmos, "")
        end
        @test !alloc_ok
        @test !atmos.is_alloc
    end

    # -----------------------------------------------------------------
    # Surface reflectance/emissivity, greybody path (atmosphere.jl:2220-2224).
    # Kirchhoff's law sets emissivity = 1 - albedo spectrally; albedo_s is
    # clamped to [0,1] by setup!() (atmosphere.jl:795) before allocate!() ever
    # sees it.
    # -----------------------------------------------------------------
    @testset "surface_reflectance_greybody" begin
        # mid-range albedo: check the Kirchhoff relation e = 1-r holds exactly,
        # not just that the arrays get populated with something
        atmos = _setup_with_surface("greybody", 0.3)
        @test atmosphere.allocate!(atmos, ""; check_safe_gas=false)
        @test length(atmos.surf_r_arr) == atmos.nbands == 1
        @test isapprox(atmos.surf_r_arr[1], 0.3; atol=1e-12)
        @test isapprox(atmos.surf_e_arr[1], 0.7; atol=1e-12)
        # discrimination guard: emissivity must not simply mirror albedo
        @test !isapprox(atmos.surf_e_arr[1], atmos.surf_r_arr[1]; atol=1e-3)
        atmosphere.deallocate!(atmos)

        # boundary: albedo_s above 1 must be clamped to 1, giving a perfectly
        # reflecting / non-emitting surface
        atmos_hi = _setup_with_surface("greybody", 1.5)
        @test isapprox(atmos_hi.albedo_s, 1.0; atol=1e-12)
        @test atmosphere.allocate!(atmos_hi, ""; check_safe_gas=false)
        @test isapprox(atmos_hi.surf_r_arr[1], 1.0; atol=1e-12)
        @test isapprox(atmos_hi.surf_e_arr[1], 0.0; atol=1e-12)
        atmosphere.deallocate!(atmos_hi)

        # boundary: albedo_s below 0 must be clamped to 0
        atmos_lo = _setup_with_surface("greybody", -0.5)
        @test isapprox(atmos_lo.albedo_s, 0.0; atol=1e-12)
        @test atmosphere.allocate!(atmos_lo, ""; check_safe_gas=false)
        @test isapprox(atmos_lo.surf_r_arr[1], 0.0; atol=1e-12)
        @test isapprox(atmos_lo.surf_e_arr[1], 1.0; atol=1e-12)
        atmosphere.deallocate!(atmos_lo)
    end

    # -----------------------------------------------------------------
    # Surface reflectance, spectral file, header token "r" (spherical
    # reflectance, AKA Bond albedo; atmosphere.jl:2281-2286). Values are
    # interpolated linearly in wavelength [m]; the greygas band centre is
    # bands_cen=1e-6 m = 1000 nm (atmosphere.jl:1617).
    # -----------------------------------------------------------------
    @testset "surface_reflectance_spectral_r" begin
        tmpdir = mktempdir()

        # Two points straddling the 1000 nm band centre, so the interpolated
        # value is a genuine linear interpolation, not one input value verbatim.
        r_path = joinpath(tmpdir, "surf_r.dat")
        write(r_path, """
        r
        500.0 0.2
        2000.0 0.8
        """)

        atmos = _setup_with_surface(r_path)
        @test atmosphere.allocate!(atmos, ""; check_safe_gas=false)

        # reference calculation: linear interpolation in wavelength [nm],
        # independent of the implementation under test
        expect_r = 0.2 + (0.8 - 0.2) * (1000.0 - 500.0) / (2000.0 - 500.0)
        @test isapprox(atmos.surf_r_arr[1], expect_r; atol=1e-6)
        @test isapprox(atmos.surf_e_arr[1], 1.0 - expect_r; atol=1e-6)

        # discrimination guard: nearest-endpoint (wrong interpolation basis)
        # would give a substantially different value than true linear interp
        @test abs(expect_r - 0.2) > 1e-2
        @test abs(expect_r - 0.8) > 1e-2

        # the scalar albedo_s is overwritten with the median of the spectral
        # array (atmosphere.jl:2314); with nbands==1 that's just the one value
        @test isapprox(atmos.albedo_s, expect_r; atol=1e-6)
        atmosphere.deallocate!(atmos)

        # edge case: single-row file. Both ends of the wavelength domain are
        # extrapolated with a constant copy of that one value (atmosphere.jl:
        # 2261-2267), so the interpolated result must equal it exactly.
        r1_path = joinpath(tmpdir, "surf_r_single.dat")
        write(r1_path, """
        r
        1000.0 0.42
        """)
        atmos1 = _setup_with_surface(r1_path)
        @test atmosphere.allocate!(atmos1, ""; check_safe_gas=false)
        @test isapprox(atmos1.surf_r_arr[1], 0.42; atol=1e-10)
        @test isapprox(atmos1.surf_e_arr[1], 0.58; atol=1e-10)
        atmosphere.deallocate!(atmos1)
    end

    # -----------------------------------------------------------------
    # Surface reflectance, spectral file, header token "e" (hemispherical
    # emissivity; atmosphere.jl:2288-2293). Values chosen asymmetrically so
    # that swapping the r/e roles (a plausible implementation bug) is
    # detectable rather than accidentally passing.
    # -----------------------------------------------------------------
    @testset "surface_reflectance_spectral_e" begin
        tmpdir = mktempdir()
        e_path = joinpath(tmpdir, "surf_e.dat")
        write(e_path, """
        e
        500.0 0.2
        2000.0 0.5
        """)

        atmos = _setup_with_surface(e_path)
        @test atmosphere.allocate!(atmos, ""; check_safe_gas=false)

        expect_e = 0.2 + (0.5 - 0.2) * (1000.0 - 500.0) / (2000.0 - 500.0)
        @test isapprox(atmos.surf_e_arr[1], expect_e; atol=1e-6)
        @test isapprox(atmos.surf_r_arr[1], 1.0 - expect_e; atol=1e-6)

        # discrimination guard: if "e" were mistakenly treated as "r" (roles
        # swapped), surf_r_arr would equal expect_e instead of 1-expect_e;
        # these differ substantially here by construction
        @test abs((1.0 - expect_e) - expect_e) > 0.1

        atmosphere.deallocate!(atmos)
    end

    # -----------------------------------------------------------------
    # Surface reflectance, spectral file, header token "w" (single-scattering
    # albedo; atmosphere.jl:2296-2303, Hapke2012 diffusive-reflectance form).
    # The expected r/e values are recomputed here from the same published
    # formula, independently of the source under test, as the discriminating
    # reference calculation (per repo testing standards: a wrong-formula
    # implementation cannot pass this by construction).
    # -----------------------------------------------------------------
    @testset "surface_reflectance_spectral_w" begin
        tmpdir = mktempdir()
        w_path = joinpath(tmpdir, "surf_w.dat")
        # constant single-scattering albedo: interpolation is trivial here, so
        # the only thing under test is the w -> (r,e) conversion formula
        write(w_path, """
        w
        500.0 0.5
        2000.0 0.5
        """)

        atmos = _setup_with_surface(w_path)
        @test atmosphere.allocate!(atmos, ""; check_safe_gas=false)

        w = 0.5
        γ = sqrt(1 - w)
        r0 = (1 - γ) / (1 + γ)
        expect_r = r0 * (1 - γ/(3+3γ))
        expect_e = 1.0 - expect_r

        @test isapprox(atmos.surf_r_arr[1], expect_r; atol=1e-10)
        @test isapprox(atmos.surf_e_arr[1], expect_e; atol=1e-10)

        # discrimination guard: naively treating w as if it were r directly
        # gives a substantially different value from the correct result
        @test abs(expect_r - w) > 0.2

        atmosphere.deallocate!(atmos)
    end

    # -----------------------------------------------------------------
    # Surface reflectance: missing spectral-albedo file must fail gracefully
    # (return false), not throw (atmosphere.jl:2248-2252).
    # -----------------------------------------------------------------
    @testset "surface_reflectance_missing_file" begin
        tmpdir = mktempdir()
        missing_path = joinpath(tmpdir, "does_not_exist.dat")
        @test !isfile(missing_path)

        atmos = _setup_with_surface(missing_path)
        logs, alloc_ok = Test.collect_test_logs() do
            atmosphere.allocate!(atmos, ""; check_safe_gas=false)
        end
        @test !alloc_ok
        @test !atmos.is_alloc
        @test any(occursin("Could not find surface albedo file", l.message)
                        for l in logs if l.level == Logging.Error)
    end

    # -----------------------------------------------------------------
    # Surface reflectance: header token other than r/e/w must fail gracefully
    # (return false), not throw (atmosphere.jl:2305-2311).
    # -----------------------------------------------------------------
    @testset "surface_reflectance_bad_header" begin
        tmpdir = mktempdir()
        bad_path = joinpath(tmpdir, "surf_bad_header.dat")
        write(bad_path, """
        q
        500.0 0.2
        2000.0 0.8
        """)

        atmos = _setup_with_surface(bad_path)
        logs, alloc_ok = Test.collect_test_logs() do
            atmosphere.allocate!(atmos, ""; check_safe_gas=false)
        end
        @test !alloc_ok
        @test !atmos.is_alloc
        @test any(occursin("Unexpected format for surface data", l.message)
                        for l in logs if l.level == Logging.Error)
    end

    # -----------------------------------------------------------------
    # Gap (documented, not fixed here): a surface-data file with only one
    # column (no "value" column) is not validated before use - readdlm
    # succeeds, but atmosphere.jl:2259 (_srf_data[:,2]) then throws
    # BoundsError instead of returning false gracefully like the two cases
    # above. This pins the *current* behaviour so a future change is
    # deliberate rather than a silent regression.
    # -----------------------------------------------------------------
    @testset "surface_reflectance_malformed_single_column" begin
        tmpdir = mktempdir()
        onecol_path = joinpath(tmpdir, "surf_onecol.dat")
        write(onecol_path, """
        r
        500.0
        2000.0
        """)

        atmos = _setup_with_surface(onecol_path)
        @test_throws BoundsError atmosphere.allocate!(atmos, ""; check_safe_gas=false)
    end

    # -----------------------------------------------------------------
    # _check_range: pure boundary-check helper, no Atmos_t required at all.
    # -----------------------------------------------------------------
    @testset "check_range" begin
        # value within [min,max]: success, no logging
        @test atmosphere._check_range("x", 15.0; min=10.0, max=20.0)

        # below min only
        @test atmosphere._check_range("x", 15.0; min=10.0) == true
        logs, ok = with_logger(MinLevelLogger(current_logger(), Logging.Error+1)) do
            Test.collect_test_logs() do
                atmosphere._check_range("too_small", 5.0; min=10.0)
            end
        end
        @test ok == false

        # above max only
        logs, ok = Test.collect_test_logs() do
            atmosphere._check_range("too_large", 25.0; max=20.0)
        end
        @test ok == false
        @test any(occursin("too_large", l.message) for l in logs if l.level == Logging.Error)

        # both min and max given, value out of range on the low side: must report the
        # combined "$name is out of range" form, not the single-bound "too small" form
        logs, ok = Test.collect_test_logs() do
            atmosphere._check_range("both_bounds", 5.0; min=10.0, max=20.0)
        end
        @test ok == false
        errs = [l.message for l in logs if l.level == Logging.Error]
        @test any(occursin("out of range", m) for m in errs)
        @test !any(occursin("too small", m) for m in errs)
    end

    # -----------------------------------------------------------------
    # generate_pgrid!: only needs p_boa/p_toa/nlev_c/nlev_l set on the struct, so
    # this is tested without calling setup!() or allocate!() at all.
    # -----------------------------------------------------------------
    @testset "generate_pgrid_low_pressure_ratio" begin
        atmos = atmosphere.Atmos_t()
        atmos.p_toa = 1.0
        atmos.p_boa = 1.0  # ratio = 1.0, below PRESSURE_RATIO_MIN (1.0001)
        atmos.nlev_c = 10
        atmos.nlev_l = 11
        logs, _ = Test.collect_test_logs() do
            atmosphere.generate_pgrid!(atmos)
        end

        # test that the low-pressure-ratio warning was logged
        warns = [l.message for l in logs if l.level == Logging.Warn]
        @test any(occursin("pressure ratio", m) for m in warns)

        # discrimination guard: p_boa must be corrected to exceed the original
        @test atmos.p_boa > 1.0
        @test isapprox(atmos.p_boa, atmosphere.PRESSURE_RATIO_MIN * atmos.p_toa; rtol=1e-12)
        @test length(atmos.p) == 10
        @test length(atmos.pl) == 11
        @test issorted(atmos.pl; rev=false)  # increasing index -> increasing pressure
    end

    # -----------------------------------------------------------------
    # set_deep_heating!: argument-validation branches only (the physics of the
    # heating profile itself is covered by test_deep_heating.jl, which is excluded
    # from the fast tier). Only needs setup!() (reads atmos.p_toa/p_boa).
    # -----------------------------------------------------------------
    @testset "set_deep_heating_validation" begin
        atmos = _setup_only()

        logs, ok = Test.collect_test_logs() do
            atmosphere.set_deep_heating!(atmos, 1e5, 1.0, 0.1, 0.0, "mass", "clamp", "badmode")
        end
        @test ok == false
        @test any(occursin("Invalid deep heating power mode", l.message)
                        for l in logs if l.level == Logging.Error)

        logs, ok = Test.collect_test_logs() do
            atmosphere.set_deep_heating!(atmos, 1e5, 1.0, 0.1, 0.0, "badnorm", "clamp", "rel")
        end
        @test ok == false
        @test any(occursin("Invalid deep heating normalisation", l.message)
                        for l in logs if l.level == Logging.Error)

        logs, ok = Test.collect_test_logs() do
            atmosphere.set_deep_heating!(atmos, 1e5, 1.0, 0.1, 0.0, "mass", "baddomain", "rel")
        end
        @test ok == false
        @test any(occursin("Invalid deep heating domain treatment", l.message)
                        for l in logs if l.level == Logging.Error)

        # power_mode="off" is a valid no-op success path, distinct from the failure
        # modes above and from the verbose success case below
        @test atmosphere.set_deep_heating!(atmos, 1e5, 1.0, 0.1, 0.0, "mass", "clamp", "off")
        @test atmos.deepheat_power_mode == "off"

        # verbose success path logs a description of the configured profile
        logs, ok = Test.collect_test_logs() do
            atmosphere.set_deep_heating!(atmos, 1e5, 1.0, 0.1, 0.0, "mass", "clamp", "rel";
                                            verbose=true)
        end
        @test ok == true
        infos = [l.message for l in logs if l.level == Logging.Info]
        @test any(occursin("power_mode=rel", m) for m in infos)
        @test any(occursin("norm_method=mass", m) for m in infos)
    end

    # -----------------------------------------------------------------
    # calc_layer_props! and calc_profile_radius!
    # -----------------------------------------------------------------
    @testset "calc_layer_props_without_allocate" begin
        atmos = _setup_only()

        # baseline: a physically-reasonable surface gravity produces a fully-bound
        # atmosphere with positive density/scale-height everywhere
        @test atmosphere.calc_layer_props!(atmos)
        @test all(atmos.layer_isbound)
        @test all(atmos.layer_ρ .> 0.0)
        @test all(atmos.layer_Hp .> 0.0)

        # radius decreases from surface (i=nlev_c) to TOA (i=1)
        @test issorted(atmos.r; rev=true)

        # gravity below the internal ming floor (1e-4 m/s^2), so the function
        # must report failure and flag the affected layers
        atmos_lowg = _setup_only(; gravity=1e-7)
        ok = with_logger(MinLevelLogger(current_logger(), Logging.Error+1)) do
            atmosphere.calc_layer_props!(atmos_lowg)
        end
        @test ok == false
        @test any(.!atmos_lowg.layer_isbound)

        # discrimination guard: this must not be trivially "all layers unbound"
        @test all(atmos_lowg.g .<= atmos.hydrograv_ming)
    end

    # -----------------------------------------------------------------
    # set_cloud!: only needs condensates/gas_sat from setup!(), no allocate!().
    # -----------------------------------------------------------------
    @testset "set_cloud" begin
        atmos = _setup_only(; condensates=["H2O"])

        # from_yield=true (default), but no condensation has occurred yet (uniform
        # hot initial profile from setup!()), so there is nothing to form clouds from
        @test atmosphere.set_cloud!(atmos) == false
        @test all(atmos.cloud_arr_l .== 0.0)

        # from_yield=false uses the saturation mask (gas_sat) instead of from yield
        fill!(atmos.gas_sat["H2O"], false)
        atmos.gas_sat["H2O"][5]  = true
        atmos.gas_sat["H2O"][10] = true
        any_cloud = atmosphere.set_cloud!(atmos; from_yield=false)
        @test any_cloud == true
        @test isapprox(atmos.cloud_arr_l[5],  atmos.cloud_val_l; rtol=1e-10)
        @test isapprox(atmos.cloud_arr_l[10], atmos.cloud_val_l; rtol=1e-10)

        # discrimination guard: an un-saturated layer is set above the numerical floor
        @test atmos.cloud_arr_l[1] > 0.0
        @test atmos.cloud_arr_l[1] < atmos.cloud_val_l

        # negative cloud particle size
        atmos.cloud_val_r = -1.0
        logs, _ = Test.collect_test_logs() do
            atmosphere.set_cloud!(atmos; from_yield=false)
        end
        @test any(occursin("Negative cloud particle size", l.message)
                        for l in logs if l.level == Logging.Warn)
    end

    # -----------------------------------------------------------------
    # set_aerosol! and set_aerosols!
    # -----------------------------------------------------------------
    @testset "set_aerosol_and_aerosols" begin
        atmos = _setup_only(; condensates=["H2O"])
        atmos.aerosol_arr_l["testaer"] = zeros(Float64, atmos.nlev_c)
        atmos.aerosol_arr_r["testaer"] = zeros(Float64, atmos.nlev_c)

        # 1D-array mmr profile branch
        mmr_profile = fill(0.5, atmos.nlev_c)
        @test atmosphere.set_aerosol!(atmos, "testaer", mmr_profile)
        @test all(isapprox.(atmos.aerosol_arr_l["testaer"], 0.5; rtol=1e-10))

        # populate layer_σ, needed by calc_cond_mmr
        atmosphere.calc_layer_props!(atmos)

        # populate with a small nonzero condensation yield
        atmos.cond_yield["H2O"][5] = 0.5
        atmos.aerosol_names = ["testaer", "orphanaer"]
        atmos.aerosol_setby["testaer"] = "H2O"

        # "orphanaer" deliberately has no aerosol_setby entry and no array allocated
        any_aerosol = atmosphere.set_aerosols!(atmos)
        @test any_aerosol == true
        expect5 = atmosphere.calc_cond_mmr(atmos, "H2O", 5)
        @test expect5 > 0.0  # sanity: the discrimination below is non-trivial
        @test isapprox(atmos.aerosol_arr_l["testaer"][5], expect5; rtol=1e-10)
        @test !haskey(atmos.aerosol_arr_l, "orphanaer")
    end

    # -----------------------------------------------------------------
    # _iphot_from_prs! to test pressure index search
    # -----------------------------------------------------------------
    @testset "iphot_from_prs_guard_and_search" begin
        atmos = _setup_only()

        # not allocated
        @test !atmos.is_alloc
        logs, idx = Test.collect_test_logs() do
            atmosphere._iphot_from_prs!(atmos, 1e3)
        end
        @test idx == 1
        @test any(occursin("not been allocated", l.message)
                        for l in logs if l.level == Logging.Warn)

        # is_alloc hand-set true, then check we get the correct layer
        atmos.is_alloc = true
        idx2 = atmosphere._iphot_from_prs!(atmos, atmos.pl[7])
        @test idx2 == 7

        # discrimination guard: a reference pressure roughly half-way down the grid
        idx3 = atmosphere._iphot_from_prs!(atmos, sqrt(atmos.pl[1]*atmos.pl[end]))
        @test 1 < idx3 < atmos.nlev_l
    end

    # -----------------------------------------------------------------
    # list_available_aerosols
    # -----------------------------------------------------------------
    @testset "list_available_aerosols_disabled" begin
        atmos = _setup_only()
        @test !atmos.control.l_aerosol
        names = atmosphere.list_available_aerosols(atmos)
        @test names == String[]
    end
end
