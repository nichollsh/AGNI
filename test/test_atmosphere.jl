# Tests for src/state/atmosphere.jl
# Covers:
#   - atmosphere.setup!() composition loading when mf_source==1 (VMR set from a CSV file),
#     including the pressure-domain extension guards and malformed-column handling.
#   - atmosphere.allocate!() rejecting a missing SOCRATES spectral file.

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
end
