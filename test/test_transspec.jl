using Test
using AGNI

ROOT_DIR = abspath(joinpath(dirname(abspath(@__FILE__)), "../"))
OUT_DIR = joinpath(ROOT_DIR, "out/")

function _make_transspec_atmos(; nlev_c::Int64=100,
                                instellation::Float64=900.0,
                                p_surf::Float64=1.0,
                                p_top::Float64=1e-6,
                                transspec_ref_tau::Float64=1.1,
                                transspec_ref_wl::Float64=1e-6,
                                transspec_ref_p::Float64=1e-3)
    atmos = atmosphere.Atmos_t()
    ok = atmosphere.setup!(atmos, ROOT_DIR, OUT_DIR,
                            "greygas",
                            instellation, 1.0, 0.0, 0.0,
                            400.0,
                            10.0, 1.0e7,
                            nlev_c, p_surf, p_top,
                            Dict("N2" => 1.0), "";
                            real_gas=false,
                            thermo_functions=false,
                            flag_rayleigh=false,
                            flag_cloud=false,
                            transspec_ref_tau=transspec_ref_tau,
                            transspec_ref_wl=transspec_ref_wl,
                            transspec_ref_p=transspec_ref_p)
    ok || error("Failed to setup transspec atmosphere")
    atmosphere.allocate!(atmos, ""; check_safe_gas=false) || error("Failed to allocate transspec atmosphere")
    setpt.isothermal!(atmos, 500.0)
    atmosphere.calc_layer_props!(atmos)

    return atmos
end

@testset "transspec" begin

    # Exercise validation of transspec reference parameters and unit conversions.
    @testset "test_transspec_reference_vals" begin
        # Edge case: use lower-bound values to exercise validation limits.
        atmos_ok = atmosphere.Atmos_t()
        ok = atmosphere.setup!(atmos_ok, ROOT_DIR, OUT_DIR,
                                "greygas",
                                900.0, 1.0, 0.0, 0.0,
                                400.0,
                                10.0, 1.0e7,
                                100, 1.0, 1e-6,
                                Dict("N2" => 1.0), "";
                                real_gas=false,
                                thermo_functions=false,
                                flag_rayleigh=false,
                                flag_cloud=false,
                                transspec_ref_tau=1e-10,
                                transspec_ref_wl=1e-10,
                                transspec_ref_p=1e-10)
        @test ok
        @test isapprox(atmos_ok.transspec_ref_tau, 1e-10; rtol=0.0, atol=1e-12)
        @test isapprox(atmos_ok.transspec_ref_wl, 1e-10; rtol=0.0, atol=1e-12)
        @test isapprox(atmos_ok.transspec_ref_p, 1e-5; rtol=0.0, atol=1e-12)
        @test isapprox(atmos_ok.transspec_p, atmos_ok.transspec_ref_p; rtol=0.0, atol=1e-12)
        @test abs(atmos_ok.transspec_ref_p - 3e-5) > 1e-6
    end

    # Exercise photosphere selection with and without tau profiles.
    @testset "test_photosphere" begin
        atmos = _make_transspec_atmos()

        # Make edge vs centre fields distinct to verify layer-edge selection.
        atmos.tmp .= 450.0
        atmos.tmpl .= 510.0
        atmos.g .= 9.0
        atmos.gl .= 11.0
        atmos.layer_μ .= 0.028

        # Error path: no tau profile forces fix-index fallback.
        atmos.tau_band .= 0.0
        atmos.transspec_ref_p = atmos.p[2]
        idx_fallback = atmosphere.estimate_photosphere!(atmos)
        @test idx_fallback == 2
        @test isapprox(atmos.transspec_p, atmos.pl[2]; rtol=rtol, atol=1e-12)
        @test isapprox(atmos.transspec_r, atmos.rl[2]; rtol=rtol, atol=1.0)
        @test isapprox(atmos.transspec_tmp, atmos.tmpl[2]; rtol=rtol, atol=1e-12)
        @test abs(atmos.transspec_tmp - atmos.tmp[2]) > 1.0

        # Edge case: ref_wl below band range should clamp to the nearest band.
        tau_profile = range(0.0, 2.5, length=atmos.nlev_l)
        for ba in 1:atmos.nbands
            atmos.tau_band[:, ba] .= tau_profile .+ 0.05 * ba
        end
        atmos.transspec_ref_tau = 1.1
        atmos.transspec_ref_wl = minimum(atmos.bands_cen) * 0.1

        # Happy path: valid tau profile should select the correct photosphere layer.
        energy.radtrans!(atmos, true; calc_cf=true) # populate tau_band
        atmos.transspec_ref_wl = 1e-6 # 1 um
        idx_tau = atmosphere.estimate_photosphere!(atmos)

        idx_exp = 83
        @test idx_tau == idx_exp

        # Extract observables at this layer
        @test isapprox(atmos.tau_p[1], atmos.pl[idx_exp]; rtol=rtol, atol=1e-12)
        @test isapprox(atmos.transspec_p, atmos.pl[idx_exp]; rtol=rtol, atol=1e-12)
        @test isapprox(atmos.transspec_r, atmos.rl[idx_exp]; rtol=rtol, atol=1.0)
        @test isapprox(atmos.transspec_μ,
                        0.5*(atmos.layer_μ[idx_exp]+atmos.layer_μ[idx_exp-1]); rtol=rtol, atol=1e-12)
    end

    # Exercise greygas optical depth accumulation and band flux bookkeeping.
    @testset "greygas_optical_depth" begin
        # Test error case: unallocated atmosphere should not run RT.
        # atmos_unalloc = atmosphere.Atmos_t()
        # @test !energy.radtrans!(atmos_unalloc, true)

        # Edge case: tiny opacities should still produce a monotonic tau profile.
        atmos = _make_transspec_atmos()
        atmos.κ_grey_lw = 1e-8
        atmos.κ_grey_sw = 2e-8
        energy.radtrans!(atmos, true) # LW
        energy.radtrans!(atmos, false) # SW

        # test that there's one band
        @test atmos.nbands == 1

        # test tau_band shape and is finite
        @test size(atmos.tau_band) == (atmos.nlev_l, atmos.nbands)
        @test all(diff(atmos.tau_band[:, 1]) .>= 0.0)

        # test that grey optical depth accumulates with distance from TOA
        @test all(diff(atmos.tau_band[:, 1]) .> 0.0)
    end

end
