using Test
using AGNI

# Tests for src/energy/aerosol_optics.jl (runtime aerosol optical properties) and the
# condensate density table in src/phys/density.jl.
# Invariants and contract clauses exercised:
#   - refractive index files are parsed robustly (headers, CR line endings, extra columns,
#     unsorted/duplicated wavelengths) and unphysical rows are rejected
#   - interpolation reproduces tabulated values, is monotonic between nodes, and flags
#     extrapolation outside the table
#   - stellar-weighted band averaging returns the constant for constant properties, weights
#     g by scattering, falls back to uniform weighting where the star has no flux, and
#     reports the extrapolated fraction
#   - Mie band properties are non-negative, have |g| ≤ 1, and scale with 1/ρ
#   - every material with a density has a refractive index file (when data are present)

const aeropt = AGNI.aerosol_optics

# Write a temporary refractive index file, returning its path
function _write_nk(content::String)::String
    path = tempname() * ".txt"
    open(path, "w") do f
        write(f, content)
    end
    return path
end

@testset "aerosol_optics" begin

    # -------------
    # Parsing: '#' comments, un-commented header lines, CR-only line endings (as in some
    # WS15 files), a leading numeric-but-invalid header (as in gCMCRT files), 4-column files
    # with an index column (as in ExoHaze files), unsorted rows, and duplicated wavelengths.
    # -------------
    @testset "read_nk_parses_heterogeneous_formats" begin
        # CR line endings, un-commented header, unsorted
        p1 = _write_nk("KCl - temp_cond=740\rWav (microns)\tn\tk\r2.0\t1.5\t0.02\r1.0\t1.4\t0.01\r")
        λ, n, k = aeropt.read_nk(p1)
        @test length(λ) == 2
        @test isapprox(λ[1], 1.0e-6; rtol=1e-12)   # converted from micron, sorted
        @test isapprox(n[2], 1.5; rtol=1e-12)
        @test isapprox(k[1], 0.01; rtol=1e-12)

        # gCMCRT-style first line and 4-column (index, λ, n, k) rows use the last 3 columns
        p2 = _write_nk("2001 .False. # Number of lines\n# header\n1 0.5 1.6 0.001\n2 0.6 1.7 0.002\n")
        λ, n, k = aeropt.read_nk(p2)
        @test isapprox(λ, [0.5e-6, 0.6e-6]; rtol=1e-12)
        @test isapprox(n, [1.6, 1.7]; rtol=1e-12)

        # duplicated wavelength is averaged
        p3 = _write_nk("1.0 1.4 0.01\n1.0 1.6 0.03\n2.0 1.5 0.02\n")
        λ, n, k = aeropt.read_nk(p3)
        @test length(λ) == 2
        @test isapprox(n[1], 1.5; rtol=1e-12)
        @test isapprox(k[1], 0.02; rtol=1e-12)

        # Unphysical rows (n=0 as in the WS15 Fe2O3 file; negative k) are skipped
        p4 = _write_nk("1.0 1.4 0.01\n2.0 1.5 0.02\n3.0 0.0 0.5\n4.0 1.6 -1e-8\n")
        λ, n, k = aeropt.read_nk(p4)
        @test length(λ) == 2
        @test all(n .> 0.0)
        @test all(k .>= 0.0)

        # Error contract: missing file, or fewer than two valid rows
        @test_throws ErrorException aeropt.read_nk(tempname())
        p5 = _write_nk("# only a header\n1.0 1.4 0.01\n2.0 -1.0 0.0\n")
        @test_throws ErrorException aeropt.read_nk(p5)

        foreach(p->rm(p; force=true), (p1, p2, p3, p4, p5))
    end

    # -------------
    # Interpolation: exact at nodes, log-linear in k between positive nodes, linear where a
    # node has k=0, and clamped with an extrapolation flag outside the table.
    # -------------
    @testset "interp_nk_is_exact_at_nodes_and_flags_extrapolation" begin
        λt = [1e-6, 4e-6, 16e-6]
        nt = [1.2, 1.8, 1.4]
        kt = [1e-4, 1e-2, 0.0]

        m, ext = aeropt.interp_nk(λt, nt, kt, 4e-6)
        @test isapprox(real(m), 1.8; rtol=1e-12)
        @test isapprox(imag(m), 1e-2; rtol=1e-12)
        @test !ext

        # geometric midpoint of first interval: n is the arithmetic mean, k the geometric mean
        m, ext = aeropt.interp_nk(λt, nt, kt, 2e-6)
        @test isapprox(real(m), 1.5; rtol=1e-12)
        @test isapprox(imag(m), 1e-3; rtol=1e-12)
        @test !ext
        # Discrimination guard: linear interpolation of k would give 5.05e-3, not 1e-3
        @test abs(imag(m) - 5.05e-3) > 1e-3

        # interval with k=0 at one end: linear interpolation, stays non-negative
        m, _ = aeropt.interp_nk(λt, nt, kt, 8e-6)
        @test isapprox(imag(m), 0.5e-2; rtol=1e-12)

        # outside the table: edge values, flagged
        m, ext = aeropt.interp_nk(λt, nt, kt, 0.5e-6)
        @test ext
        @test isapprox(real(m), 1.2; rtol=1e-12)
        m, ext = aeropt.interp_nk(λt, nt, kt, 100e-6)
        @test ext
        @test isapprox(imag(m), 0.0; atol=1e-15)
    end

    # -------------
    # Band averaging with synthetic properties and a flat stellar spectrum.
    #   - constant properties average to the same constant
    #   - g is weighted by the scattering coefficient
    #   - bands beyond the stellar spectrum fall back to uniform weighting
    #   - the extrapolated fraction is the weighted fraction of flagged points
    # -------------
    @testset "band_average_weights_by_star_and_scattering" begin
        bands = [0.5e-6 1.0e-6;
                 1.0e-6 2.0e-6;
                 50e-6  100e-6]     # third band lies beyond the stellar spectrum
        star_wl = collect(range(100.0, 10000.0, length=2000))   # [nm]
        star_fl = ones(length(star_wl))
        star = aeropt.star_cumulative(star_wl, star_fl)

        # constant properties
        constprops(λ) = (fill(3.0, length(λ)), fill(7.0, length(λ)), fill(0.4, length(λ)),
                            fill(false, length(λ)))
        ka, ks, g, fx = aeropt.band_average(bands, star, Float64[], constprops)
        @test all(isapprox.(ka, 3.0; rtol=1e-12))
        @test all(isapprox.(ks, 7.0; rtol=1e-12))
        @test all(isapprox.(g, 0.4; rtol=1e-12))
        @test all(fx .< 1e-15)

        # g weighted by k_sca: half the band has k_sca=1, g=0.8; half has k_sca=3, g=0.0
        # (split at the band's centre in λ, with uniform flux per unit λ)
        function stepprops(λ)
            left = λ .< 0.75e-6
            return (zeros(length(λ)), ifelse.(left, 1.0, 3.0), ifelse.(left, 0.8, 0.0), left)
        end
        ka, ks, g, fx = aeropt.band_average(bands[1:1,:], star, [0.75e-6], stepprops)
        @test isapprox(ks[1], 2.0; rtol=0.05)
        @test isapprox(g[1], 0.8 * 1.0 / (1.0 + 3.0); rtol=0.05)   # = 0.2
        # Discrimination guard: an unweighted mean of g would give 0.4
        @test abs(g[1] - 0.4) > 0.1
        # half of the weight was flagged as extrapolated
        @test isapprox(fx[1], 0.5; rtol=0.05)

        # a property linear in λ is averaged over the band with flat weighting
        linprops(λ) = (λ ./ 1e-6, zeros(length(λ)), zeros(length(λ)), fill(false, length(λ)))
        ka, _, _, _ = aeropt.band_average(bands, star, Float64[], linprops)
        @test isapprox(ka[1], 0.75; rtol=1e-3)
        @test isapprox(ka[2], 1.5;  rtol=1e-3)
        @test isapprox(ka[3], 75.0; rtol=1e-3)   # uniform fallback beyond the star

        # Error contract: invalid band edges and malformed stellar spectrum
        @test_throws ErrorException aeropt.band_average([1e-6 1e-6], star, Float64[], constprops)
        @test_throws ErrorException aeropt.band_average(zeros(2, 3), star, Float64[], constprops)
        @test_throws ErrorException aeropt.star_cumulative([500.0], [1.0])
    end

    # -------------
    # Stellar weighting shifts band means towards the part of the band with more flux.
    # A star with all its flux in the short-λ half of a band must return the short-λ value.
    # -------------
    @testset "band_average_follows_stellar_flux" begin
        bands = reshape([1.0e-6, 2.0e-6], 1, 2)
        star_wl = collect(range(900.0, 2100.0, length=4000))
        star_fl = ifelse.(star_wl .< 1400.0, 1.0, 0.0)
        star = aeropt.star_cumulative(star_wl, star_fl)
        stepprops(λ) = (ifelse.(λ .< 1.45e-6, 1.0, 100.0), zeros(length(λ)), zeros(length(λ)),
                        fill(false, length(λ)))
        ka, _, _, _ = aeropt.band_average(bands, star, [1.4e-6, 1.45e-6], stepprops)
        @test isapprox(ka[1], 1.0; rtol=0.02)
        # Discrimination guard: uniform weighting would give ≈ 55
        @test ka[1] < 5.0
    end

    # -------------
    # Stellar spectra are often sparsely sampled in the far-infrared Rayleigh-Jeans tail
    # (e.g. sun.txt has points at 300, 400, 1000 μm). The cumulative integral uses power-law
    # segments, which are exact for f ∝ λ^-4:  ∫ λ^-4 dλ = (λ0^-3 - x^-3)/3.
    # -------------
    @testset "star_integral_exact_for_power_law_tail" begin
        wl_nm = [1.0e5, 1.0e6]                 # 100 and 1000 μm, two points only
        fl    = (wl_nm .* 1e-9) .^ -4
        star  = aeropt.star_cumulative(wl_nm, fl)
        λ0 = 1e-4
        for x in (2e-4, 3e-4, 7e-4, 1e-3)
            exact = (λ0^-3 - x^-3) / 3
            @test isapprox(aeropt.star_integral(star, x), exact; rtol=1e-10)
        end
        # clamped outside the tabulated range
        @test isapprox(aeropt.star_integral(star, 5e-5), 0.0; atol=1e-30)
        @test isapprox(aeropt.star_integral(star, 1.0), (λ0^-3 - 1e-3^-3)/3; rtol=1e-10)
        # Discrimination guard: a trapezoid over the two points overestimates by ~45x
        trap = 0.5 * (fl[1] + fl[2]) * (1e-3 - 1e-4)
        @test trap / aeropt.star_integral(star, 1e-3) > 10.0

        # f ∝ λ^-1 (p = -1) uses the logarithmic special case
        star1 = aeropt.star_cumulative(wl_nm, (wl_nm .* 1e-9) .^ -1)
        @test isapprox(aeropt.star_integral(star1, 5e-4), log(5.0); rtol=1e-10)
        # duplicated wavelengths in the input are removed
        star2 = aeropt.star_cumulative([1e5, 1e5, 1e6], [1.0, 1.0, 1.0])
        @test length(star2[1]) == 2
    end

    # -------------
    # Band grid includes the band edges, is sorted and unique, and subsamples very dense
    # refractive index tables.
    # -------------
    @testset "band_grid_includes_edges_and_limits_density" begin
        g = aeropt.band_grid(1e-6, 2e-6, collect(range(0.5e-6, 3e-6, length=100000)))
        @test isapprox(g[1], 1e-6; rtol=1e-12)
        @test isapprox(g[end], 2e-6; rtol=1e-12)
        @test issorted(g)
        @test allunique(g)
        @test length(g) <= aeropt.NLAM_BAND_MIN + 1 + aeropt.NLAM_BAND_NK
        @test length(g) > aeropt.NLAM_BAND_NK   # dense table points were included
    end

    # -------------
    # Density lookup: known materials return positive values; unknown materials raise an
    # error (a placeholder density would silently give wrong opacities).
    # -------------
    @testset "condensate_density_lookup" begin
        @test isapprox(AGNI.density.condensate_rho("SiO2_amorph"), 2201.0; rtol=1e-12)
        @test isapprox(AGNI.density.condensate_rho("FeO"), 5970.0; rtol=1e-12)
        @test all(AGNI.density.condensate_rho(m) > 1000.0 for m in AGNI.density.list_condensate_rho())
        # iron is the densest material in the table
        @test AGNI.density.condensate_rho("Fe") >= maximum(AGNI.density.condensate_rho.(AGNI.density.list_condensate_rho()))
        @test_throws ErrorException AGNI.density.condensate_rho("Unobtainium")
    end

    # -------------
    # Full Mie band averaging with real refractive index data (skipped if the data have not
    # been downloaded). SiO2 glass is transparent in the visible and absorbs strongly in the
    # 9 μm Si-O stretching band.
    # -------------
    @testset "mie_optics_for_real_material" begin
        if !isdir(AGNI.paths.get_dir("refractive")) || isempty(aeropt.list_materials())
            @test_skip "refractive index data not available"
        else
            # every material with a density has a refractive index file
            @test Set(aeropt.list_materials()) == Set(AGNI.density.list_condensate_rho())

            bands = [0.5e-6 0.6e-6;
                     8.5e-6 9.5e-6;
                     1.0e-3 2.0e-3]      # beyond the tabulated data (487 μm)
            star_wl = collect(range(100.0, 3.0e6, length=20000))
            star_fl = ones(length(star_wl))
            ka, ks, g, fx = aeropt.compute_mie_optics("SiO2_amorph", 1e-6, 1.65,
                                                        bands, star_wl, star_fl)
            @test all(ka .>= 0.0)
            @test all(ks .>= 0.0)
            @test all(-1.0 .<= g .<= 1.0)
            @test ks[1] / (ka[1] + ks[1]) > 0.999     # visible: scattering
            @test ka[2] / (ka[2] + ks[2]) > 0.3       # 9 μm: absorbing
            @test fx[1] < 1e-12
            @test isapprox(fx[3], 1.0; rtol=1e-12)    # fully extrapolated band

            # large-particle scale: k_ext ~ 3 Q_ext/(4 ρ r_eff) with Q_ext ~ 2-3
            kext_scale = 3 * 2.0 / (4 * 2201.0 * 1e-6)
            @test 0.5 < (ka[1] + ks[1]) / kext_scale < 2.0

            # Unknown material: error
            @test_throws ErrorException aeropt.compute_mie_optics("Unobtainium", 1e-6, 1.65,
                                                                    bands, star_wl, star_fl)
        end
    end

end
