using Test
using AGNI

# Tests for src/energy/mie.jl: Mie scattering by homogeneous spheres, and integration over
# log-normal size distributions.
# Invariants and contract clauses exercised:
#   - efficiencies match an independent Mie implementation (miepython 3.3.0) across regimes
#   - small-particle (Rayleigh) limits: Q_sca ∝ x⁴, Q_abs ∝ x, g → 0
#   - energy conservation: non-absorbing spheres have Q_abs = 0
#   - extinction paradox: Q_ext → 2 for large, non-absorbing spheres
#   - log-normal quadrature reproduces analytic moments (mean volume, effective radius)
#   - mass coefficients scale as 1/ρ, and σ_g = 1 recovers the monodisperse result

const mie = AGNI.mie

@testset "mie" begin

    # -------------
    # Reference efficiencies from an independent implementation.
    # Computed with miepython 3.3.0 (efficiencies_mx, which uses m = n - ik; the sign of k is
    # flipped here to match the n + ik convention). Cases span dielectric, weakly absorbing,
    # metallic-like, and anomalous-dispersion (n<1, as in the SiO2 9 μm band) regimes, and
    # size parameters from 0.3 to 1000.
    # -------------
    @testset "efficiencies_match_independent_reference" begin
        # (n, k, x, Qext, Qsca, g)
        ref = [
            (1.55, 0.0,  5.213,  3.104995915080e+00, 3.104995915080e+00, 6.331044159947e-01),
            (1.5,  0.1,  1.0,    4.823704563470e-01, 2.087400183148e-01, 2.055966885409e-01),
            (1.5,  0.1,  10.0,   2.459790528444e+00, 1.235144209371e+00, 9.223496060998e-01),
            (1.5,  0.1,  100.0,  2.089821842804e+00, 1.132133971125e+00, 9.503916728872e-01),
            (1.33, 1e-8, 0.5,    6.773152104268e-03, 6.773139877345e-03, 4.546478178482e-02),
            (3.0,  4.0,  2.0,    2.970350954710e+00, 2.116087265144e+00, 4.289852262779e-01),
            (3.0,  4.0,  20.0,   2.316136364908e+00, 1.784228872057e+00, 6.371298199406e-01),
            (0.5,  1.5,  0.3,    2.588088598003e+00, 1.087196473010e-01, 4.469283026149e-03),
            (0.5,  1.5,  3.0,    3.392951756084e+00, 2.204559020100e+00, 6.164382232826e-01),
            (2.0,  1.0,  1000.0, 2.020999454279e+00, 1.259452936071e+00, 8.315570203015e-01),
        ]
        for (n, k, x, Qe_ref, Qs_ref, g_ref) in ref
            Qe, Qs, g = mie.mie_sphere(x, complex(n, k))
            @test isapprox(Qe, Qe_ref; rtol=1e-8)
            @test isapprox(Qs, Qs_ref; rtol=1e-8)
            @test isapprox(g,  g_ref;  atol=1e-8)
        end

        # Discrimination guard: neglecting absorption (k=0.1 → 0) at x=10 changes Qext by
        # far more than the tolerance, so the pins above are sensitive to k.
        Qe_noabs, _, _ = mie.mie_sphere(10.0, complex(1.5, 0.0))
        @test abs(Qe_noabs - 2.459790528444) > 1e-2
    end

    # -------------
    # Small particle limit. For x ≪ 1, with L = (m²-1)/(m²+2):
    #   Q_sca = (8/3) x⁴ |L|²,  Q_abs = 4x Im(L),  g → 0
    # (Bohren & Huffman 1983, Sec. 5.2). x values straddle the switch to the analytic
    # branch at x = 1e-3, so both code paths are exercised.
    # -------------
    @testset "rayleigh_limit_scaling" begin
        m = complex(1.5, 0.01)
        L = (m^2 - 1) / (m^2 + 2)
        for x in (5e-4, 2e-3, 1e-2)
            Qe, Qs, g = mie.mie_sphere(x, m)
            @test isapprox(Qs, 8/3 * x^4 * abs2(L); rtol=1e-3)
            @test isapprox(Qe - Qs, 4x * imag(L); rtol=1e-3)
            @test abs(g) < 1e-3
        end

        # Continuity across the analytic/series switch at X_RAYLEIGH
        xs = mie.X_RAYLEIGH
        Qa, _, _ = mie.mie_sphere(xs * (1 - 1e-9), m)
        Qb, _, _ = mie.mie_sphere(xs * (1 + 1e-9), m)
        @test isapprox(Qa, Qb; rtol=1e-5)

        # Discrimination guard: Q_sca must scale as x⁴, not x³ or x² (factor 10 in x gives
        # factor 1e4 in Q_sca, far from 1e3 or 1e2).
        _, Qs1, _ = mie.mie_sphere(1e-3 * 1.0001, m)
        _, Qs2, _ = mie.mie_sphere(1e-2, m)
        @test isapprox(Qs2 / Qs1, (1e-2 / (1e-3 * 1.0001))^4; rtol=1e-2)
    end

    # -------------
    # Energy conservation: a non-absorbing sphere (k = 0) scatters all extinguished light,
    # across the Rayleigh, resonance, and geometric regimes.
    # -------------
    @testset "nonabsorbing_sphere_conserves_energy" begin
        for x in (0.05, 1.0, 7.3, 60.0, 500.0)
            Qe, Qs, g = mie.mie_sphere(x, complex(1.33, 0.0))
            @test isapprox(Qe, Qs; rtol=1e-10)
            @test -1.0 <= g <= 1.0
        end
        # Absorbing sphere must have Q_abs > 0
        Qe, Qs, _ = mie.mie_sphere(7.3, complex(1.33, 0.01))
        @test Qe - Qs > 1e-3
    end

    # -------------
    # Extinction paradox: Q_ext → 2 as x → ∞ for a sphere, with diffraction contributing
    # half. At x = 5000 the residual oscillation is small.
    # -------------
    @testset "extinction_efficiency_tends_to_two" begin
        Qe, Qs, g = mie.mie_sphere(5000.0, complex(1.5, 0.0))
        @test isapprox(Qe, 2.0; atol=0.02)
        # large particles strongly forward scatter
        @test g > 0.7
        # Discrimination guard: the geometric cross-section alone would give Q_ext = 1
        @test abs(Qe - 1.0) > 0.5
    end

    # -------------
    # Gauss-Hermite nodes reproduce moments of the normal distribution N(0, 1/2):
    #   E[t²] = 1/2, E[t⁴] = 3/4. Weights sum to one.
    # -------------
    @testset "gauss_hermite_reproduces_normal_moments" begin
        for n in (8, 48)
            t, w = mie.gauss_hermite(n)
            @test isapprox(sum(w), 1.0; rtol=1e-12)
            @test isapprox(sum(w .* t.^2), 0.5; rtol=1e-10)
            @test isapprox(sum(w .* t.^4), 0.75; rtol=1e-10)
            @test abs(sum(w .* t)) < 1e-12
        end
        # single node degenerates to the mean
        t, w = mie.gauss_hermite(1)
        @test isapprox(t[1], 0.0; atol=1e-15)
        @test isapprox(w[1], 1.0; rtol=1e-15)
    end

    # -------------
    # Log-normal quadrature reproduces analytic moments:
    #   ⟨V⟩ = (4/3)π r_g³ exp(4.5 ln²σ_g),  r_eff = ⟨r³⟩/⟨r²⟩.
    # σ_g from 1 (monodisperse) to 2.5 (broad) covers typical cloud and haze widths.
    # -------------
    @testset "lognormal_quadrature_matches_analytic_moments" begin
        r_eff = 1.3e-6
        for σ in (1.0, 1.2, 1.65, 2.5)
            r, w = mie.lognormal_nodes(r_eff, σ)
            V_quad = sum(w .* (4/3 * π .* r.^3))
            @test isapprox(V_quad, mie.lognormal_mean_volume(r_eff, σ); rtol=1e-10)
            @test isapprox(sum(w .* r.^3) / sum(w .* r.^2), r_eff; rtol=1e-10)
            @test isapprox(sum(w), 1.0; rtol=1e-12)
        end
        # monodisperse limit returns a single node at r_eff
        r, w = mie.lognormal_nodes(r_eff, 1.0)
        @test length(r) == 1
        @test isapprox(r[1], r_eff; rtol=1e-15)

        # Discrimination guard: confusing r_eff with r_g would give a volume that differs
        # by exp(7.5 ln²σ) ≈ 2.3 at σ_g = 1.65
        V_wrong = 4/3 * π * r_eff^3 * exp(4.5 * log(1.65)^2)
        @test V_wrong / mie.lognormal_mean_volume(r_eff, 1.65) > 2.0
    end

    # -------------
    # Mass coefficients: k = ⟨σ⟩/(ρ⟨V⟩). For a monodisperse distribution this is
    #   k_ext = 3 Q_ext / (4 ρ r).
    # Coefficients must scale as 1/ρ, be non-negative, and give g in [-1, 1].
    # Densities chosen to match SiO2 glass (2201) and wustite (5970).
    # -------------
    @testset "mass_coefficients_scale_with_density" begin
        λ = [0.5e-6, 9.0e-6]
        m = [complex(1.46, 1e-6), complex(0.6, 1.5)]
        r = 1.0e-6

        ka, ks, g = mie.mass_coefficients(λ, m, r, 1.0, 2201.0)
        for i in 1:2
            Qe, Qs, gi = mie.mie_sphere(2π * r / λ[i], m[i])
            @test isapprox(ka[i] + ks[i], 3 * Qe / (4 * 2201.0 * r); rtol=1e-10)
            @test isapprox(ks[i], 3 * Qs / (4 * 2201.0 * r); rtol=1e-10)
            @test isapprox(g[i], gi; rtol=1e-10)
        end

        ka2, ks2, g2 = mie.mass_coefficients(λ, m, r, 1.0, 5970.0)
        @test all(isapprox.(ka2 .* 5970.0, ka .* 2201.0; rtol=1e-12))
        @test all(isapprox.(ks2 .* 5970.0, ks .* 2201.0; rtol=1e-12))
        @test all(isapprox.(g2, g; rtol=1e-12))
        # Discrimination guard: density ratio is far outside tolerance
        @test ka2[2] / ka[2] < 0.5

        # Broad distribution: physical bounds hold
        ka3, ks3, g3 = mie.mass_coefficients(λ, m, r, 2.0, 2201.0)
        @test all(ka3 .>= 0.0)
        @test all(ks3 .>= 0.0)
        @test all(-1.0 .<= g3 .<= 1.0)
        # visible light is almost purely scattered by weakly absorbing silica
        @test ks3[1] / (ka3[1] + ks3[1]) > 0.999
        # 9 μm band is strongly absorbing
        @test ka3[2] / (ka3[2] + ks3[2]) > 0.3
    end

end
