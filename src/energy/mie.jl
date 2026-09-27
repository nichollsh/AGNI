# This file is part of AGNI. License is Apache-2.0: https://apache.org/licenses/LICENSE-2.0

"""
**Contains module for Mie scattering calculations**

Scattering and absorption by homogeneous spheres, following the BHMIE algorithm of
Bohren & Huffman (1983), with the series truncation criterion of Wiscombe (1980).
Polydisperse properties are integrated over a log-normal size distribution using
Gauss-Hermite quadrature in ln(r).

https://www.astro.princeton.edu/~draine/code/bhmie.f

References:
- Bohren, C. F. & Huffman, D. R. (1983), Absorption and Scattering of Light by Small
    Particles, Wiley. https://doi.org/10.1002/9783527618156
- Wiscombe, W. J. (1980), Improved Mie scattering algorithms, Applied Optics 19, 1505.
    https://doi.org/10.1364/AO.19.001505
"""
module mie

    import LinearAlgebra: SymTridiagonal, eigen

    # Below this size parameter, the Rayleigh (small particle) limit is used
    const X_RAYLEIGH::Float64 = 1.0e-3

    # Above this size parameter, the size parameter is clamped (efficiencies are
    #    asymptotically independent of x in the geometric optics limit)
    const X_MAX::Float64 = 2.0e4

    # Default number of quadrature nodes for size distribution integration
    const N_QUAD_DEFAULT::Int = 48

    # Quadrature nodes with normalised weights below this value are skipped
    const W_QUAD_MIN::Float64 = 1.0e-30

    """
    **Mie efficiencies of a homogeneous sphere in the Rayleigh limit.**

    Valid for size parameter x << 1 and |m|x << 1.

    Arguments:
    - `x::Float64`          size parameter, 2πr/λ
    - `m::ComplexF64`       complex refractive index n + ik, with k ≥ 0 for absorption

    Returns:
    - `Qext::Float64`       extinction efficiency
    - `Qsca::Float64`       scattering efficiency
    - `g::Float64`          asymmetry parameter
    """
    function mie_rayleigh(x::Float64, m::ComplexF64)::Tuple{Float64,Float64,Float64}
        L = (m^2 - 1) / (m^2 + 2)
        Qsca = 8.0/3.0 * x^4 * abs2(L)
        Qabs = 4.0 * x * imag(L)
        return (Qabs + Qsca, Qsca, 0.0)
    end

    """
    **Mie efficiencies of a homogeneous sphere.**

    Implementation of the BHMIE algorithm (Bohren & Huffman 1983, Appendix A), with the
    number of terms set by Wiscombe (1980). Falls back to the Rayleigh limit for very small
    size parameters, and clamps very large size parameters to `X_MAX`.

    Arguments:
    - `x::Float64`          size parameter, 2πr/λ
    - `m::ComplexF64`       complex refractive index n + ik, with k ≥ 0 for absorption

    Returns:
    - `Qext::Float64`       extinction efficiency
    - `Qsca::Float64`       scattering efficiency
    - `g::Float64`          asymmetry parameter
    """
    function mie_sphere(x::Float64, m::ComplexF64)::Tuple{Float64,Float64,Float64}

        # Small particle limit
        if x < X_RAYLEIGH
            return mie_rayleigh(x, m)
        end

        # Large particle limit
        x = min(x, X_MAX)

        y = m * x
        nstop = floor(Int, x + 4.05 * cbrt(x) + 2.0)
        nmx   = floor(Int, max(nstop, abs(y))) + 15

        # Logarithmic derivative D_n(mx), by downward recurrence
        D = zeros(ComplexF64, nmx)
        for n in nmx:-1:2
            D[n-1] = n/y - 1.0/(D[n] + n/y)
        end

        # Riccati-Bessel functions, by upward recurrence
        psi0 = cos(x)
        psi1 = sin(x)
        chi0 = -sin(x)
        chi1 = cos(x)
        xi1  = complex(psi1, -chi1)

        qsca = 0.0   # efficiency factor for scattering
        qext = 0.0   # efficiency factor for extinction
        gsum = 0.0
        an1  = zero(ComplexF64)
        bn1  = zero(ComplexF64)

        for n in 1:nstop
            fn  = Float64(n)
            psi = (2fn - 1.0)/x * psi1 - psi0
            chi = (2fn - 1.0)/x * chi1 - chi0
            xi  = complex(psi, -chi)

            ta = D[n]/m + fn/x
            tb = m*D[n] + fn/x
            an = (ta*psi - psi1) / (ta*xi - xi1)
            bn = (tb*psi - psi1) / (tb*xi - xi1)

            qsca += (2fn + 1.0) * (abs2(an) + abs2(bn))
            qext += (2fn + 1.0) * real(an + bn)
            gsum += (2fn + 1.0)/(fn*(fn + 1.0)) * real(an*conj(bn))
            if n > 1
                gsum += (fn - 1.0)*(fn + 1.0)/fn * real(an1*conj(an) + bn1*conj(bn))
            end

            an1  = an
            bn1  = bn
            psi0 = psi1
            psi1 = psi
            chi0 = chi1
            chi1 = chi
            xi1  = complex(psi1, -chi1)
        end

        # g = (4/(x² Qsca)) Σ[...] = 2 Σ[...] / (unnormalised Qsca sum)
        g    = 2.0 * gsum / qsca
        qsca = 2.0 / x^2 * qsca
        qext = 2.0 / x^2 * qext  

        return (qext, qsca, g)
    end

    """
    **Gauss-Hermite quadrature nodes and weights.**

    Computed with the Golub-Welsch algorithm. Weights are normalised such that
    `sum(w .* f.(t))` approximates the expectation of f(t) for t ~ N(0, 1/2).

    Arguments:
    - `n::Int`                  number of nodes

    Returns:
    - `t::Vector{Float64}`      nodes
    - `w::Vector{Float64}`      normalised weights (sum to unity)
    """
    function gauss_hermite(n::Int)::Tuple{Vector{Float64},Vector{Float64}}
        n = max(n, 1)
        if n == 1
            return ([0.0], [1.0])
        end
        J = SymTridiagonal(zeros(n), [sqrt(i/2.0) for i in 1:n-1])
        F = eigen(J)
        t = F.values
        w = F.vectors[1,:].^2
        return (t, w ./ sum(w))
    end

    """
    **Geometric mean radius of a log-normal distribution from its effective radius.**

    For a log-normal number distribution, r_eff = r_g exp(2.5 ln²σ_g).

    Arguments:
    - `r_eff::Float64`      effective (area-weighted mean) radius [m]
    - `σ_g::Float64`        geometric standard deviation (≥ 1)

    Returns:
    - `r_g::Float64`        geometric mean radius [m]
    """
    function geometric_radius(r_eff::Float64, σ_g::Float64)::Float64
        return r_eff * exp(-2.5 * log(σ_g)^2)
    end

    """
    **Analytic mean particle volume of a log-normal distribution.**

    Arguments:
    - `r_eff::Float64`      effective radius [m]
    - `σ_g::Float64`        geometric standard deviation (≥ 1)

    Returns:
    - `V::Float64`          number-weighted mean volume [m3]
    """
    function lognormal_mean_volume(r_eff::Float64, σ_g::Float64)::Float64
        r_g = geometric_radius(r_eff, σ_g)
        return 4.0/3.0 * π * r_g^3 * exp(4.5 * log(σ_g)^2)
    end

    """
    **Quadrature nodes for a log-normal size distribution.**

    Nodes with negligible weight are removed, and the remaining weights renormalised.
    A monodisperse distribution (σ_g = 1) returns a single node at r_eff.

    Arguments:
    - `r_eff::Float64`      effective radius [m]
    - `σ_g::Float64`        geometric standard deviation (≥ 1)

    Optional arguments:
    - `n::Int`              number of Gauss-Hermite nodes

    Returns:
    - `r::Vector{Float64}`  radii [m]
    - `w::Vector{Float64}`  number weights (sum to unity)
    """
    function lognormal_nodes(r_eff::Float64, σ_g::Float64;
                                n::Int=N_QUAD_DEFAULT)::Tuple{Vector{Float64},Vector{Float64}}

        # Convert σ_g to the standard deviation of ln(r), s = ln(σ_g)
        s = log(σ_g)
        if s < 1e-8
            return ([r_eff], [1.0])
        end

        # Calculate nodes and weights
        t, w = gauss_hermite(n)

        # Skip nodes with negligible weight
        mask = w .> W_QUAD_MIN
        t = t[mask]

        # Renormalise weights and convert to radii using r = r_g exp(sqrt(2) s t)
        w = w[mask] ./ sum(w[mask])
        r = geometric_radius(r_eff, σ_g) .* exp.(sqrt(2.0) * s .* t)
        return (r, w)
    end

    """
    **Optical properties averaged over a log-normal size distribution.**

    Calculated for each wavelength and refractive index, using Gauss-Hermite quadrature.

    Arguments:
    - `λ::Float64`          wavelength [m]
    - `m::ComplexF64`       complex refractive index n + ik
    - `r_eff::Float64`      effective radius [m]
    - `σ_g::Float64`        geometric standard deviation (≥ 1)

    Optional arguments:
    - `n::Int`              number of quadrature nodes

    Returns:
    - `σ_ext::Float64`      mean extinction cross-section per particle [m2]
    - `σ_sca::Float64`      mean scattering cross-section per particle [m2]
    - `g::Float64`          scattering-weighted asymmetry parameter
    - `V::Float64`          mean particle volume [m3]
    """
    function polydisperse(λ::Float64, m::ComplexF64, r_eff::Float64, σ_g::Float64;
                            n::Int=N_QUAD_DEFAULT)::NTuple{4,Float64}
        λ = max(λ, 1e-12)
        r, w = lognormal_nodes(r_eff, σ_g; n=n) # nodes and weights
        return _polydisperse(λ, m, r, w)
    end

    function _polydisperse(λ::Float64, m::ComplexF64,
                            r::Vector{Float64}, w::Vector{Float64})::NTuple{4,Float64}
        σ_ext = 0.0 # mean extinction cross-section per particle
        σ_sca = 0.0 # mean scattering cross-section per particle
        gsca  = 0.0 # mean scattering cross-section weighted by asymmetry parameter
        V     = 0.0
        for i in eachindex(r)
            area = π * r[i]^2
            Qe, Qs, g = mie_sphere(2π * r[i] / λ, m)
            σ_ext += w[i] * area * Qe
            σ_sca += w[i] * area * Qs
            gsca  += w[i] * area * Qs * g
            V     += w[i] * 4.0/3.0 * π * r[i]^3
        end
        g = σ_sca > 0.0 ? gsca / σ_sca : 0.0
        return (σ_ext, σ_sca, g, V)
    end

    """
    **Mass absorption and scattering coefficients over a log-normal size distribution.**

    Coefficients are per unit mass of condensate, as required by SOCRATES:
    k = ⟨σ⟩ / (ρ ⟨V⟩).

    Arguments:
    - `λ::Vector{Float64}`      wavelengths [m]
    - `m::Vector{ComplexF64}`   complex refractive index at each wavelength
    - `r_eff::Float64`          effective radius [m]
    - `σ_g::Float64`            geometric standard deviation (≥ 1)
    - `ρ::Float64`              bulk density of the particle material [kg m-3]

    Optional arguments:
    - `n::Int`                  number of quadrature nodes

    Returns:
    - `k_abs::Vector{Float64}`  mass absorption coefficient [m2 kg-1]
    - `k_sca::Vector{Float64}`  mass scattering coefficient [m2 kg-1]
    - `g::Vector{Float64}`      asymmetry parameter
    """
    function mass_coefficients(λ::Vector{Float64}, m::Vector{ComplexF64},
                                r_eff::Float64, σ_g::Float64, ρ::Float64;
                                n::Int=N_QUAD_DEFAULT
                                )::Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}

        r, w = lognormal_nodes(r_eff, σ_g; n=n)

        nλ = length(λ)
        k_abs = zeros(Float64, nλ)
        k_sca = zeros(Float64, nλ)
        g     = zeros(Float64, nλ)
        for i in 1:nλ
            σe, σs, gi, V = _polydisperse(λ[i], m[i], r, w)
            k_abs[i] = max(σe - σs, 0.0) / (ρ * V)
            k_sca[i] = σs / (ρ * V)
            g[i]     = gi
        end
        return (k_abs, k_sca, g)
    end

end # end module
