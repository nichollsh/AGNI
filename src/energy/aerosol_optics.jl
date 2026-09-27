# This file is part of AGNI. License is Apache-2.0: https://apache.org/licenses/LICENSE-2.0

"""
**Contains module for calculating aerosol optical properties at runtime**

Aerosol and cloud particle optical properties are calculated from the complex refractive
index (n, k) of the particle material using Mie theory (see the `mie` module), for a
log-normal size distribution. The resultant spectral mass absorption coefficient, mass
scattering coefficient, and asymmetry parameter are averaged over each band of a SOCRATES
spectral file, weighted by the stellar spectrum ('thin' averaging).
"""
module aerosol_optics

    import ..paths
    import ..density
    import ..mie

    using LoggingExtras

    # Minimum number of log-spaced wavelength points per band
    const NLAM_BAND_MIN::Int = 16

    # Maximum number of refractive index data points to include per band
    const NLAM_BAND_NK::Int = 128

    # Warn when more than this fraction of the weight in a band is extrapolated
    const EXTRAP_WARN::Float64 = 0.01

    """
    **Path to the refractive index file for a given material.**

    Arguments:
    - `material::String`    material name (e.g. "SiO2_amorph")

    Returns:
    - `path::String`        path to file
    """
    function nk_path(material::String)::String
        return joinpath(paths.get_dir("refractive"), material*".txt")
    end

    """
    **List materials which have both a refractive index file and a known density.**

    Returns:
    - `materials::Vector{String}`   sorted list of supported materials
    """
    function list_materials()::Vector{String}
        return sort([m for m in density.list_condensate_rho() if isfile(nk_path(m))])
    end

    """
    **Read refractive index data from a text file.**

    Accepts files with three numeric columns (wavelength [micron], n, k), or more columns in
    which case the last three are used. Data are sorted by wavelength, and duplicated
    wavelengths are averaged.

    Arguments:
    - `path::String`            path to txt file

    Returns on success:
    - `λ::Vector{Float64}`      wavelength [m]
    - `n::Vector{Float64}`      real part of the refractive index
    - `k::Vector{Float64}`      imaginary part of the refractive index (≥ 0)

    Returns on failure:
    - `false`                   failure to read file
    """
    function read_nk(path::String)::Union{NTuple{3,Vector{Float64}},Bool}
        if !isfile(path)
            @warn("Refractive index file not found: '$path'")
            return false
        end

        text = replace(read(path, String), "\r\n"=>"\n", "\r"=>"\n")
        rows = Tuple{Float64,Float64,Float64}[]
        for line in split(text, '\n')
            line = strip(split(line, '#')[1])
            isempty(line) && continue
            vals = tryparse.(Float64, split(line))
            if (length(vals) >= 3) && !any(isnothing, vals)
                row = (vals[end-2], vals[end-1], vals[end])
                if all(isfinite, row) && (row[1] > 0.0) && (row[2] > 0.0) && (row[3] >= 0.0)
                    push!(rows, row)
                end
            end
        end

        if length(rows) < 2
            @warn("Refractive index file '$path' contains fewer than two valid data rows")
            return false
        end
        if length(rows) < 2
            @warn("Refractive index file '$path' contains fewer than two wavelengths")
            return false
        end

        # Sort by wavelength
        sort!(rows, by=r->r[1])

        # Average duplicated wavelengths
        λ = Float64[]
        n = Float64[]
        k = Float64[]
        i = 1
        while i <= length(rows)
            j = i
            while (j < length(rows)) && (rows[j+1][1] == rows[i][1])
                j += 1
            end
            push!(λ, rows[i][1] * 1e-6) # convert from microns to metres
            push!(n, sum(r[2] for r in rows[i:j]) / (j-i+1))
            push!(k, sum(r[3] for r in rows[i:j]) / (j-i+1))
            i = j + 1
        end

        if length(λ) < 2
            @warn("Refractive index file '$path' contains fewer than two wavelengths")
            return false
        end

        return (λ, n, k)
    end

    """
    **Interpolate refractive index to a given wavelength.**

    The real part is interpolated linearly in log(λ). The imaginary part is interpolated
    linearly in log(k) against log(λ) where both neighbouring values are positive, and
    linearly otherwise. Outside of the tabulated range, the edge values are used.

    Arguments:
    - `λ_tab, n_tab, k_tab`     tabulated data, as returned by `read_nk`
    - `λ::Float64`              wavelength at which to evaluate [m]

    Returns:
    - `m::ComplexF64`           complex refractive index n + ik
    - `extrap::Bool`            whether the wavelength is outside the tabulated range
    """
    function interp_nk(λ_tab::Vector{Float64}, n_tab::Vector{Float64}, k_tab::Vector{Float64},
                        λ::Float64)::Tuple{ComplexF64,Bool}
        if λ <= λ_tab[1]
            return (complex(n_tab[1], k_tab[1]), λ < λ_tab[1])
        elseif λ >= λ_tab[end]
            return (complex(n_tab[end], k_tab[end]), λ > λ_tab[end])
        end

        i = searchsortedlast(λ_tab, λ)
        f = (log(λ) - log(λ_tab[i])) / (log(λ_tab[i+1]) - log(λ_tab[i]))
        n = n_tab[i] + f * (n_tab[i+1] - n_tab[i])
        if (k_tab[i] > 0.0) && (k_tab[i+1] > 0.0)
            k = exp(log(k_tab[i]) + f * (log(k_tab[i+1]) - log(k_tab[i])))
        else
            k = k_tab[i] + f * (k_tab[i+1] - k_tab[i])
        end
        return (complex(n, k), false)
    end

    """
    **Wavelength grid within a single correlated-k band.**

    Union of log-spaced points spanning the band and (a subset of) the tabulated refractive
    index wavelengths inside the band.

    Arguments:
    - `λ_lo::Float64`           short wavelength edge of band [m]
    - `λ_hi::Float64`           long wavelength edge of band [m]
    - `λ_tab::Vector{Float64}`  tabulated refractive index wavelengths [m]

    Returns:
    - `λ::Vector{Float64}`      sorted wavelength grid including the band edges [m]
    """
    function band_grid(λ_lo::Float64, λ_hi::Float64, λ_tab::Vector{Float64})::Vector{Float64}

        # Log-spaced grid spanning the band, with NLAM_BAND_MIN points
        grid = collect(exp.(range(log(λ_lo), log(λ_hi), length=NLAM_BAND_MIN+1)))

        # get wavelengths within the band (not including the edges)
        inside = λ_tab[(λ_tab .> λ_lo) .& (λ_tab .< λ_hi)]

        # if there are too many points, downsample to NLAM_BAND_NK
        if length(inside) > NLAM_BAND_NK
            idx = round.(Int, range(1, length(inside), length=NLAM_BAND_NK))
            inside = inside[unique(idx)]
        end

        # return the total set of points
        return unique(sort(vcat(grid, inside)))
    end

    """
    **Integration weights of the stellar spectrum over a wavelength grid.**

    Each grid point is assigned the integral of the stellar flux over its cell, with cell
    boundaries at the midpoints between grid points. The stellar spectrum is integrated at
    its native resolution, and is taken to be zero outside of its tabulated range.

    Arguments:
    - `λ::Vector{Float64}`      sorted wavelength grid [m]
    - `star::NTuple{3,...}`     stellar spectrum, as returned by `star_cumulative`

    Returns:
    - `w::Vector{Float64}`      weights (arbitrary units)
    - `dλ::Vector{Float64}`     cell widths [m]
    """
    function cell_weights(λ::Vector{Float64},
                            star::NTuple{3,Vector{Float64}})::NTuple{2,Vector{Float64}}
        N = length(λ)
        edges = zeros(Float64, N+1)
        edges[1]   = λ[1]
        edges[end] = λ[end]
        for j in 2:N
            # The edge between λ[j-1] and λ[j] is the midpoint in log-space
            edges[j] = 0.5 * (λ[j-1] + λ[j])
        end

        # Integrate the stellar spectrum over each cell to get the weights
        C = [star_integral(star, e) for e in edges]

        # Return weights and widths
        return (max.(diff(C), 0.0), diff(edges))
    end

    """
    Integral of the flux across one segment of the stellar spectrum, from λ0 to x.
    Uses a power law between the end points where both fluxes are positive (exact for the
    λ^-4 Rayleigh-Jeans tail, which is often sparsely sampled), and linear otherwise.
    """
    function _segment_integral(λ0::Float64, f0::Float64, λ1::Float64, f1::Float64,
                                x::Float64)::Float64
        if (f0 > 0.0) && (f1 > 0.0)
            p = log(f1/f0) / log(λ1/λ0)
            if abs(p + 1.0) < 1e-8
                return f0 * λ0 * log(x/λ0)
            end
            return f0 * λ0 * ((x/λ0)^(p + 1.0) - 1.0) / (p + 1.0)
        end
        fx = f0 + (x - λ0) / (λ1 - λ0) * (f1 - f0)
        return 0.5 * (f0 + fx) * (x - λ0)
    end

    """
    **Cumulative integral of the stellar spectrum up to a given wavelength.**

    Arguments:
    - `star::NTuple{3,...}`     stellar spectrum, as returned by `star_cumulative`
    - `x::Float64`              wavelength [m]

    Returns:
    - `C::Float64`              integral of the flux from the shortest wavelength up to `x`
    """
    function star_integral(star::NTuple{3,Vector{Float64}}, x::Float64)::Float64
        λ, f, C = star
        if x <= λ[1]
            return 0.0
        elseif x >= λ[end]
            return C[end]
        end
        i = searchsortedlast(λ, x)
        return C[i] + _segment_integral(λ[i], f[i], λ[i+1], f[i+1], x)
    end

    """
    **Prepare a stellar spectrum for use as a weighting function.**

    Sorts the spectrum, converts wavelengths to metres, and tabulates its cumulative
    integral using power-law segments within each band (see `_segment_integral`).

    Arguments:
    - `star_wl::Vector{Float64}`    wavelengths [nm]
    - `star_fl::Vector{Float64}`    spectral flux (any units, per unit wavelength)

    Returns on success:
    - `λ::Vector{Float64}`          sorted, unique wavelengths [m]
    - `f::Vector{Float64}`          non-negative flux at each wavelength
    - `C::Vector{Float64}`          cumulative integral of flux at each wavelength

    Returns on failure:
    - `false`                       failure to prepare the stellar spectrum
    """
    function star_cumulative(star_wl::Vector{Float64},
                                star_fl::Vector{Float64})::Union{Bool,NTuple{3,Vector{Float64}}}

        # Sort the spectrum by wavelength and remove duplicates
        perm = sortperm(star_wl)

        # Convert to metres, and remove non-positive fluxes and duplicate wavelengths
        λ = star_wl[perm] .* 1e-9
        f = max.(star_fl[perm], 0.0)
        keep = vcat(true, diff(λ) .> 0.0)
        λ = λ[keep]
        f = f[keep]

        # Set cumulative integral to zero at first wavelength
        C = zeros(Float64, length(λ))

        # Loop over each segment of the spectrum,
        # integrating from the previous wavelength to the current one.
        # This provides weights for the band-averaging of aerosol optical properties.
        for i in 2:length(λ)
            C[i] = C[i-1] + _segment_integral(λ[i-1], f[i-1], λ[i], f[i], λ[i])
        end

        # Return the sorted wavelengths, fluxes, and cumulative integral
        return (λ, f, C)
    end

    """
    **Band-averaged aerosol optical properties.**

    Performs 'thin' averaging as in SOCRATES `scatter_average`, weighted by the stellar
    spectrum:
        k̄ = ∫ k w dλ / ∫ w dλ
        ḡ = ∫ g k_sca w dλ / ∫ k_sca w dλ.
    Where the stellar spectrum has no flux within a band, uniform weighting is used.

    Arguments:
    - `bands::Matrix{Float64}`      band edges [m], size (nbands, 2)
    - `star::NTuple{3,...}`         stellar spectrum, as returned by `star_cumulative`
    - `λ_tab::Vector{Float64}`      wavelengths at which to include extra grid points [m]
    - `props::Function`             function mapping a wavelength vector [m] to a tuple
                                    `(k_abs, k_sca, g, extrap)` of equal-length vectors

    Returns:
    - `k_abs::Vector{Float64}`      band-mean mass absorption coefficient [m2 kg-1]
    - `k_sca::Vector{Float64}`      band-mean mass scattering coefficient [m2 kg-1]
    - `g::Vector{Float64}`          band-mean asymmetry parameter
    - `f_ext::Vector{Float64}`      fraction of weight in each band which was extrapolated
    """
    function band_average(bands::Matrix{Float64}, star::NTuple{3,Vector{Float64}},
                            λ_tab::Vector{Float64},
                            props::Function)::NTuple{4,Vector{Float64}}

        # Number of bands (`bands` has shape (nb, 2))
        nb = size(bands, 1)

        # Set values in bands to zero
        k_abs = zeros(Float64, nb)
        k_sca = zeros(Float64, nb)
        g     = zeros(Float64, nb)
        f_ext = zeros(Float64, nb)

        # Loop through bands
        for b in 1:nb

            # Get band edges
            λ_lo, λ_hi = minmax(bands[b,1], bands[b,2])

            # Get wavelength grid *within* this band
            λ = band_grid(λ_lo, λ_hi, λ_tab)

            # Determine the weights and widths of each cell in the band,
            # using the stellar spectrum, so that the properties are weighted by where
            # the stellar spectrum has increased flux.
            w, dλ = cell_weights(λ, star)
            if !(sum(w) > 0.0)
                w = dλ
            end
            w = w ./ sum(w) # normalize weights to sum to 1

            # Compute the optical properties for this band
            ka, ks, gg, ext = props(λ)

            k_abs[b] = sum(w .* ka)
            k_sca[b] = sum(w .* ks)
            g[b]     = k_sca[b] > 0.0 ? sum(w .* ks .* gg) / k_sca[b] : 0.0
            f_ext[b] = sum(w[ext])
        end

        return (k_abs, k_sca, g, f_ext)
    end

    """
    **Calculate band-averaged optical properties of an aerosol from Mie theory.**

    Arguments:
    - `material::String`            material name, selecting refractive index and density
    - `r_eff::Float64`              effective radius [m]
    - `σ_g::Float64`                geometric standard deviation of log-normal distribution
    - `bands::Matrix{Float64}`      band edges [m], size (nbands, 2)
    - `star_wl::Vector{Float64}`    stellar spectrum wavelengths [nm]
    - `star_fl::Vector{Float64}`    stellar spectral flux (any units, per unit wavelength)

    Returns on success:
    - `k_abs::Vector{Float64}`      band-mean mass absorption coefficient [m2 kg-1]
    - `k_sca::Vector{Float64}`      band-mean mass scattering coefficient [m2 kg-1]
    - `g::Vector{Float64}`          band-mean asymmetry parameter
    - `f_ext::Vector{Float64}`      fraction of weight in each band which was extrapolated

    Returns on failure:
    - `false`                       failure to read refractive index file
    """
    function compute_mie_optics(material::String, r_eff::Float64, σ_g::Float64,
                                bands::Matrix{Float64},
                                star_wl::Vector{Float64},
                                star_fl::Vector{Float64})::Union{NTuple{4,Vector{Float64}},Bool}

        # Get the density of this material
        ρ = density.condensate_rho(material)

        # Read the refractive index data for this material
        read_nk_return = read_nk(nk_path(material))
        if read_nk_return === false
            @warn("Failed to read refractive index for material '$material'")
            return false
        else
            λ_tab, n_tab, k_tab = read_nk_return
        end

        # Prepare a stellar spectrum for use as a weighting function
        star = star_cumulative(star_wl, star_fl)

        function props(λ::Vector{Float64})
            mx  = Vector{ComplexF64}(undef, length(λ))
            ext = Vector{Bool}(undef, length(λ))
            for i in eachindex(λ)
                mx[i], ext[i] = interp_nk(λ_tab, n_tab, k_tab, λ[i])
            end
            ka, ks, gg = mie.mass_coefficients(λ, mx, r_eff, σ_g, ρ)
            return (ka, ks, gg, ext)
        end

        k_abs, k_sca, g, f_ext = band_average(bands, star, λ_tab, props)

        if maximum(f_ext) > EXTRAP_WARN
            nb_ext = count(f_ext .> EXTRAP_WARN)
            @warn "Refractive index of '$material' extrapolated in $nb_ext bands " *
                    "(tabulated $(round(λ_tab[1]*1e6, sigdigits=3))-" *
                    "$(round(λ_tab[end]*1e6, sigdigits=3)) μm)"
        end

        return (k_abs, k_sca, g, f_ext)
    end

end # end module
