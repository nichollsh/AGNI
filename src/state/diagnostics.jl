module diagnostics

    import ..phys
    import ..atmosphere

    """
    **Get pressure at top and bottom of convective zone**

    Arguments:
        - `atmos::Atmos_t`          the atmosphere struct instance to be used.

    Returns:
        - `p_top::Float64`          pressure [Pa] at top of convective zone
        - `p_bot::Float64`          pressure [Pa] at bottom of convective zone
    """
    function estimate_convective_zone(atmos::atmosphere.Atmos_t)::Tuple{Float64,Float64}

        # Defaults to zero, if there's no convection
        p_top::Float64 = 0.0
        p_bot::Float64 = 0.0

        # Loop from top-down to find p_top
        for i in 1:atmos.nlev_l
            if atmos.mask_c[i]
                p_top = atmos.pl[i]
                break
            end
        end

        # Loop from bottom-up to find p_bot
        for i in range(start=atmos.nlev_l, stop=1, step=-1)
            if atmos.mask_c[i]
                p_bot = atmos.pl[i]
                break
            end
        end

        # Return top, bot
        return (p_top, p_bot)
    end


    """
    **Estimate photosphere.**

    Estimates the location of the photosphere by finding the median of the contribution
    function in each band (0.2 um to 150 um), and then finding the pressure level at which
    these median values are maximised.

    Arguments:
        - `atmos::Atmos_t`          the atmosphere struct instance to be used.

    Returns:
        - `p_ref::Float64`          pressure level of photosphere [Pa]
    """
    function estimate_photosphere!(atmos::atmosphere.Atmos_t)::Float64

        # Params
        wl_min::Float64  = 0.2 * 1e-6 # 200 nm
        wl_max::Float64  = 150 * 1e-6 # 150 um
        p_min::Float64   = 10.0 # 1e-4 bar

        # tracking
        cff_max::Float64 = 0.0
        cff_try::Float64 = 0.0
        atmos.transspec_p = p_min

        # get band indices
        wl_imin = findmin(abs.(atmos.bands_cen .- wl_min))[2]
        wl_imax = findmin(abs.(atmos.bands_cen .- wl_max))[2]

        # reversed?
        if wl_imin > wl_imax
            wl_imin, wl_imax = wl_imax, wl_imin
        end

        # loop over levels
        for i in 1:atmos.nlev_c
            if atmos.p[i] < p_min
                continue
            end

            # maximum contfunc in this band
            cff_try = Statistics.median(atmos.contfunc_band[i,wl_imin:wl_imax])

            # is this more than the existing maximum?
            if cff_try > cff_max
                cff_max = cff_try
                atmos.transspec_p = atmos.p[i]
            end
        end

        return atmos.transspec_p
    end

    """
    **Calculate observed radius and bulk density.**

    This is done at the layer probed in transmission, which is set to a fixed pressure.

    Arguments:
        - `atmos::Atmos_t`          the atmosphere struct instance to be used.

    Returns:
        - `transspec_rho::Float64`  the bulk density observed in transmission
    """
    function calc_observed_rho!(atmos::atmosphere.Atmos_t)::Float64

        # transspec_r::Float64            # planet radius probed in transmission [m]
        # transspec_m::Float64            # mass [kg] of atmosphere + interior
        # transspec_rho::Float64          # bulk density [kg m-3] implied by r and m

        # Store reference pressure in atmos struct
        # estimate_photosphere!(atmos)

        # get the observed height
        idx::Int64 = findmin(abs.(atmos.p .- atmos.transspec_p))[2]
        atmos.transspec_r    = atmos.r[idx]
        atmos.transspec_μ    = atmos.layer_μ[idx]
        atmos.transspec_tmp  = atmos.tmp[idx]
        atmos.transspec_grav = atmos.g[idx]

        # get mass of whole atmosphere, assuming hydrostatic
        atmos.transspec_m = atmos.p_boa * 4 * pi * atmos.rp^2 / atmos.grav_surf

        # add mass of the interior component
        atmos.transspec_m += atmos.interior_mass

        # get density of all enclosed by observed layer
        # store this in the atmosphere struct
        atmos.transspec_rho = 3.0 * atmos.transspec_m / (4.0 * pi * atmos.transspec_r^3)

        # also return the value
        return atmos.transspec_rho
    end

    """
    **Estimate a diagnostic Rayleigh number in each layer.**

    Assuming that the Rayleigh number scales like `Ra ~ (wλ/κ)^(1/β)`
    Where `κ` is the thermal diffusivity and `β` is the convective beta parameter.

    This quantity must be taken lightly.

    Arguments:
    - `atmos::Atmos_t`      the atmosphere struct instance to be used.
    """
    function estimate_Ra!(atmos::atmosphere.Atmos_t)

        # Thermal diffusivity array
        κ::Array{Float64,1} = zero(atmos.layer_cp)
        @. κ = phys.calc_therm_diffus(atmos.layer_kc, atmos.layer_ρ, atmos.layer_cp)

        # One over beta
        ooβ::Float64 = 1.0 / phys.βRa

        # Estimate Rayleigh number
        @inbounds for i in 1:atmos.nlev_c
            atmos.diagnostic_Ra[i] = ( atmos.w_conv[i] * atmos.λ_conv[i] / κ[i]) ^ ooβ
        end

        return nothing
    end

    """
    **Estimate a diagnostic radiative timescale in each layer.**

    This quantity must be taken lightly.

    Arguments:
    - `atmos::Atmos_t`      the atmosphere struct instance to be used.
    """
    function estimate_timescale_rad!(atmos::atmosphere.Atmos_t)

        # Equation 10.1 from Seager textbook
        @inbounds for i in 1:atmos.nlev_c
            atmos.timescale_rad[i] = atmos.layer_cp[i] * (atmos.pl[i+1] - atmos.pl[i]) /
                                     (atmos.g[i] * 4 * phys.σSB * atmos.tmp[i]^3)
        end

        return nothing
    end

    """
    **Estimate a diagnostic convective timescale in each layer.**

    This quantity must be taken lightly.

    Arguments:
    - `atmos::Atmos_t`      the atmosphere struct instance to be used.
    """
    function estimate_timescale_conv!(atmos::atmosphere.Atmos_t)

        @inbounds for i in 1:atmos.nlev_c
            atmos.timescale_conv[i] = atmos.λ_conv[i] / max(atmos.w_conv[i], 1e-300)
        end

        return nothing
    end

end

