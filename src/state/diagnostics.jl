# This file is part of AGNI. License is Apache-2.0: https://apache.org/licenses/LICENSE-2.0
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
            atmos.timescale_conv[i] = atmos.λ_conv[i] / max(atmos.w_conv[i], eps(Float32))
        end

        return nothing
    end

end

