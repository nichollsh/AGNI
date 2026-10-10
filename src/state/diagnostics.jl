# This file is part of AGNI. License is Apache-2.0: https://apache.org/licenses/LICENSE-2.0
module diagnostics

    import ..phys
    import ..atmosphere
    import ..consts: k_B, R_gas

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

    """
    **Ratio of mean free path to scale height, in one layer.**

    Uses the Maxwell mean free path `l = 1 / (sqrt(2) n σ)`, with `n = p / (k_B T)`

    Mole-fraction-weighted hard-sphere cross-section `σ = Σ x_j π d_j²` from collision diameter.

    Arguments:
        - `atmos::Atmos_t`      the atmosphere struct instance to be used.
        - `i::Int64`            index of the layer

    Returns:
        - `ratio::Float64`      l / H, dimensionless (Inf if no gas is present)
        - `sigma::Float64`      mixture collision cross-section [m2]
        - `m_bar::Float64`      mean particle mass [kg]
    """
    function _mfp_over_H(atmos::atmosphere.Atmos_t, i::Int64)::Tuple{Float64,Float64,Float64}

        # Compute mixture collision cross-section and mean particle mass
        x_tot::Float64 = 0.0
        sigma::Float64 = 0.0
        for gas in atmos.gas_names
            x = atmos.gas_vmr[gas][i]
            x > 0.0 || continue # skip vmr=0
            x_tot += x
            sigma += x * pi * atmos.gas_dat[gas].particle_d^2 # add up xsec weighted by VMR
        end

        # if xsec=0, return Inf for ratio, and 0 for sigma and m_bar
        if (x_tot <= 0.0) || (sigma <= 0.0)
            return (Inf, 0.0, 0.0)
        end
        sigma /= x_tot

        # average mmw
        m_bar::Float64 = atmos.layer_μ[i] * k_B / R_gas

        # particle number density
        n::Float64     = atmos.p[i] / (k_B * atmos.tmp[i])

        # mean free path and scale height
        mfp::Float64   = 1.0 / (sqrt(2.0) * n * sigma)
        H::Float64     = k_B * atmos.tmp[i] / (m_bar * atmos.g[i])

        # return ratio
        return (mfp / H, sigma, m_bar)
    end

    """
    **Locate the exobase (where the MFP equals scale height).**

    The exobase is the first crossing of l/H = 1.
    If the modelled column does not reach the exobase, its pressure is set to the TOA.

    Arguments:
        - `atmos::Atmos_t`      the atmosphere struct instance to be used.

    Returns:
        - `in_domain::Bool`     whether the exobase lies within the modelled column
    """
    function estimate_exobase!(atmos::atmosphere.Atmos_t)::Bool

        # default to above TOA
        atmos.exobase_p         = atmos.pl[1]
        atmos.exobase_r         = atmos.rl[1]
        atmos.exobase_tmp       = atmos.tmpl[1]
        atmos.exobase_in_domain = false

        # Scan upwards from the bottom layer
        for i in range(start=atmos.nlev_c, stop=1, step=-1)
            ratio, _, _ = _mfp_over_H(atmos, i)
            if ratio >= 1.0
                atmos.exobase_in_domain = true
                atmos.exobase_p         = atmos.p[i]
                atmos.exobase_r         = atmos.r[i]
                atmos.exobase_tmp       = atmos.tmp[i]
                return true
            end
        end

        # Did not reach exobase, so leave at TOA
        return false
    end

end
