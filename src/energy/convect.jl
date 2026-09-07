# This file is part of AGNI. License is Apache-2.0: https://apache.org/licenses/LICENSE-2.0

"""
**Convection and Kzz submodule of `energy`.**
"""
module convect

    # System libraries
    using Logging

    # Local files
    import ..atmosphere
    import ..phys

    # Constants
    const CONVECT_MIN_PRESSURE::Float64   = 1e-9      # lowest pressure at which convection is allowed [bar]
    const CONVECT_REAL_GAS::Bool          = false     # use real gas EOS in convection scheme, if RG EOS enabled

    """
    **Calculate dry convective fluxes using mixing length theory.**

    Convective energy transport fluxes are calculated at every level edge, just
    like the radiative fluxes. This is not compatible with moist convection. By
    using MLT to parameterise convection, we can also calculate Kzz directly.

    Uses the mixing length formulation outlined by Joyce & Tayar (2023), which
    was also implemented in Lee et al. (2024), and also partially outlined in the review by
    Robinson & Marley (2014).
    https://arxiv.org/abs/2303.09596
    https://doi.org/10.1093/mnras/stae537
    https://ui.adsabs.harvard.edu/abs/1962JGR....67.3095B/abstract

     The adiabatic lapse rate is formulated as:
        `∇_ad = dln(T)/dln(P) = (P/T)*(dT/dP) = (P/T)*(1/[ρ c_p])`
    for an ideal gas, this becomes:
        `∇_ad = R / (μ c_p)`

    The mixing length is set to asymptotically approach H (for z>>H) or z (for
    z<H) as per Blackadar (1962). Alternatively, it can be set equal to H.
    https://doi.org/10.1029/JZ067i008p03095

    The scale height is formulated as:
        `Hp = P / (ρ g)`
    Where ρ is obtained from the equation of state.

    To account for convective stability due to compositional gradients, we can use the
    Ledoux criterion rather than the Schwarzschild criterion. This is described nicely in
    Gabriel et al. (2014), as well as Salaris & Cassisi (2017).
    http://dx.doi.org/10.1051/0004-6361/201423442
    https://doi.org/10.1098/rsos.170192

    In the ideal gas regime, the Ledoux criterion can be simply written as:
        `∇_ld = ∇_ad + dln(μ)/dln(P)`
    Using Equation (10) of Gabriel+14. Taking β=Pg/P=1 means the gas pressure equals the
    total pressure, neglecting pressure contributions from ions and electrons.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.

    Returns:
    - `Bool`                    whether the calculation succeeded
    """
    function convection!(atmos::atmosphere.Atmos_t)::Bool

        # Reset arrays
        fill!(atmos.mask_c,     false)
        fill!(atmos.flux_cdry,  0.0)
        fill!(atmos.Kzz,        0.0)
        fill!(atmos.λ_conv,     0.0)
        fill!(atmos.w_conv,     0.0)

        # Work variables
        Hp::Float64 = 0.0; hgt::Float64 = 0.0
        m1::Float64 = 0.0; m2::Float64 = 0.0; mt::Float64 = 0.0
        mu::Float64 = 0.0; c_p::Float64 = 0.0; rho::Float64 = 0.0
        ∇_ad::Float64 = 0.0; ∇_pr::Float64 = 0.0; ∇_μ::Float64 = 0.0; staby::Float64 = 0.0


        # Loop from bottom upwards (over cell-edges)
        for i in range(start=atmos.nlev_l-1, step=-1, stop=2)

            # Optionally skip low pressures
            if atmos.pl[i] <= CONVECT_MIN_PRESSURE * 1.0e5  # convert bar to Pa
                break
            end

            # Profile lapse rate: d(ln T)/d(ln P) = (P/T)*(dT/dP)
            ∇_pr = log(atmos.tmp[i-1]/atmos.tmp[i]) / log(atmos.p[i-1]/atmos.p[i])

            # Mass weights
            m1 = atmos.layer_σ[i-1]
            m2 = atmos.layer_σ[i]
            mt = m1+m2

            # Normalise weights
            m1 = m1/mt
            m2 = m2/mt

            # Properties interpolated to layer edge
            mu   = atmos.layer_μ[i]    * m2 + atmos.layer_μ[i-1]    * m1
            c_p  = atmos.layer_cp[i]   * m2 + atmos.layer_cp[i-1]   * m1
            rho  = atmos.layer_ρ[i]    * m2 + atmos.layer_ρ[i-1]    * m1
            Hp   = atmos.layer_Hp[i]   * m2 + atmos.layer_Hp[i-1]   * m1

            # Dry convective lapse rate
            if atmos.real_gas && CONVECT_REAL_GAS
                # general solution
                ∇_ad = atmos.pl[i] / (atmos.tmpl[i] * rho * c_p)
            else
                # ideal gas solution
                ∇_ad = (phys.R_gas / mu) / c_p
            end

            # Calculate lapse rate deviation from stability
            if atmos.mlt_criterion == 's'
                # Schwarzschild
                staby = ∇_pr - ∇_ad
            else
                # Ledoux is the only other option, for now
                ∇_μ = log(atmos.layer_μ[i-1]/atmos.layer_μ[i]) / log(atmos.p[i-1]/atmos.p[i])
                staby = ∇_pr - ∇_ad - ∇_μ
            end

            # Check instability
            if staby > 0

                atmos.mask_c[i] = true

                # Calculate the mixing length
                if !atmos.mlt_asymptotic
                    # Fixed
                    atmos.λ_conv[i] = phys.αMLT * Hp
                else
                    # Asymptotic
                    hgt = atmos.rl[i] - atmos.rp # height above the ground
                    atmos.λ_conv[i] = phys.k_vk * hgt / (1 + phys.k_vk * hgt/(phys.αMLT*Hp))
                end

                # Characteristic velocity (from Brunt-Vasalla frequency of parcel)
                atmos.w_conv[i] = atmos.λ_conv[i] * sqrt(atmos.gl[i]/Hp * staby)

                # Dry convective flux
                atmos.flux_cdry[i] = 0.5*rho*c_p*atmos.w_conv[i] * atmos.tmpl[i] * (atmos.λ_conv[i]/Hp) * staby

                # Kzz calculation
                if atmos.Kzz_type == 1
                    # Constant value
                    atmos.Kzz[i] = atmos.Kzz_kbreak
                elseif atmos.Kzz_type == 2
                    # Simple scaling
                    atmos.Kzz[i] = atmos.λ_conv[i] * atmos.w_conv[i]
                elseif atmos.Kzz_type == 3
                    # Eq16 from Charnay+15
                    atmos.Kzz[i] = (Hp/3.0) * (atmos.λ_conv[i]/Hp)^(4.0/3.0) * (phys.R_gas*atmos.flux_cdry[i]/(mu*rho*c_p))^(1.0/3.0)
                else
                    @warn "Invalid Kzz_type parameter: $(atmos.Kzz_type)"
                    return false
                end
            end
        end

        # Set surface quantities
        atmos.w_conv[end]    = 0.0
        atmos.λ_conv[end]    = 0.0
        atmos.flux_cdry[end] = 0.0

        return true
    end # end of mlt

    """
    **Fill Kzz values for remaining regions of profile.**

    This function is called after the convective fluxes have been calculated.
    The Kzz value in the convective regions are already calculated in the MLT scheme.

    This function calculates Kzz in the non-convective regions, by extending
    the Kzz values from the convective region through various parameters.

    Arguments:
    - `atmos::Atmos_t`      the atmosphere struct instance to be used.

    Returns:
    - `Bool`                function executed successfully
    """
    function fill_Kzz!(atmos::atmosphere.Atmos_t)::Bool

        # Temporary value
        Kzz_min::Float64  = atmos.Kzz_floor

        # Near-zero value
        Kzz_eps::Float64 = 1.0e-20

        # Find reference index for extension of Kzz, starting from convective regions
        i_Kzz_top::Int64 = atmos.nlev_l # default
        i_Kzz_bot::Int64 = atmos.nlev_l # default
        if any(atmos.Kzz .> Kzz_eps)
            # set to top of convective region
            i_Kzz_top = findfirst(x -> x > Kzz_eps, atmos.Kzz)
            # set to bottom of convective region
            i_Kzz_bot = findlast(x -> x > Kzz_eps, atmos.Kzz)
            # get minimum value, for filling intermediate zones
            Kzz_min = max(minimum(atmos.Kzz[atmos.Kzz .> Kzz_eps]), Kzz_min)
        else
            # otherwise, set to reference pressure
            i_Kzz_top = findmin(abs.(atmos.pl .- atmos.Kzz_pbreak))[2]
            i_Kzz_bot = i_Kzz_top
            atmos.Kzz[i_Kzz_top] = atmos.Kzz_kbreak
        end

        # Set zero-regions to minimum finite value
        # This covers the scenarios where multiple detatched convective regions exist
        atmos.Kzz[atmos.Kzz .<= Kzz_eps] .= Kzz_min

        # In regions above reference point, extend with power-law scaling.
        #   See equation 28 in Tsai+2020
        #       https://iopscience.iop.org/article/10.3847/1538-4357/ac29bc/pdf
        #   See also Charnay+15, Moses+16
        atmos.Kzz[1:i_Kzz_top] .= atmos.Kzz[i_Kzz_top] .* (  atmos.pl[1:i_Kzz_top] ./ atmos.pl[i_Kzz_top]) .^ atmos.Kzz_power

        # Extend Kzz downwards with constant value
        atmos.Kzz[i_Kzz_bot:end] .= atmos.Kzz[i_Kzz_bot]

        return true
    end

    export convection!, fill_Kzz!

end
