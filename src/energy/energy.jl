# This file is part of AGNI. License is Apache-2.0: https://apache.org/licenses/LICENSE-2.0


"""
**Contains the energy module, for everything relating to energy transport**

Calculates radiative, convective, conductive, etc, flux terms. Combines them together
self-consistently to update the temperature profile.
"""
module energy

    # System libraries
    using Printf
    using LinearAlgebra
    using Logging

    # Local files
    import ..atmosphere
    import ..phys
    import ..chemistry
    import ..spectrum
    import ..species

    # Constants
    const MIN_SKIN_D::Float64             = 1e-6      # minimum skin depth for conductive flux calculation [m]
    const MAX_SKIN_D::Float64             = 1e6       # maximum skin depth for conductive flux calculation [m]
    const ROUGHNESS_EPS::Float64          = 1e-3      # avoid blow-up of exchange coefficient when height ≈ roughness

    include("radtrans.jl"); import .radtrans: _make_finite!, radtrans!
    include("convect.jl"); import .convect: convection!, fill_Kzz!

    """
    **Calculate turbulent kinetic energy (TKE) exchange coefficient**.

    Based on Monin-Obukhov similarity theory, from roughness length scale.
    See eq 9 in Nicholson & Benn (2006). Added small epsilon-factor to avoid function
    blowing-up around regime where height ≈ roughness.

    A reasonable length scale can be found here: https://arxiv.org/pdf/2608.21549

    Arguments:
    - `height::Float64`     Height above surface [m]
    - `roughness::Float64`  Roughness length scale [m]

    Returns:
    - `C_d::Float64`        TKE exchange coefficient [dimensionless]
    """
    function eval_exchange_coeff(height::Float64, roughness::Float64)::Float64
        return phys.k_vk^2 / log(max(height, roughness+ROUGHNESS_EPS)/roughness)
    end

    """
    **Calculate sensible heat flux from turbulent kinetic energy (TKE)**

    Updates the values of `atmos.C_d` and `atmos.flux_sens`.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used
    """
    function sensible!(atmos::atmosphere.Atmos_t)::Bool

        # Set TKE exchange coefficient
        atmos.C_d = eval_exchange_coeff(atmos.r[end]-atmos.rp, atmos.surf_roughness)


        # TKE scheme for this 1D case
        # transports energy from the surface to the bottom node
        atmos.flux_sens = atmos.layer_cp[end]*atmos.layer_μ[end]*
                            atmos.p[end]/(phys.R_gas*atmos.tmp[end]) *
                            atmos.C_d * atmos.surf_windspeed *
                            (atmos.tmp_surf-atmos.tmp[end])
        return true
    end


    """
    **Calculate conductive heat fluxes using Fourier's law**

    Updates array of `atmos.flux_cdct` at each layer of the atmosphere.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used
    """
    function conduct!(atmos::atmosphere.Atmos_t)::Bool
        # top layer (to space)
        atmos.flux_cdct[1] = 0.0

        # bulk layers
        @inbounds for i in 2:atmos.nlev_l-1
            atmos.flux_cdct[i] = atmos.layer_kc[i] * (atmos.tmp[i]-atmos.tmp[i-1]) /
                                                        atmos.layer_thick[i]
        end

        # bottom layer (from surface)
        atmos.flux_cdct[end] = atmos.layer_kc[end] * (atmos.tmp[end]-atmos.tmp_surf) /
                                                      (atmos.r[end] - atmos.rp)
        return true
    end


    """
    **Calculate deep atmospheric heating flux.**

    The heating is deposited as a Gaussian distribution in log-pressure space,
    centered at `P_dep` with width `sigma_P`.

    Two power modes are supported:
    - `"rel"`   total flux = `deepheat_flux_rel * instellation` (stellar efficiency)
    - `"abs"`   total flux = `deepheat_flux_abs` (fixed radiative flux in W m⁻²)

    The flux gradient is defined as:
    dF_deep/dP = F_total / (sqrt(2π) * σ_P * P) * exp(-(ln(P) - ln(P_dep))² / (2 * σ_P²))

    This flux is integrated from the TOA downwards to obtain the cumulative
    flux at each cell edge, representing energy being deposited into the atmosphere.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    """
    function deep_heating!(atmos::atmosphere.Atmos_t)::Bool

        # Reset flux array
        fill!(atmos.flux_deep, 0.0)

        # Extract parameters
        sigma_P::Float64        = atmos.deepheat_Pwid
        below_domain::Bool      = atmos.deepheat_domain == "boundary_flux"

        # Determine total deposited flux [W m-2]
        F_total::Float64 = 0.0
        if atmos.deepheat_power_mode == "rel"
            # Stellar efficiency: define heating as a fraction of instellation (per unit area).
            F_total = atmos.deepheat_flux_rel * atmos.instellation
        elseif atmos.deepheat_power_mode == "abs"
            # Fixed radiative flux [W m-2]
            F_total = atmos.deepheat_flux_abs
        elseif atmos.deepheat_power_mode == "off"
            # No deep heating
            return true
        else
            @warn "Invalid deep heating power_mode: $(atmos.deepheat_power_mode)"
            return false
        end

        # Invert total flux for this internal-part of the calculation.
        #     This means that the deep heating flux will represent an additional source.
        #     Increasing the deep_heating then means that radiative (etc) fluxes will have
        #     to increase, to achieve the same total flux.
        F_total *= -1.0

        # If deposition is outside the domain and requested, apply as a bottom boundary flux.
        if below_domain && !( atmos.p_boa > atmos.deepheat_Pmid > atmos.p_toa )
            fill!(atmos.flux_deep, F_total)
            return true
        end

        # Prepare log-pressure variables
        log_Pmid::Float64 = log10(atmos.deepheat_Pmid)
        log_P::Array{Float64,1} = log10.(atmos.p)

        # Integrate from TOA downwards to get cumulative flux at each level edge
        atmos.flux_deep[1] = 0.0

        if atmos.deepheat_norm_method == "pressure"
            # Legacy: pressure-normalised dF/dP profile
            norm_factor::Float64 = 1.0 / (sqrt(2.0 * π) * sigma_P)
            dF_dP::Float64 = 0.0
            gaussian::Float64 = 0.0

            @inbounds for i in 1:atmos.nlev_c
                gaussian = exp(-(log_P[i] - log_Pmid)^2 / (2.0 * sigma_P^2))
                dF_dP = F_total * norm_factor * gaussian / atmos.p[i]
                atmos.flux_deep[i+1] = atmos.flux_deep[i] + dF_dP * (atmos.pl[i+1] - atmos.pl[i])
            end

        elseif atmos.deepheat_norm_method == "mass"
            # dm-weighted normalisation: ensures Σ(ε_dep*dm) = F_total
            # Column mass per unit area for layer i: dm_i = dp_i / g_i  [kg m⁻²]
            denom::Float64 = 0.0
            G::Float64 = 0.0
            dm_i::Float64 = 0.0

            @inbounds for i in 1:atmos.nlev_c
                G = exp(-(log_P[i] - log_Pmid)^2 / (2.0 * sigma_P^2))
                dm_i = (atmos.pl[i+1] - atmos.pl[i]) / atmos.g[i]
                denom += G * dm_i
            end

            if denom <= 1e-10
                @warn "Deep heating normalisation factor has non-positive denominator: $denom"
                return false
            end

            scale::Float64 = F_total / denom
            @inbounds for i in 1:atmos.nlev_c
                G = exp(-(log_P[i] - log_Pmid)^2 / (2.0 * sigma_P^2))
                dm_i = (atmos.pl[i+1] - atmos.pl[i]) / atmos.g[i]
                atmos.flux_deep[i+1] = atmos.flux_deep[i] + (scale * G * dm_i)
            end

        else
            @warn "Invalid deep heating normalisation: $(atmos.deepheat_norm_method)"
            return false
        end

        return true
    end


    """
    **Analytical diffusion scheme for condensation and evaporation energy.**

    Updates fluxes. Requires `chemistry._sat_aloft` to be called first.

    Integrates from bottom of model upwards. Based on the amount of
    phase change at each level, a phase change flux is calculated by assuming
    a fixed condensation timescale.

    If evaporation is enabled, then integrates from top downwards to determine flux from
    re-evaporation of droplets. Any droplets which reach the ground go towards forming an ocean.

    Should ideally perform a microphysical treatment; e.g. by following this paper:
    https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2020JE006653

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.

    Returns:
    - `Bool`                    whether the calculation succeeded
    """
    function latent!(atmos::atmosphere.Atmos_t)::Bool

        # Check if there are no condensates enabled
        if !atmos.condense_any
            return true
        end

        fill!(atmos.flux_l, 0.0)
        fill!(atmos.mask_l, false)

        # For each condensable
        for c in atmos.condensates

            # reset df,fl for this condensable
            fill!(atmos.phs_wrk_df,0.0)
            fill!(atmos.phs_wrk_fl,0.0)

            # Loop from top to bottom
            for i in 1:atmos.nlev_c-1

                # Skip bottom-most layer. Condensation at the surface is assumed to be in
                #    eqm with a surface ocean, so easier to assume there's no significant
                #    energy exchange, otherwise we get weird behaviour in the energy balance.

                # Calculate latent heat release at this level from the contributions
                #   of condensation (+) and evaporation (-), and a fixed timescale.
                atmos.phs_wrk_df[i] += species.get_Lv(atmos.gas_dat[c], atmos.tmp[i]) *
                                    (atmos.cond_yield[c][i] / atmos.phs_timescale)

            end # go to next level

            # Convert divergence to cell-edge fluxes.
            #     Assuming zero condensation at TOA, integrating downwards
            for i in 1:atmos.nlev_c
                atmos.phs_wrk_fl[i+1] = max(atmos.phs_wrk_df[i] + atmos.phs_wrk_fl[i], 0.0)
            end

            # Ensure that flux is zero at bottom of dry region.
            for i in 1:atmos.nlev_c
                # check for where no phase change is occuring below this level
                if maximum(abs.(atmos.phs_wrk_df[i:end])) < 1.0e-3
                    # if so, set all phase change fluxes to zero in that region
                    atmos.phs_wrk_fl[i+1:end] .= 0.0
                    break
                end
            end

            # add energy from this condesable to total energy from all condensables
            @. atmos.flux_l += atmos.phs_wrk_fl

            # calculate mask
            @. atmos.mask_l = (abs(atmos.flux_l) > eps(Float32))

        end # go to next condensable

        return true
    end

    """
    **Calculate conductive flux carried by conductive skin boundary layer.**

    This is a simple implementation of fourier's conduction law, with fixed conductivity
    and thickness of the boundary layer. Parameters are set in the atmos struct.

    F = k * ΔT / d

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.

    Returns:
    - `flux::Float64`           conductive flux through the skin layer [W m-2].
    """
    function skin_flux(atmos::atmosphere.Atmos_t)::Float64
        return (atmos.tmp_magma - atmos.tmp_surf) * atmos.skin_k / atmos.skin_d
    end

    """
    **Calculate conductive skin boundary layer thickness, given a conductive flux.**

    This is effectively the inverse of the `skin_flux` function, and can be used to
    determine the thickness of the boundary layer.

    d = k * ΔT / F

    Arguments:
    - `atmos::Atmos_t`      the atmosphere struct instance to be used.
    - `flux_skn::Float64`   the conductive flux through the skin layer [W m-2].

    Returns:
    - `skin_d::Float64`     conductive skin boundary layer thickness [m].
    """
    function skin_depth(atmos::atmosphere.Atmos_t, flux_skn::Float64)::Float64
        return clamp(
                    (atmos.tmp_magma - atmos.tmp_surf) * atmos.skin_k / flux_skn,
                    MIN_SKIN_D, MAX_SKIN_D
                )
    end

    """
    **Reset energy fluxes to zero.**
    """
    function reset_fluxes!(atmos::atmosphere.Atmos_t)::Bool

        # sensible heating
        atmos.flux_sens = 0.0

        # conduct
        fill!(atmos.flux_cdct, 0.0)

        # convect
        fill!(atmos.flux_cdry, 0.0)

        # latent heating
        fill!(atmos.flux_l, 0.0)

        # radiative (bolometric)
        atmos.is_out_sw = false
        atmos.is_out_lw = false
        fill!(atmos.flux_u, 0.0)
        fill!(atmos.flux_d, 0.0)
        fill!(atmos.flux_n, 0.0)
        fill!(atmos.flux_l, 0.0)
        fill!(atmos.flux_n_lw, 0.0)
        fill!(atmos.flux_n_sw, 0.0)
        fill!(atmos.flux_u_lw, 0.0)
        fill!(atmos.flux_u_sw, 0.0)
        fill!(atmos.flux_d_sw, 0.0)
        fill!(atmos.flux_d_lw, 0.0)

        # radiative (per band)
        fill!(atmos.band_u_lw, 0.0)
        fill!(atmos.band_d_lw, 0.0)
        fill!(atmos.band_n_lw, 0.0)
        fill!(atmos.band_u_sw, 0.0)
        fill!(atmos.band_d_sw, 0.0)
        fill!(atmos.band_n_sw, 0.0)

        # deep heating
        fill!(atmos.flux_deep, 0.0)

        # total fluxes, and difference across each layer
        fill!(atmos.flux_tot, 0.0)
        fill!(atmos.flux_dif, 0.0)

        return true
    end


    """
    **Calculate energy flux at each level.**

    Calculates flux components (radtrans, convection, etc.) and sums them to get total flux.
    Also updates thermodynamic properties (heat capacity, density, etc.) at each layer.

    Assumes that chemistry functions have already been called, if wanted. Does not call
    fastchem here.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere instance to be used.

    Optional arguments:
    - `radiative::Bool`                 include radiation fluxes
    - `latent_heat::Bool`               include condensation flux
    - `convective::Bool`                include MLT convection flux
    - `sens_heat::Bool`                 include TKE sensible heat flux
    - `conductive::Bool`                include conductive heat flux
    - `deep::Bool`                      include deep heating flux (layer-internal production)
    - `convect_sf::Float64`             scale factor applied to convection fluxes
    - `latent_sf::Float64`              scale factor applied to phase change fluxes
    - `calc_cf::Bool`                   calculate LW contribution function?
    - `calc_hr::Bool`                   calculate heating rates from fluxes?

    Returns:
    - `Bool`                            calculation succeeded
    """
    function calc_fluxes!(atmos::atmosphere.Atmos_t;
                          radiative::Bool=false, latent_heat::Bool=false, convective::Bool=false,
                          sens_heat::Bool=false, conductive::Bool=false, deep::Bool=false,
                          convect_sf::Float64=1.0, latent_sf::Float64=1.0,
                          calc_cf::Bool=false, calc_hr::Bool=false)::Bool


        # Reset fluxes
        reset_fluxes!(atmos)
        ok::Bool = true

        # Warn if no flux terms are enabled
        if !(radiative || latent_heat || convective || sens_heat || conductive || deep)
            @warn "No flux terms enabled in call to `calc_fluxes!`"
            ok = false
        end

        # +Latent heating
        if latent_heat
            ok &= latent!(atmos)           # Calculate latent heat fluxes
            atmos.flux_l *= latent_sf           # Modulate for stability?
            @. atmos.flux_tot += atmos.flux_l   # Add to total flux
        end

        # +Radiation
        if radiative
            ok &= radtrans!(atmos, true, calc_cf=calc_cf)   # Longwave
            ok &= radtrans!(atmos, false)                   # Shortwave
            @. atmos.flux_tot += atmos.flux_n  # Add to total flux
        end

        # +Dry convection
        if convective
            ok &= convection!(atmos)                          # Calc dry convection heat flux
            atmos.flux_cdry *= convect_sf               # Modulate for stability?
            @. atmos.flux_tot += atmos.flux_cdry # Add to total flux
        end

        # Calculate Kzz in non-convective regions
        fill_Kzz!(atmos)

        # +Surface turbulence
        if sens_heat
            ok &= sensible!(atmos)
            atmos.flux_tot[end] += atmos.flux_sens
        end

        # +Conduction
        if conductive
            ok &= conduct!(atmos)
            @. atmos.flux_tot += atmos.flux_cdct
        end

        # +Deep atmospheric heating
        if deep
            ok &= deep_heating!(atmos)
            @. atmos.flux_tot += atmos.flux_deep
        end

        # Flux difference across each level
        # Positive value => heating
        atmos.flux_dif[1:end] .= atmos.flux_tot[2:end] .- atmos.flux_tot[1:end-1]

        # Heating rate
        if calc_hr
            ok &= calc_hrates!(atmos)
        end

        return ok
    end

    """
    **Calculate heating rates at cell-centres from the total flux.**

    Requires the total flux to have already been set (atmos.flux_dif). Heating
    rates are calculated in units of kelvin per day.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.

    Returns:
    - `Bool`                            whether the calculation succeeded
    """
    function calc_hrates!(atmos::atmosphere.Atmos_t)::Bool
        # Ensure flux difference has been calculated
        atmos.flux_dif[1:end] .= atmos.flux_tot[2:end] .- atmos.flux_tot[1:end-1]

        # Evaluate heating rate
        for i in 1:atmos.nlev_c
            atmos.heating_rate[i] = (atmos.a[i] / atmos.layer_cp[i]) *
                                        atmos.flux_dif[i] / (atmos.pl[i+1] - atmos.pl[i])
        end

        atmos.heating_rate *= 86400.0 # K/day

        return true
    end

end # end module
