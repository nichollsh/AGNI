# Contains the energy module, for everything relating to energy transport

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module energy

    # System libraries
    using Printf
    using LinearAlgebra
    using Logging
    using LoopVectorization

    # Local files
    import ..atmosphere
    import ..phys
    import ..chemistry
    import ..spectrum

    """
    **Calculate radiative fluxes using SOCRATES.**

    Uses the configuration inside the atmos struct. Can either do LW or SW
    calculation as required. Imports SOCRATES wrapper from the atmosphere
    module, rather than loading it twice.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `lw::Bool`                        longwave calculation? Else: shortwave
    - `calc_cf::Bool=false`             also calculate contribution function?
    """
    function radtrans!(atmos::atmosphere.Atmos_t, lw::Bool; calc_cf::Bool=false)
        if !atmos.is_alloc
            error("atmosphere arrays have not been allocated")
        end
        if !atmos.is_param
            error("atmosphere parameters have not been set")
        end

        atmos.num_rt_eval += 1

        # Longwave or shortwave calculation?
        # Set the two-stream approximation to be used (-t f)
        if lw
            atmos.control.isolir = atmosphere.SOCRATES.rad_pcf.ip_infra_red
            atmos.control.i_2stream = atmosphere.SOCRATES.rad_pcf.ip_elsasser
            # Practical improved flux method (1985) with Elsasser's diffusivity (D=1.66)
        else
            atmos.control.isolir = atmosphere.SOCRATES.rad_pcf.ip_solar
            atmos.control.i_2stream = atmosphere.SOCRATES.rad_pcf.ip_pifm80
            # Practical improved flux method (original form of 1980)
        end

        # Check files are acceptable and set instellation if doing SW flux
        if lw
            if !Bool(atmos.spectrum.Basic.l_present[6])
                error("The spectral file contains no data for the Planck function.
                       Check that the file contains a stellar spectrum.")
            end

            if Bool(atmos.spectrum.Basic.l_present[2])
                atmos.control.l_solar_tail_flux = true
            end

        else
            if !Bool(atmos.spectrum.Basic.l_present[2])
                error("The spectral file contains no solar spectral data.")
            end

            # Downward SW flux in atmosphere at TOA stored by AGNI
            atmos.toa_heating = atmos.instellation * (1.0 - atmos.albedo_b) *
                                    atmos.s0_fact * cosd(atmos.zenith_degrees)

            # SOCRATES requires this to be passed as two variables, since it
            #     needs to know the angle of the direct beam.

            # Convert the zenith angles to secants.
            atmos.bound.zen_0[1] = 1.0/cosd(atmos.zenith_degrees)

            # Pass effective solar constant
            atmos.bound.solar_irrad[1] = atmos.instellation *
                                            (1.0 - atmos.albedo_b) * atmos.s0_fact
        end

        #####################################
        # Angular integration
        # see src/aux/angular_control_cdf.f
        #####################################

        # Cl_run_cdf +R flag
        atmos.control.l_rescale = false
        if atmos.control.l_rescale
            atmos.control.l_henyey_greenstein_pf = true
        end

        # The internal SOCRATES solver used for the two-stream calculations (-v flag)
        if atmos.control.l_cloud
            # 16 is recommended for cloudy-sky (ip_solver_mix_direct_hogan)
            # 17 is recommended for cloud with separate stratiform and convective regions
            atmos.control.i_solver = atmosphere.SOCRATES.rad_pcf.ip_solver_mix_direct_hogan
        else
            # 13 is recommended for clear-sky (ip_solver_homogen_direct)
            atmos.control.i_solver = atmosphere.SOCRATES.rad_pcf.ip_solver_homogen_direct
        end


        #      Arrays of fluxes must be of the full size.
        atmos.dimen.nd_2sg_profile =        atmos.dimen.nd_profile
        atmos.dimen.nd_flux_profile =       atmos.dimen.nd_profile
        atmos.dimen.nd_radiance_profile =   1
        atmos.dimen.nd_j_profile =          1
        atmos.dimen.nd_viewing_level =      1
        atmos.dimen.nd_sph_coeff =          1

        # Reset dimen.nd_max_order to reduce memory requirements
        atmos.dimen.nd_max_order = 1

        #####################################
        # Surface albedo
        #####################################

        fill!(atmos.bound.rho_alb, 0.0)
        atmos.bound.rho_alb[1, atmosphere.SOCRATES.rad_pcf.ip_surf_alb_diff, :] .= atmos.surf_r_arr
        atmos.bound.rho_alb[1, atmosphere.SOCRATES.rad_pcf.ip_surf_alb_dir,  :] .= atmos.surf_r_arr

        ###################################################
        # Cloud information
        ###################################################

        atmos.cld.w_cloud[1,:]               .= atmos.cloud_arr_f[:]
        atmos.cld.frac_cloud[1,:,1]          .= 1.0
        atmos.cld.condensed_mix_ratio[1,:,1] .= atmos.cloud_arr_l[:]
        atmos.cld.condensed_dim_char[1,:,1]  .= atmos.cloud_arr_r[:]

        ###################################################
        # Treatment of scattering
        ###################################################

        atmos.control.i_scatter_method = atmosphere.SOCRATES.rad_pcf.ip_scatter_full
        for i in atmos.control.first_band:atmos.control.last_band
            atmos.control.i_scatter_method_band[i] = atmos.control.i_scatter_method
        end

        ####################################################
        # Temperature
        ###################################################

        atmos.atm.p[1, :] .= atmos.p[:]
        atmos.atm.t[1, :] .= atmos.tmp[:]
        atmos.atm.p_level[1, 0:end] .= atmos.pl[:]
        atmos.atm.t_level[1, 0:end] .= atmos.tmpl[:]

        if lw
            atmos.bound.t_ground[1] = atmos.tmp_surf
        end

        if lw
            atmos.control.l_ir_source_quad = true
        end

        # Set flux in surface emission, by band
        #     Equal to integral of planck function over band width, which in
        #     this case is done by simply evaluating at the midpoint and
        #     multiplying by band width. Scaled by the emissivity.
        @. atmos.surf_flux = phys.evaluate_planck(atmos.bands_cen, atmos.tmp_surf) *
                                atmos.bands_wid * 1e9 * atmos.surf_e_arr

        # Pass to socrates array
        #     I would argue that the 1-albedo term shouldn't be here, but it is to correct
        #     for it also (strangely) appearing inside diff_planck_source_mod.f90 on
        #     line 129. Having this 1-albedo term (and using this low-order integration)
        #     gives the correct results from my tests versus SOCRATES's native function.
        @inbounds for i in 1:atmos.nbands
            atmos.bound.flux_ground[1,i] = atmos.surf_flux[i] * (1.0 - atmos.surf_r_arr[i])
        end

        ######################################################
        # Run radiative transfer model
        ######################################################

        # Calculate contribution function?
        atmos.control.l_contrib_func_band = calc_cf

        # Set composition for each gas,level
        for (i_gas,s_gas) in enumerate(atmos.gas_soc_names)
            for i in 1:atmos.nlev_c
                # skip unspecified gases
                if (s_gas in atmos.gas_names)
                    # convert VOLUME mixing ratio to MASS mixing ratio
                    atmos.atm.gas_mix_ratio[1, i, i_gas] = atmos.gas_vmr[s_gas][i] *
                                                            atmos.gas_dat[s_gas].mmw /
                                                            atmos.layer_μ[i]
                else
                    atmos.atm.gas_mix_ratio[1, i, i_gas] = 0.0
                end
                # do not normalise MMRs to 1
            end
        end


        # Do radiative transfer
        atmosphere.atmosphere.SOCRATES.radiance_calc(atmos.control,
                                                     atmos.dimen, atmos.spectrum,
                                                     atmos.atm, atmos.cld, atmos.aer,
                                                     atmos.bound, atmos.radout)

        # Check finite
        if !all(isfinite, atmos.radout.flux_down)
            if lw
                @error "Non-finite value in LW DN flux array"
            else
                @error "Non-finite value in SW DN flux array"
            end
        end
        if !all(isfinite, atmos.radout.flux_up)
            if lw
                @error "Non-finite value in LW UP flux array"
            else
                @error "Non-finite value in SW UP flux array"
            end
        end

        # Store new fluxes in atmos struct
        idx::Int = 1
        if lw
            # LW case
            for lv in 1:atmos.nlev_l      # sum over levels
                for ba in 1:atmos.dimen.nd_channel  # sum over bands
                    idx = lv+(ba-1)*atmos.nlev_l
                    atmos.band_d_lw[lv,ba] = max(0.0, atmos.radout.flux_down[idx])
                    atmos.band_u_lw[lv,ba] = max(0.0, atmos.radout.flux_up[idx])
                end
                atmos.flux_d_lw[lv] = sum(atmos.band_d_lw[lv,:])
                atmos.flux_u_lw[lv] = sum(atmos.band_u_lw[lv,:])
            end
            atmos.band_n_lw = atmos.band_u_lw - atmos.band_d_lw
            atmos.flux_n_lw = atmos.flux_u_lw - atmos.flux_d_lw

            # Contribution function (only LW stream contributes)
            fill!(atmos.contfunc_band,0.0)
            if calc_cf
                for ba in 1:atmos.dimen.nd_channel
                    for lv in 1:atmos.nlev_c
                        atmos.contfunc_band[lv,ba] = atmos.radout.contrib_funcf_band[1,lv,ba]
                    end
                end
            end
            atmos.is_out_lw = true
        else
            # SW case
            for lv in 1:atmos.nlev_l                # sum over levels
                for ba in 1:atmos.dimen.nd_channel  # sum over bands
                    idx = lv+(ba-1)*atmos.nlev_l
                    atmos.band_d_sw[lv,ba] = max(0.0,atmos.radout.flux_down[idx])
                    atmos.band_u_sw[lv,ba] = max(0.0,atmos.radout.flux_up[idx])
                end
                atmos.flux_d_sw[lv] = sum(atmos.band_d_sw[lv,:])
                atmos.flux_u_sw[lv] = sum(atmos.band_u_sw[lv,:])
            end
            atmos.band_n_sw = atmos.band_u_sw - atmos.band_d_sw
            atmos.flux_n_sw = atmos.flux_u_sw - atmos.flux_d_sw
            atmos.is_out_sw = true
        end

        # Store net fluxes when we have both SW and LW components
        if atmos.is_out_lw && atmos.is_out_sw
            atmos.flux_d = atmos.flux_d_lw + atmos.flux_d_sw
            atmos.flux_u = atmos.flux_u_lw + atmos.flux_u_sw
            atmos.flux_n = atmos.flux_n_lw + atmos.flux_n_sw
        end

        return nothing
    end # end of radtrans


    # Calculate sensible heat flux (turbulence at surface boundary)
    function sensible!(atmos::atmosphere.Atmos_t)
        # TKE scheme for this 1D case
        # transports energy from the surface to the bottom node
        atmos.flux_sens = atmos.layer_cp[end]*atmos.layer_μ[end]*
                            atmos.p[end]/(phys.R_gas*atmos.tmp[end]) *
                            atmos.C_d * atmos.U *
                            (atmos.tmp_surf-atmos.tmp[end])
        return nothing
    end


    # Calculate conductive fluxes
    function conduct!(atmos::atmosphere.Atmos_t)
        # top layer
        atmos.flux_cdct[1] = 0.0

        # bulk layers
        @inbounds for i in 2:atmos.nlev_l-1
            atmos.flux_cdct[i] = atmos.layer_kc[i] * (atmos.tmp[i]-atmos.tmp[i-1]) /
                                                        (atmos.r[i-1]-atmos.r[i])
        end

        # bottom layer
        atmos.flux_cdct[end] = atmos.layer_kc[end] * (atmos.tmp[end]-atmos.tmp_surf) /
                                                      (atmos.r[end] - atmos.rp)
        return nothing
    end


    """
    **Calculate dry convective fluxes using mixing length theory.**

    Uses the mixing length formulation outlined by Joyce & Tayar (2023), which
    was also implemented in Lee et al. (2024), and partially outlined in an
    earlier paper by Robinson & Marley (2014).

    https://arxiv.org/abs/2303.09596
    https://doi.org/10.1093/mnras/stae537
    https://ui.adsabs.harvard.edu/abs/1962JGR....67.3095B/abstract

    Convective energy transport fluxes are calculated at every level edge, just
    like the radiative fluxes. This is not compatible with moist convection. By
    using MLT to parameterise convection, we can also calculate Kzz directly.

    The mixing length is set to asymptotically approach H (for z>>H) or z (for
    z<H) as per Blackadar (1962). Alternatively, it can be set equal to H.

    The scale height is formulated as:
        `Hp = P / (ρ g)`
    The adiabatic lapse rate is formulated as:
        `∇_ad = dln(T)/dln(P) = (P/T)*(dT/dP) = (P/T)*(1/[ρ c_p])`
    for an ideal gas, this becomes:
        `∇_ad = R / (μ c_p)`


    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    - `pmin::Float64`           pressure [bar] below which convection is disabled
    - `mltype::Int`             mixing length value (1: scale height, 2: asymptotic)
    """
    function convection!(atmos::atmosphere.Atmos_t; pmin::Float64=1.0e-4, mltype::Int=2)

        # Reset arrays
        fill!(atmos.mask_c,     false)
        fill!(atmos.flux_cdry,  0.0)
        fill!(atmos.Kzz,        0.0)

        # Work variables
        Hp::Float64 = 0.0; λ::Float64 = 0.0; w::Float64 = 0.0
        m1::Float64 = 0.0; m2::Float64 = 0.0; mt::Float64 = 0.0
        grav::Float64 = 0.0; mu::Float64 = 0.0; c_p::Float64 = 0.0; rho::Float64 = 0.0
        ∇_ad::Float64 = 0.0; ∇_pr::Float64 = 0.0; hgt::Float64 = 0.0

        # Loop from bottom upwards (over cell-edges)
        for i in range(start=atmos.nlev_l-1, step=-1, stop=2)

            # Profile lapse rate: d(ln T)/d(ln P) = (P/T)*(dT/dP)
            ∇_pr = ( log(atmos.tmp[i-1]/atmos.tmp[i]) )/( log(atmos.p[i-1]/atmos.p[i]) )

            # Optionally skip low pressures
            if atmos.pl[i] <= pmin * 1.0e5  # convert bar to Pa
                break
            end

            # Mass weights
            m1 = atmos.layer_mass[i-1]
            m2 = atmos.layer_mass[i]
            mt = m1+m2

            # Normalise weights
            m1 = m1/mt
            m2 = m2/mt

            # Properties interpolated to layer edge
            grav = atmos.layer_grav[i] * m2 + atmos.layer_grav[i-1] * m1
            mu   = atmos.layer_μ[i]    * m2 + atmos.layer_μ[i-1]    * m1
            c_p  = atmos.layer_cp[i]   * m2 + atmos.layer_cp[i-1]   * m1
            rho  = atmos.layer_ρ[i]    * m2 + atmos.layer_ρ[i-1]    * m1

            # Dry convective lapse rate, and pressure scale height
            if atmos.real_gas
                # general solution
                ∇_ad = atmos.pl[i] / (atmos.tmpl[i] * rho * c_p)
                Hp = atmos.pl[i] / (rho * grav)
            else
                # ideal gas solution
                ∇_ad = (phys.R_gas / mu) / c_p
                Hp = phys.R_gas * atmos.tmpl[i] / (mu * grav)
            end

            # Check instability
            if ∇_pr > ∇_ad

                atmos.mask_c[i] = true

                # Calculate the mixing length
                if mltype == 1
                    # Fixed
                    λ = phys.αMLT * Hp
                elseif mltype == 2
                    # Asymptotic
                    hgt = atmos.rl[i] - atmos.rp # height above the ground
                    λ = phys.k_vk * hgt / (1 + phys.k_vk * hgt/Hp)
                else
                    # Otherwise
                    error("Invalid mixing length type selected: $mltype")
                end

                # Characteristic velocity (from Brunt-Vasalla frequency of parcel)
                w = λ * sqrt(grav/Hp * (∇_pr-∇_ad))

                # Dry convective flux
                atmos.flux_cdry[i] = 0.5*rho*c_p*w * atmos.tmpl[i] * (λ/Hp)*(∇_pr-∇_ad)

                # Convection eddy diffusion coefficient [m2 s-1]
                atmos.Kzz[i] = w * λ

            end
        end

        # Set surface Kzz
        atmos.Kzz[end] = atmos.Kzz[end-1]

        # Check for spurious shallow convection occuring ABOVE condensing regions
        #    If found, reset convective flux to zero AT THIS LAYER ONLY.
        #    This is okay because this shouldn't physically happen, and will only occur
        #    because of weird numerical issues which only act to make solving difficult.
        @inbounds for i in 1:atmos.nlev_l-1
            if (!atmos.mask_l[i] && any(atmos.mask_l[i+1:end])) #|| (atmos.mask_l[i] && !atmos.mask_c[i-1] && !atmos.mask_c[i+1])
                atmos.mask_c[i] = false
                atmos.flux_cdry[i] = 0.0
            end
        end

        return nothing
    end # end of mlt

    """
    **Fill Kzz values for the entire profile.**

    This function is called after the convection scheme has been run.

    The Kzz value in the convective region (and below) are set equal to the maximum value
    in the convective region, as calculated by MLT. The value increases with power-law
    scaling with pressure in the stratosphere.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    """
    function fill_Kzz!(atmos::atmosphere.Atmos_t)

        Fthresh::Float64 = 1.0e-8

        # Kzz limits
        clamp!(atmos.Kzz, atmos.Kzz_floor, atmos.Kzz_ceiling)

        # Find bottom of convective region (looping downwards)
        i_cnvct_bot::Int = atmos.nlev_l
        @inbounds for i in 1:atmos.nlev_l-1
            if (atmos.flux_cdry[i+1] < Fthresh) && (atmos.flux_cdry[i] > Fthresh)
                i_cnvct_bot = i
            end
        end

        # Find top of convective region (looping upwards from bottom)
        i_cnvct_top::Int = i_cnvct_bot
        atmos.Kzz_pbreak = min(1e5, atmos.pl[end])
        atmos.Kzz_kbreak = atmos.Kzz_floor
        @inbounds for i in range(start=i_cnvct_bot, step=-1, stop=1)
            if (atmos.flux_cdry[i] < Fthresh) && (atmos.flux_cdry[i+1] > Fthresh)
                i_cnvct_top = i
            end
        end

        # Store breakpoint values
        atmos.Kzz_pbreak = atmos.pl[i_cnvct_top]
        atmos.Kzz_kbreak = maximum(atmos.Kzz[i_cnvct_top:end])

        # Set Kzz within and below convective region to constant value. This value best
        #    represents the diffusive mixing in this region.
        atmos.Kzz[i_cnvct_top:end] .= atmos.Kzz_kbreak

        # Extend Kzz in stratosphere based on power-law scaling.
        #   See equation 28 in Tsai+2020
        #   https://iopscience.iop.org/article/10.3847/1538-4357/ac29bc/pdf
        atmos.Kzz[1:i_cnvct_top-1] .= atmos.Kzz[i_cnvct_top] .*
                                        ( atmos.Kzz_pbreak ./ atmos.pl[1:i_cnvct_top-1]).^0.4

        # Kzz limits
        clamp!(atmos.Kzz, atmos.Kzz_floor, atmos.Kzz_ceiling)

        return nothing
    end

    """
    **Analytical diffusion scheme for condensation and evaporation energy.**

    Integrates from bottom of model upwards. Based on the amount of
    phase change at each level, a phase change flux is calculated by assuming
    a fixed condensation timescale.

    Updates fluxes. Requires `chemistry.handle_saturation` to be called first in the
    multi-component case.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    """
    function condense_diffuse!(atmos::atmosphere.Atmos_t)

        # Check if there are no condensates enabled
        if !atmos.condense_any
            return nothing
        end

        fill!(atmos.flux_l, 0.0)
        fill!(atmos.mask_l, false)

        # Single-gas relaxation case
        single::Bool = Bool(atmos.gas_num == 1)

        # Variables for tracking phase change energy
        evap_eff::Float64 =  0.1    # evaporation efficiency [dimensionless]
        evap_scl::Float64 =  1.0e-3 # relative increase in evap_eff w/ pressure [K-1]
        E_accum::Float64 =   0.0    # accumulated condensational energy [J]
        i_dry_top::Int =     1      # deepest point at which condensation occurs
        i_dry_bot::Int =     2      # deepest point at which criticality occurs

        # For each condensable
        for c in atmos.condensates

            # reset df,fl for this condensable
            fill!(atmos.phs_wrk_df,0.0)
            fill!(atmos.phs_wrk_fl,0.0)

            # Loop from bottom to top
            for i in 1:atmos.nlev_c

                if single
                    # --------------------------------
                    # Single-component relaxation scheme

                    # Check criticality
                    if atmos.tmp[i] > atmos.gas_dat[c].T_crit
                        continue
                    end

                    # check saturation
                    qsat = phys.get_Psat(atmos.gas_dat[c], atmos.tmp[i])/atmos.p[i]
                    if (atmos.gas_vmr[c][i] < qsat+1.0e-10)
                        continue
                    end

                    # relaxation function
                    atmos.phs_wrk_df[i] = (atmos.gas_vmr[c][i]-qsat)* atmos.layer_thick[i]*
                                          (atmos.pl[i+1]-atmos.pl[i])/ atmos.phs_tau_sgl

                    # flag layer as set by saturation
                    if atmos.phs_wrk_df[i] > 1.0e-10
                        atmos.gas_sat[c][i] = true
                    end

                else
                    # --------------------------------
                    # Multicomponent diffusion scheme

                    # Calculate latent heat release at this level from the contributions
                    #   of condensation (+) and evaporation (-), and a fixed timescale.
                    atmos.phs_wrk_df[i] += phys.get_Lv(atmos.gas_dat[c], atmos.tmp[i]) *
                                        (atmos.gas_yield[c][i] / atmos.phs_tau_mix)

                end

            end # go to next level

            # Evaporation flux ...

            # find top of dry region
            for i in 1:atmos.nlev_c
                if abs(atmos.phs_wrk_df[i]) > 0
                    i_dry_top=i+2
                end
            end

            # find bottom of dry region
            i_dry_bot = i_dry_top + 1
            for i in i_dry_top+1:atmos.nlev_c
                for c in atmos.condensates
                    if atmos.tmp[i] < atmos.gas_dat[c].T_crit
                        # this layer is not supercritical for this gas, so
                        # evaporative flux can be dissipated here
                        i_dry_bot = i
                        break # go to next layer (below)
                    end
                end
            end

            i_dry_bot = atmos.nlev_c

            # accumulated condensational energy
            E_accum = sum(atmos.phs_wrk_df[1:i_dry_top])

            # dissipate E_accum by evaporation in the dry region
            for i in i_dry_top+1:i_dry_bot

                if E_accum < 1.0e-7
                    # dissipate all of the flux
                    atmos.phs_wrk_df[i] = -E_accum
                    E_accum = 0.0
                    break
                else
                    # dissipate some fraction of the accumuated flux
                    atmos.phs_wrk_df[i] = -1.0 * evap_eff *E_accum

                    # update total energy budget
                    E_accum += atmos.phs_wrk_df[i]

                    # evaporation becomes increasingly efficient at hotter levels
                    evap_eff = min(1.0, evap_eff * (1 + evap_scl * atmos.tmp[i]))
                end

            end

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
            atmos.phs_wrk_fl[end] = 0.0

            # add energy from this gas to total
            @turbo @. atmos.flux_l += atmos.phs_wrk_fl

            # calculate mask
            @. atmos.mask_l = (abs(atmos.flux_l) > 1.0e-30)

        end # go to next condensable

        return nothing
    end

    """
    **Calculate flux carried by conductive skin.**
    """
    function skin_flux(atmos::atmosphere.Atmos_t)::Float64
        return (atmos.tmp_magma - atmos.tmp_surf) * atmos.skin_k / atmos.skin_d
    end

    """
    **Reset energy fluxes to zero.**
    """
    function reset_fluxes!(atmos::atmosphere.Atmos_t)

        # scalar fluxes
        atmos.flux_sens = 0.0

        # conduct
        fill!(atmos.flux_cdct, 0.0)

        # convect
        fill!(atmos.flux_cdry, 0.0)

        # radiative grey
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

        # radiative band
        fill!(atmos.band_u_lw, 0.0)
        fill!(atmos.band_d_lw, 0.0)
        fill!(atmos.band_n_lw, 0.0)
        fill!(atmos.band_u_sw, 0.0)
        fill!(atmos.band_d_sw, 0.0)
        fill!(atmos.band_n_sw, 0.0)

        # total fluxes
        fill!(atmos.flux_dif, 0.0)
        fill!(atmos.flux_tot, 0.0)
    end

    """
    **Reset mixing ratios to their original values**
    """
    function restore_composition!(atmos::atmosphere.Atmos_t)
        for g in atmos.gas_names
            @turbo @. atmos.gas_vmr[g] = atmos.gas_ovmr[g]
        end
    end

    """
    **Calculate energy flux at each level.**

    Calculates flux components (radtrans, convection, etc.) and sums them to get total flux.
    Also updates thermodynamic properties (heat capacity, density, etc.) at each layer.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `latent::Bool`                    include condensation flux
    - `convect::Bool`                   include MLT convection flux
    - `sens_heat::Bool`                 include TKE sensible heat transport
    - `conduct::Bool`                   include conductive heat transport
    - `convect_sf::Float64`             scale factor applied to convection fluxes
    - `latent_sf::Float64`              scale factor applied to phase change fluxes
    - `calc_cf::Bool`                   calculate LW contribution function?
    - `rainout::Bool`                   allow rainout ( do not reset VMRs to dry values )
    """
    function calc_fluxes!(atmos::atmosphere.Atmos_t,
                          latent::Bool, convect::Bool, sens_heat::Bool, conduct::Bool;
                          convect_sf::Float64=1.0, latent_sf::Float64=1.0,
                          calc_cf::Bool=false, rainout::Bool=true)

        # Reset fluxes
        reset_fluxes!(atmos)

        # +Condensation and evaporation
        if atmos.condense_any && latent

            # Restore mixing ratios
            restore_composition!(atmos)
            atmosphere.calc_layer_props!(atmos)

            # Handle rainout
            chemistry.handle_saturation!(atmos)

            # Calculate latent heat flux
            condense_diffuse!(atmos)

            # Modulate?
            atmos.flux_l *= latent_sf

            # Add to total flux
            @turbo @. atmos.flux_tot += atmos.flux_l

            # Restore mixing ratios - do not allow rainout
            if !rainout
                restore_composition!(atmos)
            end
        end

        # Calculate layer properties
        atmosphere.calc_layer_props!(atmos)

        # +Radiation
        radtrans!(atmos, true, calc_cf=calc_cf)
        radtrans!(atmos, false)
        @turbo @. atmos.flux_tot += atmos.flux_n

        # +Dry convection
        if convect
            # Calc flux
            convection!(atmos)

            # Modulate?
            atmos.flux_cdry *= convect_sf

            # Add to total flux
            @turbo @. atmos.flux_tot += atmos.flux_cdry
        end

        # +Surface turbulence
        if sens_heat
            sensible!(atmos)
            atmos.flux_tot[end] += atmos.flux_sens
        end

        # +Conduction
        if conduct
            conduct!(atmos)
            @turbo @. atmos.flux_tot += atmos.flux_cdct
        end

        # Flux difference across each level
        # Positive value => heating
        atmos.flux_dif[1:end] .= (atmos.flux_tot[2:end] .- atmos.flux_tot[1:end-1])

        return nothing
    end

    """
    **Calculate heating rates at cell-centres from the total flux.**

    Requires the total flux to have already been set (atmos.flux_tot). Heating
    rates are calculated in units of kelvin per day.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    """
    function calc_hrates!(atmos::atmosphere.Atmos_t)

        dF::Float64 = 0.0
        dp::Float64 = 0.0

        fill!(atmos.heating_rate, 0.0)

        for i in 1:atmos.nlev_c
            dF = atmos.flux_tot[i+1] - atmos.flux_tot[i]
            dp = atmos.pl[i+1] - atmos.pl[i]
            atmos.heating_rate[i] = (atmos.layer_grav[i] / atmos.layer_cp[i]) * dF/dp # K/s
        end

        atmos.heating_rate *= 86400.0 # K/day

        return nothing
    end

end # end module

