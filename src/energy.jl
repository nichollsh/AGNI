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

    # Local files
    import ..atmosphere
    import ..phys
    import ..chemistry
    import ..spectrum

    const FILL_FINITE_FLUX::Float64 = 1.0

    """
    **Set non-finite values in an array equal to a given fill value**.

    Arguments:
    - `arr`      array potentially containing non-finite values
    - `fill`     replacement value to fill with
    """
    function make_finite!(arr, val)
        arr[findall(x -> !isfinite(x), arr)] .= val
    end

    """
    **Solve radiative transfer using SOCRATES**

    Imports SOCRATES wrapper from the atmosphere module, rather than loading it twice.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `lw::Bool`                        longwave calculation? Else: shortwave
    - `calc_cf::Bool=false`             also calculate contribution function?
    """
    function _radtrans_socrates!(atmos::atmosphere.Atmos_t, lw::Bool; calc_cf::Bool=false)

        # Longwave or shortwave calculation?
        # Set the two-stream approximation to be used (-t f)
        if lw
            atmos.control.isolir = atmosphere.SOCRATES.rad_pcf.ip_infra_red

            # Eddington's approximation
            # atmos.control.i_2stream = atmosphere.SOCRATES.rad_pcf.ip_eddington

            # Practical improved flux method (1985) with Elsasser's diffusivity (D=1.66)
            atmos.control.i_2stream = atmosphere.SOCRATES.rad_pcf.ip_elsasser
        else
            atmos.control.isolir = atmosphere.SOCRATES.rad_pcf.ip_solar

            # Eddington's approximation
            # atmos.control.i_2stream = atmosphere.SOCRATES.rad_pcf.ip_eddington

            # Practical improved flux method (original form of 1980)
            atmos.control.i_2stream = atmosphere.SOCRATES.rad_pcf.ip_pifm80
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
            atmos.control.i_solver = atmosphere.SOCRATES.rad_pcf.ip_solver_mix_direct_hogan

            # 17 is recommended for cloud with separate stratiform and convective regions
            # atmos.control.i_solver = atmosphere.SOCRATES.rad_pcf.ip_solver_triple_hogan
        else
            # 13 is recommended for clear-sky (Direct solution in a homogeneous column)
            atmos.control.i_solver = atmosphere.SOCRATES.rad_pcf.ip_solver_homogen_direct

            # 1 is also possible (Pentadiagonal solver for homogeneous column)
            # atmos.control.i_solver = atmosphere.SOCRATES.rad_pcf.ip_solver_pentadiagonal
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

        if atmos.control.l_cloud
            atmos.cld.w_cloud[1,:]               .= atmos.cloud_arr_f[:]    # Total cloud area fraction in layers
            atmos.cld.condensed_mix_ratio[1,:,1] .= atmos.cloud_arr_l[:]    # Mass mixing ratios of condensate
            atmos.cld.condensed_dim_char[1,:,1]  .= atmos.cloud_arr_r[:]    # Characteristic dimensions of condensed species
        else
            atmos.cld.w_cloud[1,:]               .= 0.0
            atmos.cld.condensed_mix_ratio[1,:,1] .= 0.0
            atmos.cld.condensed_dim_char[1,:,1]  .= 0.0
        end

        ###################################################
        # Treatment of scattering
        ###################################################

        atmos.control.i_scatter_method = atmosphere.SOCRATES.rad_pcf.ip_scatter_full
        for i in atmos.control.first_band:atmos.control.last_band
            atmos.control.i_scatter_method_band[i] = atmos.control.i_scatter_method
        end

        ####################################################
        # Temperature, pressure, radius, etc.
        ###################################################

        atmos.atm.p[1, :]           .= atmos.p[:]
        atmos.atm.r_layer[1,:]      .= atmos.r[:]
        atmos.atm.t[1, :]           .= atmos.tmp[:]

        atmos.atm.p_level[1, 0:end] .= atmos.pl[:]
        atmos.atm.r_level[1, 0:end] .= atmos.rl[:]
        atmos.atm.t_level[1, 0:end] .= atmos.tmpl[:]

        atmos.atm.mass[1, :]        .= atmos.layer_mass[:]
        atmos.atm.density[1,:]      .= atmos.layer_ρ[:]

        if lw
            atmos.bound.t_ground[1] = atmos.tmp_surf
        end

        if lw
            atmos.control.l_ir_source_quad = true
        end

        ####################################################
        # Pass surface flux to SOCRATES
        ###################################################

        # Pass to socrates array
        #     I would argue that the 1-albedo term shouldn't be here, but it is to correct
        #     for it also (strangely) appearing inside diff_planck_source_mod.f90 on
        #     line 129. Having this 1-albedo term (and using this low-order integration)
        #     gives the correct results from my tests versus SOCRATES's native function.
        @inbounds for i in 1:atmos.nbands
            atmos.bound.flux_ground[1,i] = atmos.surf_flux[i] * atmos.surf_e_arr[i]
        end

        ######################################################
        # Run SOCRATES radiative transfer calculation
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

        # Ensure all VMRs are between 0 and 1
        clamp!(atmos.atm.gas_mix_ratio, 0.0, 1.0)

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
            make_finite!(atmos.radout.flux_down, FILL_FINITE_FLUX)
        end
        if !all(isfinite, atmos.radout.flux_up)
            if lw
                @error "Non-finite value in LW UP flux array"
            else
                @error "Non-finite value in SW UP flux array"
            end
            make_finite!(atmos.radout.flux_up, FILL_FINITE_FLUX)
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
    end

    """
    **Solve RT using double grey-gas formulation**

    Simple two-stream double grey RT solver which integrates fluxes from the TOA and BOA.

    Uses two opacity values to represent the LW and SW components of the flux field.

    Loosely following this tutorial, which is based on Pierrehumbert (2010).
    https://brian-rose.github.io/ClimateLaboratoryBook/courseware/radiative-transfer/

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    """
    function _radtrans_greygas!(atmos)

        # Working layer transmissivity and emissivity
        trans::Float64 = 0.0

        # Down-directed SW and LW beams, looping from TOA downwards
        atmos.flux_d_sw[1] = atmos.toa_heating
        atmos.flux_d_lw[1] = 0.0
        for i in 1:atmos.nlev_c
            # Downward LW flux at bottom of layer
            trans = exp( (atmos.pl[i] - atmos.pl[i+1]) * atmos.κ_grey_lw / atmos.layer_grav[i] )
            atmos.flux_d_lw[i+1] = atmos.flux_d_lw[i] * trans + (phys.σSB * atmos.tmp[i]^4) * (1 - trans)

            # Downward SW flux at bottom of layer
            trans = exp( (atmos.pl[i] - atmos.pl[i+1]) * atmos.κ_grey_sw / atmos.layer_grav[i] )
            atmos.flux_d_sw[i+1] = atmos.flux_d_sw[i] * trans
        end

        # Up-directed LW beam, looping from surface upwards
        atmos.flux_u_lw[end] = phys.σSB * atmos.tmp_surf^4 * (1-atmos.albedo_s)
        for i in range(start=atmos.nlev_c, stop=1, step=-1)
            trans = exp( (atmos.pl[i] - atmos.pl[i+1]) * atmos.κ_grey_lw / atmos.layer_grav[i] )
            atmos.flux_u_lw[i] = atmos.flux_u_lw[i+1] * trans + (phys.σSB * atmos.tmp[i]^4) * (1 - trans)
        end

        # Set other arrays to zero
        fill!(atmos.flux_u_sw, 0.0)
        atmos.is_out_sw = true
        atmos.is_out_lw = true

        # Set net arrays
        atmos.flux_d    = atmos.flux_d_lw + atmos.flux_d_sw  # net down
        atmos.flux_u    = atmos.flux_u_lw + atmos.flux_u_sw  # net up
        atmos.flux_n_lw = atmos.flux_u_lw - atmos.flux_d_lw  # net lw
        atmos.flux_n_sw = atmos.flux_u_sw - atmos.flux_d_sw  # net sw
        atmos.flux_n    = atmos.flux_n_lw + atmos.flux_n_sw  # net
    end

    """
    **Calculate radiative fluxes using the desired scheme.**

    Uses the configuration inside the atmos struct. Can either do LW or SW
    calculation, set by `lw` function argument.

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

        if atmos.benchmark
            time_start::UInt64 = time_ns()
        end

        # Downward SW flux in atmosphere at TOA
        atmos.toa_heating = atmos.instellation * (1.0 - atmos.albedo_b) *
                                    atmos.s0_fact * cosd(atmos.zenith_degrees)

        # Set flux in surface emission, by band
        #     Equal to integral of planck function over band width, which in
        #     this case is done by simply evaluating at the midpoint and
        #     multiplying by band width. Scaled by the emissivity.
        @. atmos.surf_flux = phys.evaluate_planck(atmos.bands_cen, atmos.tmp_surf) *
                                atmos.bands_wid * 1e9 * atmos.surf_e_arr


        # Run the RT using the desired scheme
        if atmos.rt_scheme == atmosphere.RT_SOCRATES
            _radtrans_socrates!(atmos, lw, calc_cf=calc_cf)

        elseif atmos.rt_scheme == atmosphere.RT_GREYGAS
            _radtrans_greygas!(atmos)

        else
            @error "Invalid RT scheme: $(atmos.rt_scheme)"
        end

        # Store time
        if atmos.benchmark
            atmos.tim_rt_eval += time_ns() - time_start
        end

        return nothing
    end # end of radtrans


    # Calculate sensible heat flux (turbulence at surface boundary)
    function sensible!(atmos::atmosphere.Atmos_t)

        # Calculate exchange coefficient
        #    Based on Monin–Obukhov similarity theory, from roughness length scale.
        #    See eq 9 in Nicholson & Benn (2009)
        #    Added small epsilon-factor to ensure that this does not blow-up
        atmos.C_d = phys.k_vk^2 / log(max(atmos.r[end]-atmos.rp, 1e-30)/atmos.surf_roughness)


        # TKE scheme for this 1D case
        # transports energy from the surface to the bottom node
        atmos.flux_sens = atmos.layer_cp[end]*atmos.layer_μ[end]*
                            atmos.p[end]/(phys.R_gas*atmos.tmp[end]) *
                            atmos.C_d * atmos.surf_windspeed *
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
                                                        atmos.layer_thick[i]
        end

        # bottom layer
        atmos.flux_cdct[end] = atmos.layer_kc[end] * (atmos.tmp[end]-atmos.tmp_surf) /
                                                      (atmos.r[end] - atmos.rp)
        return nothing
    end


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
    - `pmin::Float64`           pressure [bar] below which convection is disabled
    """
    function convection!(atmos::atmosphere.Atmos_t; pmin::Float64=1.0e-4)

        # Reset arrays
        fill!(atmos.mask_c,     false)
        fill!(atmos.flux_cdry,  0.0)
        fill!(atmos.Kzz,        0.0)
        fill!(atmos.λ_conv,     0.0)
        fill!(atmos.w_conv,     0.0)

        # Work variables
        Hp::Float64 = 0.0; hgt::Float64 = 0.0
        m1::Float64 = 0.0; m2::Float64 = 0.0; mt::Float64 = 0.0
        grav::Float64 = 0.0; mu::Float64 = 0.0; c_p::Float64 = 0.0; rho::Float64 = 0.0
        ∇_ad::Float64 = 0.0; ∇_pr::Float64 = 0.0; ∇_μ::Float64 = 0.0; staby::Float64 = 0.0


        # Loop from bottom upwards (over cell-edges)
        for i in range(start=atmos.nlev_l-1, step=-1, stop=2)

            # Optionally skip low pressures
            if atmos.pl[i] <= pmin * 1.0e5  # convert bar to Pa
                break
            end

            # Profile lapse rate: d(ln T)/d(ln P) = (P/T)*(dT/dP)
            ∇_pr = log(atmos.tmp[i-1]/atmos.tmp[i]) / log(atmos.p[i-1]/atmos.p[i])

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
                atmos.w_conv[i] = atmos.λ_conv[i] * sqrt(grav/Hp * staby)

                # Dry convective flux
                atmos.flux_cdry[i] = 0.5*rho*c_p*atmos.w_conv[i] * atmos.tmpl[i] * (atmos.λ_conv[i]/Hp) * staby

                # Convection eddy diffusion coefficient [m2 s-1]
                atmos.Kzz[i] =  atmos.w_conv[i] * atmos.λ_conv[i]

            end
        end

        # Set surface quantities
        atmos.Kzz[end]    = atmos.Kzz[end-1]
        atmos.w_conv[end] = atmos.w_conv[end-1]
        atmos.λ_conv[end] = atmos.λ_conv[end-1]

        # Check for spurious shallow convection occuring ABOVE condensing regions
        #    If found, reset convective flux to zero AT THIS LAYER ONLY.
        #    This is okay because this shouldn't physically happen, and will only occur
        #    because of weird numerical issues which only act to make solving difficult.
        # @inbounds for i in 1:atmos.nlev_l-1
        #     if (!atmos.mask_l[i] && any(atmos.mask_l[i+1:end])) #|| (atmos.mask_l[i] && !atmos.mask_c[i-1] && !atmos.mask_c[i+1])
        #         atmos.mask_c[i] = false
        #         atmos.flux_cdry[i] = 0.0
        #     end
        # end

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
        @inbounds for i in 1:atmos.nlev_l
            if atmos.flux_cdry[i] > Fthresh
                i_cnvct_bot = i
            end
        end

        # Find top of convective region (looping upwards from bottom)
        i_cnvct_top::Int = i_cnvct_bot
        atmos.Kzz_pbreak = min(1e5, atmos.pl[end])
        atmos.Kzz_kbreak = atmos.Kzz_floor
        for i in range(start=i_cnvct_bot, step=-1, stop=1)
            if atmos.flux_cdry[i] > Fthresh
                i_cnvct_top = i
            end
        end

        # Store breakpoint values
        atmos.Kzz_pbreak = atmos.pl[i_cnvct_top]
        atmos.Kzz_kbreak = atmos.Kzz[i_cnvct_top]

        # Set Kzz within and below convective region to constant value. This value best
        #    represents the diffusive mixing in this region.
        atmos.Kzz[i_cnvct_bot:end] .= atmos.Kzz[i_cnvct_bot]

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

    Updates fluxes. Requires `chemistry.handle_saturation` to be called first.

    Integrates from bottom of model upwards. Based on the amount of
    phase change at each level, a phase change flux is calculated by assuming
    a fixed condensation timescale.

    If evaporation is enabled, then integrates from top downwards to determine flux from
    re-evaporation of droplets

    Any droplets which reach the ground go towards forming an ocean.

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
                atmos.phs_wrk_df[i] += phys.get_Lv(atmos.gas_dat[c], atmos.tmp[i]) *
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

            # add energy from this gas to total
            @. atmos.flux_l += atmos.phs_wrk_fl

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
            @. atmos.gas_vmr[g] = atmos.gas_ovmr[g]
        end
    end

    """
    **Calculate energy flux at each level.**

    Calculates flux components (radtrans, convection, etc.) and sums them to get total flux.
    Also updates thermodynamic properties (heat capacity, density, etc.) at each layer.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `radiative::Bool`                 include radiation fluxes
    - `latent::Bool`                    include condensation flux
    - `convect::Bool`                   include MLT convection flux
    - `sens_heat::Bool`                 include TKE sensible heat transport
    - `conduct::Bool`                   include conductive heat transport
    - `convect_sf::Float64`             scale factor applied to convection fluxes
    - `latent_sf::Float64`              scale factor applied to phase change fluxes
    - `calc_cf::Bool`                   calculate LW contribution function?
    - `rainout::Bool`                   allow rainout ( do not reset VMRs to dry values )

    Returns:
    - `ok::Bool`                        calculation performed ok?
    """
    function calc_fluxes!(atmos::atmosphere.Atmos_t,
                          radiative::Bool,
                          latent::Bool, convect::Bool, sens_heat::Bool, conduct::Bool;
                          convect_sf::Float64=1.0, latent_sf::Float64=1.0,
                          calc_cf::Bool=false, rainout::Bool=true)::Bool

        # Reset fluxes
        reset_fluxes!(atmos)

        ok::Bool = true

        # +Condensation and evaporation
        if atmos.condense_any && (latent || rainout)

            # Restore mixing ratios
            restore_composition!(atmos)
            ok &= atmosphere.calc_layer_props!(atmos)

            # Handle rainout
            chemistry.handle_saturation!(atmos)

            # Calculate latent heat flux
            if latent
                condense_diffuse!(atmos)                    # Calculate latent heat flux
                atmos.flux_l *= latent_sf                   # Modulate for stability?
                @. atmos.flux_tot += atmos.flux_l    # Add to total flux
            end

            # Restore mixing ratios - do not allow rainout
            if !rainout
                restore_composition!(atmos)
            end
        end

        # Recalculate layer properties
        #    Returns false if atmosphere becomes non-hydrostatic
        ok &= atmosphere.calc_layer_props!(atmos)

        # +Radiation
        if radiative
            radtrans!(atmos, true, calc_cf=calc_cf)   # Longwave
            radtrans!(atmos, false)                   # Shortwave
            @. atmos.flux_tot += atmos.flux_n  # Add to total flux
        end

        # +Dry convection
        if convect
            convection!(atmos)                          # Calc dry convection heat flux
            atmos.flux_cdry *= convect_sf               # Modulate for stability?
            @. atmos.flux_tot += atmos.flux_cdry # Add to total flux
        end

        # +Surface turbulence
        if sens_heat
            sensible!(atmos)
            atmos.flux_tot[end] += atmos.flux_sens
        end

        # +Conduction
        if conduct
            conduct!(atmos)
            @. atmos.flux_tot += atmos.flux_cdct
        end

        # Flux difference across each level
        # Positive value => heating
        atmos.flux_dif[1:end] .= (atmos.flux_tot[2:end] .- atmos.flux_tot[1:end-1])

        return ok
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

