# This file is part of AGNI. License is Apache-2.0: https://apache.org/licenses/LICENSE-2.0

"""
**Radiative transfer submodule of `energy`.**
"""
module radtrans

    # Local files
    import ..atmosphere
    import ..phys

    # Constants
    const SMALL_TRANS::Float64            = 1e-10     # minimum transmissivity
    const DFSVTY_FCTR::Float64            = 2.0       # diffusivity factor
    const SKIP_SW_THRESH::Float64         = 1e-9      # skip SW calculation if TOA heating is below this threshold [W m-2]
    const FILL_FINITE_FLUX::Float64       = 1.0       # filling value for NaN fluxes [W m-2]

    """
    **Set non-finite values in an array equal to a given fill value**.

    Arguments:
    - `arr`      array potentially containing non-finite values
    - `fill`     replacement value to fill with
    """
    function _make_finite!(arr, val)
        arr[findall(x -> !isfinite(x), arr)] .= val
    end

    """
    **Solve radiative transfer using SOCRATES**

    Imports SOCRATES wrapper from the atmosphere module, rather than loading it twice.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    - `lw::Bool`                True: longwave calculation. False: shortwave calculation.

    Optional arguments:
    - `calc_cf::Bool`           calculate contribution function and store optical depths?
    - `gauss_ir::Bool`          using gaussian angular integration in IR, otherwise uses two-stream approximation
    - `rescale_pf::Bool`        perform rescaling on phase function
    """
    function _radtrans_socrates!(atmos::atmosphere.Atmos_t, lw::Bool;
                                            calc_cf::Bool=false,
                                            gauss_ir::Bool=false,
                                            rescale_pf::Bool=false)::Bool


        # Longwave or shortwave calculation?
        if lw
            # Set source function
            atmos.control.isolir = atmosphere.SOCRATES.rad_pcf.ip_infra_red

            # Angular integration can be gauss or two-stream for LW
            if gauss_ir
                atmos.control.i_angular_integration = atmosphere.SOCRATES.rad_pcf.ip_ir_gauss
            else
                atmos.control.i_angular_integration = atmosphere.SOCRATES.rad_pcf.ip_two_stream
            end

            # Eddington's approximation
            # atmos.control.i_2stream = atmosphere.SOCRATES.rad_pcf.ip_eddington
            atmos.control.i_2stream = atmosphere.SOCRATES.rad_pcf.ip_elsasser

            # Check spectral file is ok
            if !Bool(atmos.spectrum.Basic.l_present[6])
                @warn("The spectral file contains no data for the Planck function. Check that the file contains a stellar spectrum.")
                return false
            end
            if Bool(atmos.spectrum.Basic.l_present[2])
                atmos.control.l_solar_tail_flux = true
            end
        else
            # Set source function
            atmos.control.isolir = atmosphere.SOCRATES.rad_pcf.ip_solar

            # Angular integration is always two-stream for SW
            atmos.control.i_angular_integration = atmosphere.SOCRATES.rad_pcf.ip_two_stream

            # Eddington's approximation
            # atmos.control.i_2stream = atmosphere.SOCRATES.rad_pcf.ip_eddington
            atmos.control.i_2stream = atmosphere.SOCRATES.rad_pcf.ip_pifm80

            # SOCRATES requires this to be passed as two variables, since it
            #       needs to know the angle of the direct beam.
            #   - Convert the zenith angles to secants.
            atmos.bound.zen_0[1] = 1.0/cosd(atmos.zenith_degrees)
            #   - Pass effective solar constant
            atmos.bound.solar_irrad[1] = atmos.instellation *
                                            (1.0 - atmos.albedo_b) * atmos.s0_fact

            # Check spectral file is ok
            if !Bool(atmos.spectrum.Basic.l_present[2])
                @warn ("The spectral file contains no solar spectral data.")
                return false
            end
        end

        #####################################
        # Angular integration
        # see src/aux/angular_control_cdf.f
        #####################################

        # Cl_run_cdf +R flag
        atmos.control.l_rescale = rescale_pf
        atmos.control.l_henyey_greenstein_pf = rescale_pf

        # The internal SOCRATES solver used for the two-stream calculations (-v flag)
        if atmos.control.l_cloud
            # 16 is recommended for cloudy-sky (ip_solver_mix_direct_hogan)
            atmos.control.i_solver = atmosphere.SOCRATES.rad_pcf.ip_solver_mix_direct_hogan

            # 17 is recommended for cloud with separate stratiform and convective regions
            # atmos.control.i_solver = atmosphere.SOCRATES.rad_pcf.ip_solver_triple_hogan
        else
            # 13 is recommended for clear-sky (Direct solution in a homogeneous column)
            atmos.control.i_solver = atmosphere.SOCRATES.rad_pcf.ip_solver_homogen_direct

            # Direct mixed column scheme for full fluxes
            # atmos.control.i_solver = atmosphere.SOCRATES.rad_pcf.ip_solver_mix_direct
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
            # SOCRATES expects:
            #   w_cloud              -> Total cloud area fraction in layers, in [0, 1] (dimensionless)
            #   condensed_mix_ratio  -> Mass mixing ratios of condensate [kg kg-1], (LWC)
            #   condensed_dim_char   -> Characteristic dimensions of condensed species [m], (radius)
            # The LWC is the mass of condensate per mass of air
            #   Refer to opt_propt_water_cloud.f90 (L234); c.f. Slingo & Schrecker (1982) Eq 15,16,17
            atmos.cld.w_cloud[1,:]               .= atmos.cloud_arr_f[:]
            atmos.cld.condensed_mix_ratio[1,:,1] .= atmos.cloud_arr_l[:]
            atmos.cld.condensed_dim_char[1,:,1]  .= atmos.cloud_arr_r[:]
        end

        ###################################################
        # Aerosol information
        ###################################################

        # Set mixing ratio profiles for aerosols
        fill!(atmos.aer.mix_ratio, 0.0)
        if atmos.control.l_aerosol
            for i = 1:atmos.spectrum.Aerosol.n_aerosol_mr
                atmos.aer.mix_ratio[1, :, i] .= atmos.aerosol_arr_l[atmos.aerosol_names[i]][:]
            end
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

        atmos.atm.mass[1, :]        .= atmos.layer_σ[:]
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
        if !lw && (atmos.toa_heating < SKIP_SW_THRESH)
            # If no stellar flux is reaching the atmosphere, skip the SW calculation
            fill!(atmos.radout.flux_down, 0.0)
            fill!(atmos.radout.flux_up, 0.0)
        else
            atmosphere.SOCRATES.radiance_calc(atmos.control,
                                                     atmos.dimen, atmos.spectrum,
                                                     atmos.atm, atmos.cld, atmos.aer,
                                                     atmos.bound, atmos.radout)
        end

        # Check finite
        if !all(isfinite, atmos.radout.flux_down)
            if lw
                @warn "Non-finite value in LW DN flux array"
            else
                @warn "Non-finite value in SW DN flux array"
            end
            _make_finite!(atmos.radout.flux_down, FILL_FINITE_FLUX)
        end
        if !all(isfinite, atmos.radout.flux_up)
            if lw
                @warn "Non-finite value in LW UP flux array"
            else
                @warn "Non-finite value in SW UP flux array"
            end
            _make_finite!(atmos.radout.flux_up, FILL_FINITE_FLUX)
        end

        # Store new fluxes in atmos struct
        idx::Int64 = 1
        trans_prev::Float64 = 1.0
        dlog10p::Float64 = 0.0
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

            # Contribution function
            fill!(atmos.contfunc_band ,0.0)
            fill!(atmos.tau_band, 0.0)
            if calc_cf
                # get contribution function
                for ba in 1:atmos.dimen.nd_channel
                    for lv in 1:atmos.nlev_c
                        atmos.contfunc_band[lv,ba] = atmos.radout.contrib_funcf_band[1,lv,ba]
                    end
                end

                # get optical depth from contribution function
                for ba in 1:atmos.dimen.nd_channel

                    # calculate source function in this band
                    atmosphere.SOCRATES.diff_planck_source(
                        atmos.control, atmos.dimen, atmos.spectrum,
                        atmos.atm, atmos.bound, ba, atmos.planck
                    )

                    # loop through layers from TOA downwards
                    trans_prev = 1.0 # reset transmissivity at TOA, to unity
                    for lv in 1:atmos.nlev_c
                        # pressure change across level
                        dlog10p = log10(atmos.atm.p_level[1, lv]) -
                                        log10(atmos.atm.p_level[1, lv-1])

                        # transmissivity change from cff
                        # dT = cff * dffsvty * d(log pressure) / (2 * source function)
                        delta_trans = atmos.contfunc_band[lv, ba] *
                                        DFSVTY_FCTR * dlog10p /
                                        (2.0 * max(atmos.planck.flux[1, lv], SMALL_TRANS))

                        # decrease transmissivity
                        trans_now = clamp(trans_prev - delta_trans, SMALL_TRANS, 1.0)
                        trans_prev = trans_now

                        # set optical depth at bottom of this layer from transmissivity
                        #   trans = exp(-tau * diffusivity_factor)
                        atmos.tau_band[lv+1, ba] = -log(trans_now) / DFSVTY_FCTR

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

        return true
    end

    """
    **Solve RT using double grey-gas formulation**

    Simple two-stream double grey RT solver which integrates fluxes from the TOA and BOA.
    Uses two opacity values to represent the LW and SW components of the flux field.
    Loosely following this tutorial, which is based on Pierrehumbert (2010).
    https://brian-rose.github.io/ClimateLaboratoryBook/courseware/radiative-transfer/

    * Optical depth across each layer, `τ = κ * (p_bot - p_top) / g`
    * Transmissivity of each layer, `T = exp(-τ)`
    * Emissivity of each layer, `ε = 1 - T`
    * Where κ is the opacity, p is pressure, and g is gravity.


    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    """
    function _radtrans_greygas!(atmos::atmosphere.Atmos_t)::Bool

        # Working layer's transmissivity and optical depth
        trans::Float64 = 0.0  # reused for LW and SW
        tau_lw::Float64 = 0.0
        tau_sw::Float64 = 0.0

        # Set TOA boundary conditions
        atmos.flux_d_sw[1] = atmos.toa_heating
        atmos.flux_d_lw[1] = 0.0
        fill!(atmos.tau_band, 0.0)

        # Down-directed SW and LW beams, looping from TOA downwards
        for i in 1:atmos.nlev_c
            # Downward LW flux at bottom of layer
            tau_lw = (atmos.pl[i+1] - atmos.pl[i]) * atmos.κ_grey_lw / atmos.g[i]
            trans = exp( -tau_lw )
            atmos.flux_d_lw[i+1] = atmos.flux_d_lw[i] * trans + (phys.σSB * atmos.tmp[i]^4) * (1 - trans)

            # Downward SW flux at bottom of layer
            tau_sw = (atmos.pl[i+1] - atmos.pl[i]) * atmos.κ_grey_sw / atmos.g[i]
            trans = exp( -tau_sw )
            atmos.flux_d_sw[i+1] = atmos.flux_d_sw[i] * trans

            # Store LW+SW optical depth across this layer
            atmos.tau_band[i+1, 1] = atmos.tau_band[i, 1] + (tau_lw + tau_sw)
        end

        # Up-directed LW beam, looping from surface upwards
        atmos.flux_u_lw[end] = phys.σSB * atmos.tmp_surf^4 * (1-atmos.albedo_s)
        for i in range(start=atmos.nlev_c, stop=1, step=-1)
            tau_lw = (atmos.pl[i+1] - atmos.pl[i]) * atmos.κ_grey_lw / atmos.g[i]
            trans = exp( -tau_lw )
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

        # Set band arrays (just one band for greygas)
        for i in 1:atmos.nlev_l
            atmos.band_d_lw[i,1] = atmos.flux_d_lw[i]
            atmos.band_u_lw[i,1] = atmos.flux_u_lw[i]
            atmos.band_d_sw[i,1] = atmos.flux_d_sw[i]
            atmos.band_u_sw[i,1] = atmos.flux_u_sw[i]
        end

        return true
    end

    """
    **Calculate radiative fluxes using the desired scheme.**

    Uses the configuration inside the atmos struct. Can either do LW or SW
    calculation, set by `lw` function argument.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `lw::Bool`                        longwave calculation? Else: shortwave
    - `calc_cf::Bool=false`             calculate contribution function and optical depths?

    Returns:
    - `Bool`                            whether the calculation succeeded
    """
    function radtrans!(atmos::atmosphere.Atmos_t, lw::Bool; calc_cf::Bool=false)::Bool
        if !atmos.is_alloc
            @warn "Atmosphere arrays have not been allocated"
            return false
        end
        if !atmos.is_param
            @warn "Atmosphere parameters have not been set"
            return false
        end

        atmos.num_rt_eval += 1

        if atmos.benchmark
            time_start::UInt64 = time_ns()
        end

        # Downward SW flux in atmosphere at TOA
        # atmos.toa_heating = atmosphere.calc_toa_heating(atmos)

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
            return false
        end

        # Store time
        if atmos.benchmark
            atmos.tim_rt_eval += time_ns() - time_start
        end

        return true
    end # end of radtrans

    export radtrans!

end
