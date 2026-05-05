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

    # Constants
    SKIP_SW_THRESH::Float64         = 1e-9      # skip SW calculation if TOA heating is below this threshold [W m-2]
    FILL_FINITE_FLUX::Float64       = 1.0       # filling value for NaN fluxes [W m-2]
    CONVECT_MIN_PRESSURE::Float64   = 1e-9      # lowest pressure at which convection is allowed [bar]
    CONVECT_REAL_GAS::Bool          = false     # use real gas EOS in convection scheme, if RG EOS enabled
    MIN_SKIN_D::Float64             = 1e-6      # minimum skin depth for conductive flux calculation [m]
    MAX_SKIN_D::Float64             = 1e6       # maximum skin depth for conductive flux calculation [m]
    ROUGHNESS_EPS::Float64          = 1e-3      # avoid blow-up of exchange coefficient when height ≈ roughness

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
    - `calc_cf::Bool`           also calculate contribution function?
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
            # Practical improved flux method (1985) with Elsasser's diffusivity (D=1.66)
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
            # Practical improved flux method (original form of 1980)
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
            # SOCRATES expects:
            #   w_cloud              -> Total cloud area fraction in layers, in [0, 1] (dimensionless)
            #   condensed_mix_ratio  -> Mass mixing ratios of condensate [kg kg-1], (LWC)
            #   condensed_dim_char   -> Characteristic dimensions of condensed species [m], (radius)
            # The LWC is the mass of condensate per mass of dry air
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

        return true
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
    function _radtrans_greygas!(atmos::atmosphere.Atmos_t)::Bool

        # Working layer transmissivity and emissivity
        trans::Float64 = 0.0

        # Down-directed SW and LW beams, looping from TOA downwards
        atmos.flux_d_sw[1] = atmos.toa_heating
        atmos.flux_d_lw[1] = 0.0
        for i in 1:atmos.nlev_c
            # Downward LW flux at bottom of layer
            trans = exp( (atmos.pl[i] - atmos.pl[i+1]) * atmos.κ_grey_lw / atmos.g[i] )
            atmos.flux_d_lw[i+1] = atmos.flux_d_lw[i] * trans + (phys.σSB * atmos.tmp[i]^4) * (1 - trans)

            # Downward SW flux at bottom of layer
            trans = exp( (atmos.pl[i] - atmos.pl[i+1]) * atmos.κ_grey_sw / atmos.g[i] )
            atmos.flux_d_sw[i+1] = atmos.flux_d_sw[i] * trans
        end

        # Up-directed LW beam, looping from surface upwards
        atmos.flux_u_lw[end] = phys.σSB * atmos.tmp_surf^4 * (1-atmos.albedo_s)
        for i in range(start=atmos.nlev_c, stop=1, step=-1)
            trans = exp( (atmos.pl[i] - atmos.pl[i+1]) * atmos.κ_grey_lw / atmos.g[i] )
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

        return true
    end

    """
    **Calculate radiative fluxes using the desired scheme.**

    Uses the configuration inside the atmos struct. Can either do LW or SW
    calculation, set by `lw` function argument.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `lw::Bool`                        longwave calculation? Else: shortwave
    - `calc_cf::Bool=false`             also calculate contribution function?

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

    """
    **Calculate turbulent kinetic energy (TKE) exchange coefficient**.

    Based on Monin–Obukhov similarity theory, from roughness length scale.
    See eq 9 in Nicholson & Benn (2006). Added small epsilon-factor to avoid function
    blowing-up around regime where height ≈ roughness.

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
        Kzz_min::Float64  = 0.0

        # Near-zero value
        Kzz_eps::Float64 = 1.0e-10

        # Find reference index for extension of Kzz, starting from convective regions
        i_Kzz_top::Int64 = atmos.nlev_l # default
        i_Kzz_bot::Int64 = atmos.nlev_l # default
        if any(atmos.Kzz .> Kzz_eps)
            # set to top of convective region
            i_Kzz_top = findfirst(x -> x > Kzz_eps, atmos.Kzz)
            # set to bottom of convective region
            i_Kzz_bot = findlast(x -> x > Kzz_eps, atmos.Kzz)
            # get minimum value, for filling intermediate zones
            Kzz_min = minimum(atmos.Kzz[atmos.Kzz .> Kzz_eps])
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

    """
    **Analytical diffusion scheme for condensation and evaporation energy.**

    Updates fluxes. Requires `chemistry.rainout_and_evaporate` to be called first.

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

            # add energy from this condesable to total energy from all condensables
            @. atmos.flux_l += atmos.phs_wrk_fl

            # calculate mask
            @. atmos.mask_l = (abs(atmos.flux_l) > 1.0e-30)

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
        for i in 1:atmos.nlev_c
            atmos.heating_rate[i] = (atmos.g[i] / atmos.layer_cp[i]) *
                                        atmos.flux_dif[i] / (atmos.pl[i+1] - atmos.pl[i])
        end

        atmos.heating_rate *= 86400.0 # K/day

        return true
    end

end # end module
