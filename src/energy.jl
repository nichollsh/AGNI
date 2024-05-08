# Contains the energy module, for everything relating to energy transport 

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 


module energy  

    # System libraries 
    using Revise
    using Printf
    using LinearAlgebra
    using Logging
    
    # Local files
    import atmosphere
    import moving_average
    import phys
    import spectrum

    """
    **Calculate radiative fluxes using SOCRATES.**

    Uses the configuration inside the atmos struct. Can either do LW or SW
    calculation, as required. Imports SOCRATES wrapper from the atmosphere
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
        if lw
            atmos.control.isolir = atmosphere.SOCRATES.rad_pcf.ip_infra_red
        else
            atmos.control.isolir = atmosphere.SOCRATES.rad_pcf.ip_solar
        end

        # Check files are acceptable and set instellation if doing SW flux
        if lw
            if !Bool(atmos.spectrum.Basic.l_present[6])
                error("The spectral file contains no data for the Planckian function.\nCheck that the file contains a stellar spectrum.")
            end

            if Bool(atmos.spectrum.Basic.l_present[2])
                atmos.control.l_solar_tail_flux = true
            end

        else
            if !Bool(atmos.spectrum.Basic.l_present[2])
                error("The spectral file contains no solar spectral data.")
            end
            
            atmos.bound.zen_0[1] = 1.0/cosd(atmos.zenith_degrees)   #   Convert the zenith angles to secants.
            atmos.bound.solar_irrad[1] = atmos.toa_heating / cosd(atmos.zenith_degrees)
        end

        atmos.bound.rho_alb[:, atmosphere.SOCRATES.rad_pcf.ip_surf_alb_diff, :] .= atmos.albedo_s

        # Set the two-stream approximation to be used (-t f)
        if lw
            atmos.control.i_2stream = atmosphere.SOCRATES.rad_pcf.ip_elsasser
            # Practical improved flux method (1985) with Elsasser's diffusivity (D=1.66)
        else
            atmos.control.i_2stream = atmosphere.SOCRATES.rad_pcf.ip_pifm80
            # Practical improved flux method (original form of 1980)
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
        # surface properties
        # see src/aux/assign_surface_char_cdf.f
        # IP_surface_char  = 51, file suffix 'surf'
        #####################################

        if atmos.control.i_angular_integration == atmosphere.SOCRATES.rad_pcf.ip_two_stream
            if !lw
                atmos.bound.rho_alb[:, atmosphere.SOCRATES.rad_pcf.ip_surf_alb_dir, :] .= atmos.bound.rho_alb[:, atmosphere.SOCRATES.rad_pcf.ip_surf_alb_diff, :]
            end
        end

        ###################################################
        # Cloud information
        ###################################################

        atmos.cld.w_cloud[1,:]               .= atmos.cloud_arr_f[:]   
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

        ######################################################
        # Run radiative transfer model
        ######################################################   

        # Calculate contribution function?
        atmos.control.l_contrib_func_band = calc_cf

        # Set composition for each gas,level
        for (i_gas,s_gas) in enumerate(atmos.gas_soc_names)
            for i in 1:atmos.nlev_c 
                # skip unspecified gases
                if s_gas in atmos.gas_all_names
                    # convert mole fraction to mass mixing ratio
                    atmos.atm.gas_mix_ratio[1, i, i_gas] = atmos.gas_all_dict[s_gas][i] * phys.lookup_safe("mmw", s_gas) / atmos.layer_mmw[i]
                else 
                    atmos.atm.gas_mix_ratio[1, i, i_gas] = 0.0
                end 
                # do not normalise MMRs to 1
            end
        end 

        # Do radiative transfer
        atmosphere.atmosphere.SOCRATES.radiance_calc(atmos.control, atmos.dimen, atmos.spectrum, 
                                atmos.atm, atmos.cld, atmos.aer, 
                                atmos.bound, atmos.radout)

        # Store new fluxes in atmos struct
        idx::Int = 1
        if lw 
            # LW case
            fill!(atmos.flux_u_lw, 0.0)
            fill!(atmos.flux_d_lw, 0.0)
            fill!(atmos.flux_n_lw, 0.0)
            fill!(atmos.band_u_lw, 0.0)
            fill!(atmos.band_d_lw, 0.0)
            fill!(atmos.band_n_lw, 0.0)
            for lv in 1:atmos.nlev_l      # sum over levels
                for ba in 1:atmos.nbands  # sum over bands
                    idx = lv+(ba-1)*atmos.nlev_l

                    atmos.band_d_lw[lv,ba] = atmos.radout.flux_down[idx]
                    atmos.band_u_lw[lv,ba] = atmos.radout.flux_up[idx]
                    atmos.band_n_lw[lv,ba] = atmos.band_u_lw[lv,ba] - atmos.band_d_lw[lv,ba]

                    atmos.flux_d_lw[lv] += max(0.0, atmos.radout.flux_down[idx])
                    atmos.flux_u_lw[lv] += max(0.0, atmos.radout.flux_up[idx])
                end 
                atmos.flux_n_lw[lv] = atmos.flux_u_lw[lv] - atmos.flux_d_lw[lv] 
            end

            # Normalised contribution function (only LW stream contributes)
            fill!(atmos.contfunc_norm,0.0)
            if calc_cf
                cf_max::Float64 = maximum(atmos.radout.contrib_funcf_band[1,:,:])
                for ba in 1:atmos.dimen.nd_channel
                    for lv in 1:atmos.nlev_c               
                        atmos.contfunc_norm[lv,ba] = atmos.radout.contrib_funcf_band[1,lv,ba]/cf_max
                    end 
                end
            end
            atmos.is_out_lw = true 
        else
            # SW case
            fill!(atmos.flux_u_sw, 0.0)
            fill!(atmos.flux_d_sw, 0.0)
            fill!(atmos.flux_n_sw, 0.0)
            fill!(atmos.band_u_sw, 0.0)
            fill!(atmos.band_d_sw, 0.0)
            fill!(atmos.band_n_sw, 0.0)
            for lv in 1:atmos.nlev_l                # sum over levels
                for ba in 1:atmos.dimen.nd_channel  # sum over bands
                    idx = lv+(ba-1)*atmos.nlev_l

                    atmos.band_d_sw[lv,ba] = atmos.radout.flux_down[idx]
                    atmos.band_u_sw[lv,ba] = atmos.radout.flux_up[idx]
                    atmos.band_n_sw[lv,ba] = atmos.band_u_sw[lv,ba] - atmos.band_d_sw[lv,ba]

                    atmos.flux_d_sw[lv] += max(0.0, atmos.radout.flux_down[idx])
                    atmos.flux_u_sw[lv] += max(0.0, atmos.radout.flux_up[idx])
                end 
                atmos.flux_n_sw[lv] = atmos.flux_u_sw[lv] - atmos.flux_d_sw[lv]
            end
            atmos.is_out_sw = true
        end

        # Store net fluxes when we have both SW and LW components
        if (atmos.is_out_lw && atmos.is_out_sw)
            atmos.flux_d = atmos.flux_d_lw .+ atmos.flux_d_sw
            atmos.flux_u = atmos.flux_u_lw .+ atmos.flux_u_sw
            atmos.flux_n = atmos.flux_n_lw .+ atmos.flux_n_sw
        end

        return nothing
    end # end of radtrans


    # Calculate sensible heat flux (turbulence at surface boundary)
    function sensible!(atmos::atmosphere.Atmos_t)
        # TKE scheme for this 1D case
        # transports energy from the surface to the bottom node
        atmos.flux_sens = atmos.layer_cp[end]*atmos.layer_mmw[end]*atmos.p[end]/(phys.R_gas*atmos.tmp[end]) * atmos.C_d * atmos.U * (atmos.tmp_surf-atmos.tmp[end])
        return nothing
    end

    # Calculate conductive fluxes 
    function conduct!(atmos::atmosphere.Atmos_t)
        # top layer 
        atmos.flux_cdct[1] = 0.0

        # bulk layers
        for i in 2:atmos.nlev_l-1  
            atmos.flux_cdct[i] = 0.5*(atmos.layer_kc[i-1]+atmos.layer_kc[i]) * (atmos.tmp[i]-atmos.tmp[i-1])/(atmos.z[i-1]-atmos.z[i])
        end 
        
        # bottom layer 
        atmos.flux_cdct[end] = atmos.layer_kc[end] * (atmos.tmp[end]-atmos.tmp_surf)/atmos.z[end]
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
    z<H) as per Blackadar (1962). Or alternatively it can be set equal to H.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `pmin::Float64=0.0`               pressure below which convection is disabled.
    - `mltype::Int=1`                   mixing length (0: fixed, 1: asymptotic)
    """
    function mlt!(atmos::atmosphere.Atmos_t; pmin::Float64=0.0, mltype::Int=1)

        # Reset arrays
        fill!(atmos.flux_cdry, 0.0)
        fill!(atmos.Kzz,    0.0)

        # Work variables 
        H::Float64 = 0.0; l::Float64 = 0.0; w::Float64 = 0.0
        m1::Float64 = 0.0; m2::Float64 = 0.0; mt::Float64 = 0.0
        grav::Float64 = 0.0; mu::Float64 = 0.0; c_p::Float64 = 0.0; rho::Float64 = 0.0
        grad_ad::Float64 = 0.0; grad_pr::Float64 = 0.0; grad_df::Float64 = 0.0

        # Loop from bottom upwards (over cell-edges)
        for i in range(start=atmos.nlev_l-1, step=-1, stop=3) 

            # Optionally skip low pressures 
            if atmos.pl[i] < pmin
                continue
            end

            # Skip condensing regions
            if (i <= atmos.nlev_c) && (atmos.mask_p[i] == atmos.mask_decay)
                continue
            end 

            m1 = atmos.layer_mass[i-1]
            m2 = atmos.layer_mass[i]
            mt = m1+m2

            grav = (atmos.layer_grav[i] * m2 + atmos.layer_grav[i-1] * m1)/mt
            mu   = (atmos.layer_mmw[i]  * m2 + atmos.layer_mmw[i-1]  * m1)/mt
            c_p  = (atmos.layer_cp[i]   * m2 + atmos.layer_cp[i-1]   * m1)/mt

            grad_ad = (phys.R_gas / mu) / c_p
            grad_pr = ( log(atmos.tmp[i-1]/atmos.tmp[i]) ) / ( log(atmos.p[i-1]/atmos.p[i]) )

            # Check instability
            if (grad_pr > grad_ad)

                # Check if this layer is condensing (this shouldn't ever be
                # true, because the condensation curve dT/dp is too shallow)
                if atmos.mask_p[i] > 0
                    @warn "Condensing region is unstable to dry convection"
                end 

                rho = (atmos.layer_density[i] * m2 + atmos.layer_density[i-1] * m1)/mt

                if i <= atmos.nlev_c-1
                    atmos.mask_c[i+1] = atmos.mask_decay
                end
                atmos.mask_c[i]   = atmos.mask_decay
                atmos.mask_c[i-1] = atmos.mask_decay
                
                # Pressure scale height
                H = phys.R_gas * atmos.tmpl[i] / (mu * grav)

                # Mixing length
                if mltype == 0
                    # Fixed
                    l = H
                elseif mltype == 1
                    # Asymptotic 
                    l = phys.k_vk * atmos.zl[i] / (1 + phys.k_vk * atmos.zl[i]/H)
                else 
                    # Otherwise
                    error("Invalid mixing length type selected")
                end

                # Characteristic velocity (from Brunt-Vasalla frequency of parcel oscillations)
                w = l * sqrt(grav/H * (grad_pr-grad_ad))

                # Dry convective flux
                atmos.flux_cdry[i] = 0.5 * rho * c_p * w * atmos.tmpl[i] * l/H * (grad_pr-grad_ad)

                # Thermal eddy diffusion coefficient
                atmos.Kzz[i] = w * l
            
            end
        end

        return nothing
    end # end of mlt



    """
    **Analytical relaxation scheme for condensation fluxes in saturated regions.**

    Calculates flux release by vertical latent heat transport using an heuristic 
    function of the temperature difference T-T_dew. 

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `condensates::Array`              list of condensates (strings)
    """
    function condense_relax!(atmos::atmosphere.Atmos_t, condensates::Array=[])

        # Reset flux
        fill!(atmos.flux_p, 0.0)

        # Check if empty
        if length(condensates) == 0
            return nothing 
        end 

        # Work variables
        a::Float64 = 2.0
        pp::Float64 = 0.0
        Psat::Float64 = 0.0
        dif::Array = zeros(Float64, atmos.nlev_c)

        # Calculate flux (negative) divergence due to latent heat release...
        # For all condensates
        for c in condensates
            # For all levels 
            for i in 1:atmos.nlev_c

                # Get partial pressure 
                pp = atmos.gas_all_dict[c][i] * atmos.p[i]
                if pp < 1.0e-10 
                    continue
                end

                # check saturation
                Psat = phys.calc_Psat(c, atmos.tmp[i])
                if (pp < Psat+1.0e-10)
                    continue 
                end 

                # relaxation function
                dif[i] += (1.0/a)*(pp-Psat)
                # dif[i] += a/Psat - a/pp

                # handle water cloud case 
                atmos.mask_p[i] = atmos.mask_decay 
                if c == "H2O"
                    atmos.cloud_arr_r[i] = atmos.cloud_val_r
                    atmos.cloud_arr_l[i] = atmos.cloud_val_l
                    atmos.cloud_arr_f[i] = atmos.cloud_val_f
                else
                    atmos.cloud_arr_r[i] = 0.0
                    atmos.cloud_arr_l[i] = 0.0
                    atmos.cloud_arr_f[i] = 0.0
                end 

            end # end levels 
        end # end condensates 

        # Convert divergence to cell-edge fluxes
        # Assuming zero condensation at surface
        for i in range(start=atmos.nlev_l-1, stop=1, step=-1)
            atmos.flux_p[i] = atmos.flux_p[i+1] - dif[i]
        end 

        return nothing
    end # end of condense_relax



    """
    **Calculate total flux at each level.**

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `condense::Bool`                  include condensation flux
    - `convect::Bool`                   include MLT convection flux
    - `sens_heat::Bool`                 include TKE sensible heat transport 
    - `conduct::Bool`                   include conductive heat transport
    - `condensates::Array=[]`           list of condensates included in relaxation scheme 
    - `convect_sf::Float64=1.0`         scale factor applied to convection fluxes
    """
    function calc_fluxes!(atmos::atmosphere.Atmos_t, 
                            condense::Bool, convect::Bool, sens_heat::Bool, conduct::Bool;
                            condensates::Array=[], convect_sf::Float64=1.0)

        # Reset fluxes
        fill!(atmos.flux_tot, 0.0)
        fill!(atmos.flux_dif, 0.0)

        # +Radiation
        energy.radtrans!(atmos, true)
        energy.radtrans!(atmos, false)
        atmos.flux_tot += atmos.flux_n

        # +Condensation
        if condense
            energy.condense_relax!(atmos, condensates)
            atmos.flux_tot += atmos.flux_p
        end

        # +Dry convection
        if convect
            # Calc flux
            energy.mlt!(atmos)

            # Stabilise?
            atmos.flux_cdry *= convect_sf

            # Add to total flux
            atmos.flux_tot += atmos.flux_cdry
        end

        # +Surface turbulence
        if sens_heat
            energy.sensible!(atmos)
            atmos.flux_tot[end] += atmos.flux_sens
        end

        # +Conduction 
        if conduct
            energy.conduct!(atmos)
            atmos.flux_tot += atmos.flux_cdct
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

        atmos.heating_rate[:] .= 0.0

        for i in 1:atmos.nlev_c
            dF = atmos.flux_tot[i+1] - atmos.flux_tot[i]
            dp = atmos.pl[i+1] - atmos.pl[i]
            atmos.heating_rate[i] = (atmos.layer_grav[i] / atmos.layer_cp[i]) * dF/dp # K/s
        end

        atmos.heating_rate *= 86400.0 # K/day

        return nothing
    end

    # Dry convective adjustment (returning the temperature tendency without modifying the atmos struct)
    function adjust_dry(atmos::atmosphere.Atmos_t, nsteps::Int)

        tmp_old = zeros(Float64, atmos.nlev_c)  # old temperatures
        tmp_tnd = zeros(Float64, atmos.nlev_c)  # temperature tendency

        tmp_old[:] .+= atmos.tmp[:]

        for i in 1:nsteps

            # Downward pass
            for i in 2:atmos.nlev_c

                T1 = atmos.tmp[i-1]    # upper layer
                p1 = atmos.p[i-1]

                T2 = atmos.tmp[i]  # lower layer
                p2 = atmos.p[i]
                
                cp = 0.5 * ( atmos.layer_cp[i-1] * atmos.layer_mmw[i-1] +  atmos.layer_cp[i] * atmos.layer_mmw[i])
                pfact = (p1/p2)^(phys.R_gas / cp)
                
                # If slope dT/dp is steeper than adiabat (unstable), adjust to adiabat
                if T1 < T2*pfact
                    Tbar = 0.5 * ( T1 + T2 )
                    T2 = 2.0 * Tbar / (1.0 + pfact)
                    T1 = T2 * pfact
                    atmos.tmp[i-1] = T1
                    atmos.tmp[i]   = T2
                end
            end

            # Upward pass
            for i in atmos.nlev_c:2

                T1 = atmos.tmp[i-1]
                p1 = atmos.p[i-1]

                T2 = atmos.tmp[i]
                p2 = atmos.p[i]
                
                cp = 0.5 * ( atmos.layer_cp[i-1] * atmos.layer_mmw[i-1] + atmos.layer_cp[i] * atmos.layer_mmw[i])
                pfact = (p1/p2)^(phys.R_gas / cp)

                if T1 < T2*pfact
                    Tbar = 0.5 * ( T1 + T2 )
                    T2 = 2.0 * Tbar / ( 1.0 + pfact)
                    T1 = T2 * pfact
                    atmos.tmp[i-1] = T1
                    atmos.tmp[i]   = T2 
                end 
            end
        end
        
        # Calculate tendency
        tmp_tnd[:] .= atmos.tmp[:] .- tmp_old[:]

        # Restore temperature array
        atmos.tmp[:] .= tmp_old[:]

        return tmp_tnd
    end

end # end module 
