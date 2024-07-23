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
    import ..spectrum

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
        # Surface albedo
        #####################################

        fill!(atmos.bound.rho_alb, 0.0)
        atmos.bound.rho_alb[1, atmosphere.SOCRATES.rad_pcf.ip_surf_alb_diff, :] .= atmos.albedo_s_arr
        atmos.bound.rho_alb[1, atmosphere.SOCRATES.rad_pcf.ip_surf_alb_dir,  :] .= atmos.albedo_s_arr

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

        # set surface emission
        pl_w::Float64 = 0.0; pl_x::Float64 = 0.0
        for i in 1:atmos.nbands
            # get band width, midpoint 
            pl_w = 1e9 * (atmos.bands_max[i] - atmos.bands_min[i])
            pl_x = 1e9 * atmos.bands_min[i] + 0.5 * pl_w 

            # Set flux in band 
            #  Equal to integral of planck function over band width, which in 
            #  this case is done by simply evaluating at the midpoint and 
            #  multiplying by band width.
            #  It is then scaled by the emissivity and 1-albedo. I'd argue
            #  that the albedo term shouldn't be here, but it's to correct for 
            #  it also (strangely) appearing inside diff_planck_source_mod.f90
            #  on line 129. 
            #  Having this 1-albedo term and using this low-order integration 
            #  gives the correct results from my tests.  
            atmos.bound.flux_ground[1,i] = phys.evaluate_planck(pl_x, atmos.tmp_surf)*pl_w*
                                                atmos.emiss_s_arr[i] *
                                                (1.0 - atmos.albedo_s_arr[i])
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
                                                            atmos.layer_mmw[i]
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
            fill!(atmos.flux_u_lw, 0.0)
            fill!(atmos.flux_d_lw, 0.0)
            fill!(atmos.flux_n_lw, 0.0)
            fill!(atmos.band_u_lw, 0.0)
            fill!(atmos.band_d_lw, 0.0)
            fill!(atmos.band_n_lw, 0.0)
            for lv in 1:atmos.nlev_l      # sum over levels
                for ba in 1:atmos.dimen.nd_channel  # sum over bands
                    idx = lv+(ba-1)*atmos.nlev_l

                    atmos.band_d_lw[lv,ba] = max(0.0, atmos.radout.flux_down[idx])
                    atmos.band_u_lw[lv,ba] = max(0.0, atmos.radout.flux_up[idx])
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
                        atmos.contfunc_norm[lv,ba] = 
                                            atmos.radout.contrib_funcf_band[1,lv,ba]/cf_max
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

                    atmos.band_d_sw[lv,ba] = max(0.0,atmos.radout.flux_down[idx])
                    atmos.band_u_sw[lv,ba] = max(0.0,atmos.radout.flux_up[idx])
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
        atmos.flux_sens = atmos.layer_cp[end]*atmos.layer_mmw[end]*
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
        for i in 2:atmos.nlev_l-1  
            atmos.flux_cdct[i] = 0.5*(atmos.layer_kc[i-1]+atmos.layer_kc[i]) * 
                                     (atmos.tmp[i]-atmos.tmp[i-1])/(atmos.z[i-1]-atmos.z[i])
        end 
        
        # bottom layer 
        atmos.flux_cdct[end] = atmos.layer_kc[end] *
                                    (atmos.tmp[end]-atmos.tmp_surf) / atmos.z[end]
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

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    - `pmin::Float64`           pressure [bar] below which convection is disabled
    - `mltype::Int`             mixing length value (0: scale height, 1: asymptotic)
    """
    function mlt!(atmos::atmosphere.Atmos_t; pmin::Float64=1.0e-9, mltype::Int=1)

        pmin *= 1.0e5 # convert bar to Pa

        # Reset arrays
        fill!(atmos.flux_cdry, 0.0)
        fill!(atmos.Kzz,    0.0)

        # Work variables 
        H::Float64 = 0.0; l::Float64 = 0.0; w::Float64 = 0.0
        m1::Float64 = 0.0; m2::Float64 = 0.0; mt::Float64 = 0.0
        grav::Float64 = 0.0; mu::Float64 = 0.0; c_p::Float64 = 0.0; rho::Float64 = 0.0
        grad_ad::Float64 = 0.0; grad_pr::Float64 = 0.0; grad_df::Float64 = 0.0
        beta::Float64 = 0.0; xv::Float64=0.0; xv_av::Float64=0.0
        inhib::Float64 = 0.0; condition::Bool = false
        cmax::String = ""; do_moist::Bool = false

        # Loop from bottom upwards (over cell-edges)
        for i in range(start=atmos.nlev_l-1, step=-1, stop=2)

            # Profile lapse rate: d(ln T)/d(ln P)
            grad_pr = ( log(atmos.tmp[i-1]/atmos.tmp[i]) )/( log(atmos.p[i-1]/atmos.p[i]) )

            # Optionally skip low pressures 
            if atmos.pl[i] <= pmin
                continue
            end

            m1 = atmos.layer_mass[i-1]
            m2 = atmos.layer_mass[i]
            mt = m1+m2

            grav = (atmos.layer_grav[i] * m2 + atmos.layer_grav[i-1] * m1)/mt
            mu   = (atmos.layer_mmw[i]  * m2 + atmos.layer_mmw[i-1]  * m1)/mt
            c_p  = (atmos.layer_cp[i]   * m2 + atmos.layer_cp[i-1]   * m1)/mt
            tmp = (atmos.tmp[i] * m2 + atmos.tmp[i-1] * m1)/mt
            

            # Dry convection
            grad_ad = (phys.R_gas / mu) / c_p
            condition = (grad_pr > grad_ad)

            # Define moist elsewhere eventually, link to condensates
            # To do:
            #  - Find regions where we are condensing (pass from the condensation scheme)
            #  - Find the condensable with the largest value of (L/RT)^2*x 
            #  - Calculate the adiabatic lapse rate
            #  - Calculate stabilisation w.r.t. the adiabat using criterion
            #  - Set the vapour contents to the new saturated value
            #  - Adjust the dry component to compensate
            
            
            if atmos.condense_any && do_moist
                # Check which condensable species has the largest (L/RT)^2*x

                cmax = ""
                for c in atmos.condensates
                    if atmos.gas_yield[c][i] > 0.0
                        cmax = c
                        break
                    end 
                end

                if !isempty(cmax)
                    mmw_v = atmos.gas_dat[cmax].mmw
                    xv = (atmos.gas_vmr[cmax][i]*m2 + atmos.gas_vmr[cmax][i-1]*m1)/mt
                    beta = phys.get_Lv(atmos.gas_dat[cmax], tmp)*mmw_v/tmp/phys.R_gas

                    grad_ad = (1 - xv + beta*xv) /
                                 (c_p*mu/phys.R_gas * (1-xv) + beta^2 * xv)

                    # Critical value (in vmr, NOT mmr form, only when mmw_v>mmw_d)
                    if mmw_v > (mu - xv*mmw_v)/(1-xv)
                        inhib = phys.R_gas * tmp/L * (1-xv)/(mmw_v/mu - 1)
                        condition = (grad_pr > grad_ad) && xv < inhib
                    else
                        condition = (grad_pr > grad_ad)
                    end
                end 

            end

            # Check instability
            if condition

                rho = (atmos.layer_density[i] * m2 + atmos.layer_density[i-1] * m1)/mt

                atmos.mask_c[i] = true
                
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

                # Characteristic velocity (from Brunt-Vasalla frequency of parcel)
                w = l * sqrt(grav/H * (grad_pr-grad_ad))

                # Dry convective flux
                atmos.flux_cdry[i] = 0.5*rho*c_p*w * atmos.tmpl[i] * (l/H)*(grad_pr-grad_ad)

                # Mixing eddy diffusion coefficient
                atmos.Kzz[i] = w * l
            
            end
        end

        # Check for spurious shallow convection occuring ABOVE condensing regions
        #    If found, reset convective flux to zero AT THIS LAYER ONLY.
        #    This is okay because this shouldn't physically happen, and will only occur 
        #    because of weird numerical issues which only act to make solving difficult.
        for i in 1:atmos.nlev_l-1
            if (!atmos.mask_l[i] && any(atmos.mask_l[i+1:end])) #|| (atmos.mask_l[i] && !atmos.mask_c[i-1] && !atmos.mask_c[i+1])
                atmos.mask_c[i] = false 
                atmos.flux_cdry[i] = 0.0
            end 
        end 

        return nothing
    end # end of mlt

    """
    **Analytical diffusion scheme for condensation and evaporation energy.**

    Integrates from bottom of model upwards. Based on the amount of 
    phase change at each level, a phase change flux is calculated by assuming 
    a fixed condensation timescale. 
    
    Updates fluxes. Requires `atmosphere.handle_saturation` to be called first in the 
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
        E_accum::Float64 = 0.0
        i_dry_top::Int = 1

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

            # Single gas 'evaporation' flux 
            if single 
                # find top of dry region
                for i in 1:atmos.nlev_c
                    if abs(atmos.phs_wrk_df[i]) > 0
                        i_dry_top=i+1
                    end 
                end 

                E_accum = sum(atmos.phs_wrk_df[1:i_dry_top])

                # evaporative flux in dry region
                for i in i_dry_top+1:atmos.nlev_c

                    if E_accum < 1.0e-4
                        # dissipate all of the flux
                        atmos.phs_wrk_df[i] = -E_accum
                        E_accum = 0.0
                        break
                    else
                        # dissipate some fraction of the accumuated flux 
                        atmos.phs_wrk_df[i] = -0.6*E_accum
                        E_accum += atmos.phs_wrk_df[i]
                    end 
                    
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
                if maximum(abs.(atmos.phs_wrk_df[i:end])) < 1.0e-6 
                    # if so, set all phase change fluxes to zero in that region
                    atmos.phs_wrk_fl[i+1:end] .= 0.0
                    break
                end  
            end 
            atmos.phs_wrk_fl[end] = 0.0

            # add energy from this gas to total 
            atmos.flux_l[:] .+= atmos.phs_wrk_fl[:]

            # calculate mask 
            atmos.mask_l[:] .= (abs.(atmos.flux_l[:]) .> 1.0e-10)

        end # go to next condensable 

        return nothing 
    end 

    """
    **Reset energy fluxes to zero.** 
    """ 
    function reset_fluxes!(atmos::atmosphere.Atmos_t)

        # masks 
        fill!(atmos.mask_l, false)
        fill!(atmos.mask_c, false)

        # scalar fluxes 
        atmos.flux_sens = 0.0

        # component fluxes 
        fill!(atmos.flux_cdct, 0.0)
        fill!(atmos.flux_cdry, 0.0)
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

        # total fluxes
        fill!(atmos.flux_dif, 0.0)
        fill!(atmos.flux_tot, 0.0)
    end 

    """
    **Calculate total flux at each level.**

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `latent::Bool`                    include condensation flux
    - `convect::Bool`                   include MLT convection flux
    - `sens_heat::Bool`                 include TKE sensible heat transport 
    - `conduct::Bool`                   include conductive heat transport
    - `convect_sf::Float64`             scale factor applied to convection fluxes
    - `latent_sf::Float64`              scale factor applied to phase change fluxes
    - `calc_cf::Bool=false`             calculate LW contribution function?
    """
    function calc_fluxes!(atmos::atmosphere.Atmos_t, 
                          latent::Bool, convect::Bool, sens_heat::Bool, conduct::Bool;
                          convect_sf::Float64=1.0, latent_sf::Float64=1.0,
                          calc_cf::Bool=false)

        # Reset fluxes
        reset_fluxes!(atmos)

        atmosphere.calc_layer_props!(atmos)

        # +Condensation and evaporation
        if atmos.condense_any && latent 
            # Handle rainout 
            if atmos.gas_num > 1
                atmosphere.handle_saturation!(atmos)
            end

            # Calc flux
            energy.condense_diffuse!(atmos)

            # Modulate?
            atmos.flux_l *= latent_sf

            # Add to total flux 
            atmos.flux_tot += atmos.flux_l

            # Restore mixing ratios
            for g in atmos.gas_names
                atmos.gas_vmr[g][:] .=  atmos.gas_ovmr[g][:]
            end 
        end 

        atmosphere.calc_layer_props!(atmos)

        # +Radiation
        energy.radtrans!(atmos, true, calc_cf=calc_cf)
        energy.radtrans!(atmos, false)
        atmos.flux_tot += atmos.flux_n

        # +Dry convection
        if convect
            # Calc flux
            energy.mlt!(atmos)

            # Modulate?
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

    # Dry convective adjustment (returning the temperature tendency)
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
                
                cp = 0.5 * ( atmos.layer_cp[i-1] * atmos.layer_mmw[i-1] + 
                             atmos.layer_cp[i] * atmos.layer_mmw[i])
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
                
                cp = 0.5 * ( atmos.layer_cp[i-1] * atmos.layer_mmw[i-1] + 
                             atmos.layer_cp[i] * atmos.layer_mmw[i])
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

