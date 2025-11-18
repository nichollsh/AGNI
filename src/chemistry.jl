# Contains things relating to atmospheric chemistry

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

"""
This module handles chemistry, condensation, and evaporation.

Note the important distinctions between the variables which store atmospheric composition.
 * `gas_ovmr` stores VMRs inputted by the user, which are usually constant in height.
 * `gas_vmr` stores the runtime gas volume mixing ratios, after all calculations are performed.

"""
module chemistry

    # System libraries
    using Printf
    using Logging
    import DelimitedFiles:readdlm

    # Local files
    import ..atmosphere
    import ..phys
    import ..ocean

    # Constants
    const COND_EPS::Float64     = 1e-10     # negligible amount of condesate [kg/m^2]
    const SMOOTH_SCALE::Float64 = 12.0      # smoothing scale for fastchem floor

    """
    **Normalise gas VMRs, keeping condensates unchanged**

    Only acts on a single model level.

    Parameters:
    - `atmos::atmosphere.Atmos_t`       atmosphere structure
    - `i::Int`                          level index to normalise
    """
    function normalise_vmrs!(atmos::atmosphere.Atmos_t, i::Int)
        # Work variables
        x_con::Float64 =    0.0
        x_dry::Float64 =    0.0
        x_dry_old::Float64= 0.0

        # Work out total VMR of all condensing gases
        for c in atmos.condensates
            if atmos.gas_sat[c][i]
                x_con += atmos.gas_vmr[c][i]
            end
        end

        # Calculate current and target dry VMRs
        x_dry = 1.0 - x_con
        x_dry_old = 0.0
        for g in atmos.gas_names
            # skip condensing gases, since their VMR is set by saturation
            if !atmos.gas_sat[g][i]
                x_dry_old += atmos.gas_vmr[g][i]
            end
        end

        # Renormalise VMR to unity, scaling DRY COMPONENTS ONLY
        for g in atmos.gas_names
            if !atmos.gas_sat[g][i]
                atmos.gas_vmr[g][i] *= x_dry / x_dry_old
            end
        end

        # Check total VMR at this level
        x_tot = 0.0
        for g in atmos.gas_names
            x_tot += atmos.gas_vmr[g][i]
        end
        if abs(x_tot - 1.0) > 1.0e-5
            @warn @sprintf("Mixing ratios sum to %.6e (level %d)",x_tot,i)
        end

        return nothing
    end


    """
    **Reset mixing ratios, pressures, and oceans to their original values**
    """
    function restore_composition!(atmos::atmosphere.Atmos_t)

        # Mixing ratios to original values
        # Oceans to original values
        for g in atmos.gas_names
            @. atmos.gas_vmr[g] = atmos.gas_ovmr[g]
            @. atmos.gas_cvmr[g] = atmos.gas_ovmr[g]

            atmos.ocean_tot[g] = atmos.ocean_ini[g]
        end

        # Pressure grid
        atmos.p_boa = atmos.p_oboa
        atmosphere.generate_pgrid!(atmos)

        # Layer properties
        atmosphere.calc_layer_props!(atmos)
    end


    """
    **Handle sub/super-saturation at the surface, forming oceans.**

    NOTE: this function does not reset to `atmos` to its original state before operating.

    If a condensable is supersaturated at the surface, its partial surface pressure is
    reduced to exactly saturation. The condensed mass is added to the ocean reservoir.
    If a condensable is subsaturated at the surface, its partial pressure is increased
    to exactly saturation, subject to the amount of ocean available. Reduces or increases
    the total surface pressure accordingly.

    Accounts for mmw ratio when converting from partial pressures to masses:
        m_i = (p_i/g) * (mu_i / mu_tot)
    where m_i has units [kg/m^2] and mu_tot is the total mixture MMW.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.

    Returns:
    - `any_changed::Bool`       did any component condense or evaporate?
    """
    function _sat_surf!(atmos::atmosphere.Atmos_t)::Bool

        # Define work variables
        any_changed::Bool = false
        dp::Float64 = 0.0  # change in surface partial pressure [Pa]

        # Populate partial pressure dictionary for ALL gases
        p_gas::Dict{String, Float64} = Dict()
        for gas in atmos.gas_names
            p_gas[gas] =  atmos.gas_vmr[gas][end] * atmos.p_boa
        end

        # For each condensable volatile...
        for c in atmos.condensates

            # If supercritical, skip
            if atmos.tmp_surf > atmos.gas_dat[c].T_crit
                continue
            end

            # Calculate partial pressure and saturation pressure for this condensable
            #     Work out amount of sub(+) or super(-) saturation
            dp = phys.get_Psat(atmos.gas_dat[c], atmos.tmp_surf) - p_gas[c]

            # Negligible
            if abs(dp) < 1e-10
                continue

            # Super-saturated at the surface...
            elseif dp < 0
                @debug @sprintf("            %s super-saturated, partial pressure += %+.3f bar", c, dp/1e5)

            # Sub-saturated at the surface...
            #     work out change in partial pressure based on initial reservoir amount
            else
                dp = min(dp, atmos.ocean_ini[c] * atmos.grav_surf / (atmos.gas_dat[c].mmw/atmos.layer_μ[end]))
                @debug @sprintf("            %s sub-saturated, partial pressure += %+.3f bar", c, dp/1e5)
            end

            # Record that at least one component has as changed
            any_changed = true

            # Reduce or increase total pressure and partial pressure
            atmos.p_boa += dp
            p_gas[c] += dp

            # Change in surface reservoir, with opposite sign to change in pressure
            atmos.ocean_tot[c] -= (dp / atmos.grav_surf) * (atmos.gas_dat[c].mmw/atmos.layer_μ[end])
        end

        # Exit now if nothing needs to be done
        if !any_changed
            return any_changed
        end

        # Do now allow p_boa to become smaller than p_toa
        if atmos.p_boa < atmos.p_toa * atmosphere.PRESSURE_RATIO_MIN
            @warn @sprintf("Surface pressure too low! Calculated %.3f bar",atmos.p_boa/1e5)
            atmos.p_boa = atmos.p_toa * atmosphere.PRESSURE_RATIO_MIN
        end

        # Recalculate VMRs from partial pressures (dalton's law)
        #     Sets to well-mixed composition
        for gas in atmos.gas_names
            fill!(atmos.gas_vmr[gas], p_gas[gas]/atmos.p_boa)
        end

        # Generate new pressure grid with updated p_boa
        @debug @sprintf("            new p_boa = %.3f bar",atmos.p_boa/1e5)
        atmosphere.generate_pgrid!(atmos)

        # Calculate new values for layer properties
        atmosphere.calc_layer_props!(atmos)

        return any_changed
    end


    """
    **Handle sub/super-saturation aloft; required for latent heat flux calcaulation**

    Adjust gas VMRs according to saturation and cold-trap requirements.

    Volatiles which are allowed to condense are rained-out at condensing levels
    until the gas is exactly saturated, not supersaturated. If evaporation is enabled here,
    it will lead to enhanced mixing ratios at deeper levels as rain is converted back
    into gas from a liquid state.

    Any rain that reaches the surface without being evaporated leads to ocean formation.
    Disabling evaporation will lead to unclosed energy budget when calculation of latent
    heat fluxes is performed.

    This function can be called *after* the fastchem calculation.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    """
    function _sat_aloft!(atmos::atmosphere.Atmos_t)

        # Work arrays
        maxvmr::Dict{String, Float64} = Dict{String, Float64}() # max running VMR for each condensable
        x_sat::Float64 = 0.0
        supcrit::Bool =  false

        # Reset mixing ratios to post-chemistry values, since we will write to gas_vmr
        # Reset phase change flags
        # Reset condensation yield values
        for g in atmos.gas_names
           fill!(atmos.gas_sat[g],    false)
           fill!(atmos.cond_yield[g], 0.0)  # kg/m2 of condensate produced at levels
        end

        # Set initial maximum value as surface composition (for cold trapping)
        for c in atmos.condensates
            maxvmr[c] = atmos.gas_vmr[c][end]
        end

        # Reset water cloud
        fill!(atmos.cloud_arr_r, 0.0)
        fill!(atmos.cloud_arr_l, 0.0)
        fill!(atmos.cloud_arr_f, 0.0)

        # Handle condensation (loop from bottom up)
        for i in range(start=atmos.nlev_c, stop=1, step=-1)

            # For each condensate
            for c in atmos.condensates

                # check criticality
                supcrit = atmos.tmp[i] > atmos.gas_dat[c].T_crit+1.0e-5

                # saturation mixing ratio
                x_sat = phys.get_Psat(atmos.gas_dat[c], atmos.tmp[i]) / atmos.p[i]

                # Apply cold trap, implicitly accounting for condensable not making it
                #   this high up by vertical dynamical-transport processes.
                if atmos.gas_vmr[c][i] > maxvmr[c]
                    atmos.gas_vmr[c][i] = maxvmr[c]
                    atmos.gas_sat[c][i] = true
                end

                # condense here if supersaturated
                if (atmos.gas_vmr[c][i] > x_sat) && !supcrit

                    # set rainout [kg/m2]
                    #    based on difference in mole fraction between
                    #    current value and saturation
                    atmos.cond_yield[c][i] = atmos.gas_dat[c].mmw*atmos.p[i]*
                                            (atmos.gas_vmr[c][i] - x_sat)/
                                            (atmos.layer_grav[i] * atmos.layer_μ[i])

                    # set new vmr to saturated value
                    #   this will always be <= to the current value
                    atmos.gas_vmr[c][i] = x_sat

                    # store vmr for cold trapping at levels above this one
                    maxvmr[c] = x_sat

                    # flag condensate as actively condensing at this level
                    atmos.gas_sat[c][i] = true

                end # end saturation check

                # recalculate layer properties after raining-out this species
                # atmosphere.calc_layer_props!(atmos)

            end # end condensate

            normalise_vmrs!(atmos, i)
        end # end i levels

        # Ensure that all yields are positive at this point
        for c in atmos.condensates
            clamp!(atmos.cond_yield[c], 0.0, Inf)
        end

        # Work out total condensate yield  and do evaporation in lower layers
        for c in atmos.condensates

            # set to zero at TOA
            atmos.cond_accum[c] = 0.0

            # no rain? go to next condensable
            if sum(atmos.cond_yield[c]) < COND_EPS
                continue
            end

            # loop from top down (rain always goes downwards)
            for j in 1:atmos.nlev_c

                # set negligible rain to zero
                if atmos.cond_yield[c][j] < COND_EPS
                    atmos.cond_yield[c][j] = 0.0
                end

                # raining in this layer...
                #     don't evaporate
                #     add condensate to total budget
                if atmos.cond_yield[c][j] > 0.0
                    atmos.cond_accum[c] += atmos.cond_yield[c][j]
                    continue
                end

                # in a dry layer...

                # skip if no rain entering from above
                if atmos.cond_accum[c] < eps(1.0)
                    continue
                end

                # exit loop if supercritical, because then condensate mixes miscibly
                if atmos.tmp[j] >= atmos.gas_dat[c].T_crit
                    atmos.cond_accum[c] = 0.0
                    break
                end

                # change in partial pressure that would saturate
                dp_sat = phys.get_Psat(atmos.gas_dat[c], atmos.tmp[j]) -
                                                atmos.gas_vmr[c][j]*atmos.p[j]

                # production of gas mass (kg/m2) that would saturate
                dm_sat = atmos.gas_dat[c].mmw * dp_sat/
                                        (atmos.layer_grav[j] * atmos.layer_μ[j])

                # Evaporation efficiency factor
                #   This is how close the layer can be brought to saturation by evap.
                #   In reality, this would depend on the microphysical processes.
                dm_sat *= atmos.evap_efficiency

                # don't evaporate more rain than the total available
                dm_sat = min(dm_sat, atmos.cond_accum[c])

                # offset condensate yield at this level by the evaporation
                atmos.cond_yield[c][j] -= dm_sat

                # convert evaporated mass back to partial pressure
                dp_sat = dm_sat * atmos.layer_grav[j] *
                                            atmos.layer_μ[j] / atmos.gas_dat[c].mmw

                # convert change in partial pressure to change in vmr
                atmos.gas_vmr[c][j] += dp_sat / atmos.p[j]

                # flag as 'saturated' - somewhat a misnomer
                atmos.gas_sat[c][j] = true

                # Recalculate layer mmw
                atmos.layer_μ[j] = 0.0
                for g in atmos.gas_names
                    atmos.layer_μ[j] += atmos.gas_vmr[g][j] * atmos.gas_dat[g].mmw
                end

                # recalculate total rain correspondingly
                atmos.cond_accum[c] -= dm_sat

            end # go to next j level (below)

        end # end loop over condensates

        # Set water clouds at levels where condensation occurs
        if "H2O" in atmos.condensates
            for i in 1:atmos.nlev_c-1
                if atmos.cond_yield["H2O"][i] > 0.0
                    # liquid water content (take ratio of mass surface densities [kg/m^2])
                    atmos.cloud_arr_l[i] = (atmos.cond_yield["H2O"][i]*atmos.cloud_alpha) /
                                                atmos.layer_mass[i]

                    # droplet radius and area fraction (fixed values)
                    atmos.cloud_arr_r[i] = atmos.cloud_val_r
                    atmos.cloud_arr_f[i] = atmos.cloud_val_f
                end
            end
        end

        # Layer properties
        atmosphere.calc_layer_props!(atmos)

        return nothing
    end


    """
    **Calculate gas-phase composition at dry thermochemical equilibrium.**

    Uses FastChem to calculate the gas composition at each level of the atmosphere.
    Volatiles are converted to bulk elemental abundances, which are then provided to
    FastChem alongside the temperature/pressure profile. FastChem is currently called as an
    executable, which is not optimal.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `write_cfg::Bool`                 write config and elements

    Returns:
    - `state::Int`                      fastchem state (0: success, 1: critical_fail, 2: elem_fail, 3: conv_fail, 4: both_fail)
    """
    function _chem_gas!(atmos::atmosphere.Atmos_t, write_cfg::Bool)::Int

        @debug "Running equilibrium chemistry"

        # Return code
        state::Int = 0

        # Check fastchem enabled
        if !atmos.flag_fastchem
            @warn "Fastchem is not enabled but chemistry was called. Have you set FC_DIR?"
            return 1
        end

        # Check minimum temperature
        if maximum(atmos.tmpl) < atmos.fastchem_floor
            @warn "The entire temperature profile is too cold for FastChem"
        end

        count_elem_nonzero::Int = 0

        # Write config (fastchem is quite particular about the format)
        if write_cfg || !isfile(atmos.fastchem_conf)
            open(atmos.fastchem_conf,"w") do f
                write(f,"#Atmospheric profile input file \n")
                write(f,atmos.fastchem_prof*" \n\n")

                write(f,"#Chemistry calculation type (gas phase only = g, equilibrium condensation = ce, rainout condensation = cr) \n")
                write(f,"g \n\n")

                write(f,"#Chemistry output file \n")
                write(f,atmos.fastchem_chem*" "*atmos.fastchem_cond*" \n\n")

                write(f,"#Monitor output file \n")
                write(f,atmos.fastchem_moni*" \n\n")

                write(f,"#FastChem console verbose level (1 - 4); 1 = almost silent, 4 = detailed console output \n")
                write(f,"1 \n\n")

                write(f,"#Output mixing ratios (MR) or particle number densities (ND, default) \n")
                write(f,"ND \n\n")

                write(f,"#Element abundance file  \n")
                write(f,atmos.fastchem_elem*" \n\n")

                write(f,"#Species data files    \n")
                logK = joinpath(atmos.FC_DIR, "input/","logK/")
                write(f,joinpath(logK,"logK.dat")*" "*joinpath(logK,"logK_condensates.dat")*" \n\n")

                write(f,"#Accuracy of chemistry iteration \n")
                write(f,@sprintf("%.3e \n\n", atmos.fastchem_xtol_chem))

                write(f,"#Accuracy of element conservation \n")
                write(f,@sprintf("%.3e \n\n", atmos.fastchem_xtol_elem))

                write(f,"#Max number of chemistry iterations  \n")
                write(f,@sprintf("%d \n\n", atmos.fastchem_maxiter_chem))

                write(f,"#Max number internal solver iterations  \n")
                write(f,@sprintf("%d \n\n", atmos.fastchem_maxiter_solv))
            end
        end # end write config

        # Write metallicites
        if write_cfg || !isfile(atmos.fastchem_elem)

            # Reset metallicities
            atmos.metal_calc = Dict{String,Float64}()

            # Metallicities provided by user
            if !isempty(atmos.metal_orig)
                @debug "Elements set by user-provided metallicities"

                # copy original to calculated; set elem to zero if it was not provided
                for e in phys.elems_standard
                    atmos.metal_calc[e] = get(atmos.metal_orig, e, 0.0)
                end
                atmos.metal_calc["H"] = 1.0

            # Not provided -- calculate from composition at surface
            else
                @debug "Elements set by surface gas composition"

                # Calculate elemental abundances from surface mixing ratios [molecules/m^3]
                #   assuming ideal gas: N/V = P*x/(Kb*T) , where x is the VMR
                N_t = zeros(Float64, length(phys.elems_standard)) # total atoms in all gases
                N_g = zeros(Float64, length(phys.elems_standard)) # atoms in current gas
                #    loop over gases
                for gas in atmos.gas_names
                    fill!(N_g, 0.0)

                    # count atoms in this gas
                    d = phys.count_atoms(gas)
                    for (i,e) in enumerate(phys.elems_standard)
                        if haskey(d, e)
                            N_g[i] += d[e] # N_g stores num of atoms in this gas
                        end
                    end

                    # Get gas abundance from original VMR value
                    #    scale number of atoms by the abundance of the gas
                    N_g *= atmos.gas_vmr[gas][end] * atmos.p[end] / (phys.k_B * atmos.tmp[end])

                    # Add atoms from this gas to total atoms in the mixture
                    N_t += N_g
                end

                # Convert elemental abundances to metallicity number ratios, rel to hydrogen
                for (i,e) in enumerate(phys.elems_standard)
                    atmos.metal_calc[e] = N_t[i]/N_t[1]
                end
            end

            # Write metallicities to FC input file in the required format
            #     number densities normalised relative to hydrogen
            #     for each element `e`, value = log10(N_e/N_H) + 12
            open(atmos.fastchem_elem,"w") do f
                write(f,"# Elemental abundances file written by AGNI \n")
                for e in phys.elems_standard

                    # skip this element if its abundance is too small
                    if atmos.metal_calc[e] < 1.0e-30
                        continue
                    end

                    # normalise abundance relative to hydrogen
                    write(f, @sprintf("%s    %.3f \n",e,log10(atmos.metal_calc[e]) + 12.0))
                    count_elem_nonzero += 1
                end
            end

        end # end write metallicities

        """
        **Smoothly transform temperature around floor**

        Approximately returns _x for _x > atmos.fastchem_floor.
        Approximately returns atmos.fastchem_floor for _x < atmos.fastchem_floor.
        Transitions discontinuously for a->Inf
        """
        function _transform_floor(_x::Float64)::Float64
            d = 1.0 / ( 1.0 + exp(-SMOOTH_SCALE*(_x-atmos.fastchem_floor)/atmos.fastchem_floor) )
            return _x*d + atmos.fastchem_floor*(1-d)
        end

        # Write PT profile every time
        open(atmos.fastchem_prof,"w") do f
            write(f,"# AGNI temperature structure \n")
            write(f,"# bar, kelvin \n")
            for i in 1:atmos.nlev_c
                write(  f,
                        @sprintf("%.6e    %.6e \n",
                            atmos.p[i]*1e-5,
                            _transform_floor(atmos.tmp[i])
                            )
                     )
            end
        end

        # Run fastchem
        run(pipeline(`$(atmos.fastchem_exec) $(atmos.fastchem_conf)`, stdout=devnull))

        # Check monitor output
        data = readdlm(atmos.fastchem_moni, '\t', String)
        fail_elem::String = ""
        fail_conv::String = ""
        for i in 1:atmos.nlev_c
            if data[i+1,6][1] == 'f'
                fail_elem *= @sprintf("%d ",i)
            end
            if data[i+1,5][1] == 'f'
                fail_conv *= @sprintf("%d ",i)
            end
        end
        if !isempty(fail_elem)
            @debug "Element conservation failed at levels  "*fail_elem
            state = 2
        end
        if !isempty(fail_conv)
            @debug "FastChem solver failed at levels  "*fail_conv
            if state == 2
                state = 4
            else
                state = 3
            end
        end

        # Get gas chemistry output
        if !isfile(atmos.fastchem_chem)
            @error "Could not find fastchem output"
            return 1
        end
        (data,head) = readdlm(atmos.fastchem_chem, '\t', Float64, header=true)
        data = transpose(data)  # convert to: gas, level

        # Clear VMRs
        for g in atmos.gas_names
            fill!(atmos.gas_vmr[g],  0.0)
        end

        # Parse gas chemistry
        g_fc::String = atmosphere.UNSET_STR
        d_fc::Dict = Dict{String, Int}()
        g_in::String = atmosphere.UNSET_STR
        match::Bool = false
        N_t = data[4,:] # at each level: sum of gas number densities

        for (i,h) in enumerate(head)  # for each column (gas)

            # skip columns (p, T, ntot, ngas, mu, and elemental abundances)
            if i <= 5+count_elem_nonzero
                continue
            end

            # parse name
            g_fc = rstrip(lstrip(h))
            if occursin("_", g_fc)
                g_fc = split(g_fc, "_")[1]
            end
            g_fc = replace(g_fc, "cis"=>"", "trans"=>"")

            match = false

            # firstly, check if we have the FC name already stored
            for g in atmos.gas_names
                if atmos.gas_dat[g].fastchem_name == g_fc
                    match = true
                    g_in = g
                    break
                end
            end

            # not stored => search based on atoms
            if !match
                d_fc = phys.count_atoms(g_fc)  # get atoms dict from FC name

                for g in atmos.gas_names
                    if phys.same_atoms(d_fc, atmos.gas_dat[g].atoms)
                        match = true
                        g_in = g
                        break
                    end
                end
            end

            # matched?
            if match
                N_g = data[i,:]  # number densities for this gas
                @. atmos.gas_vmr[g_in] += N_g / N_t    # VMR for this gas
            end
        end

        # Do not renormalise mixing ratios, since this is done by fastchem
        # If we are missing gases then that's okay.

        # Find where T(p) drops below fastchem_floor temperature.
        # Make sure that regions above that use the reasonable VMR values
        i_trunc::Int = 0
        for i in range(start=atmos.nlev_c, stop=1, step=-1)
            if atmos.tmp[i] < atmos.fastchem_floor
               i_trunc = i
               break
            end
        end
        if i_trunc > 0
            # @warn @sprintf("Temperature below FC floor, at p < %.1e Pa", atmos.p[i_trunc])
            i_trunc = min(i_trunc, atmos.nlev_c-1)
            for g in atmos.gas_names
                atmos.gas_vmr[g][1:i_trunc] .= atmos.gas_vmr[g][i_trunc+1]
            end
        end

        # Also record this result in gas_cvmr dictionary
        for g in atmos.gas_names
            @. atmos.gas_cvmr[g] = atmos.gas_vmr[g]
        end

        # recalculate layer properties
        atmosphere.calc_layer_props!(atmos)

        # See docstring for return codes
        return state
    end

    """
    **Run condensation and chemistry schemes as required.**

    This function is designed as a wrapper for appropriately handling these three
    schemes together in the correct order, so that variables are appropriately updated.
    Steps:
    1. call `_sat_surf!` to handle saturation at the surface, and adjust the pressure grid
    2. call `_chem_gas!` to handle gas-phase chemistry in the column
    3. call `_sat_aloft!` to handle saturation aloft, above the surface.

    Arguments:
    - `atmos::Atmos_t`       the atmosphere struct instance to be used.
    - `do_surf::Bool`        do saturation cond/evap at surface
    - `do_chem::Bool`        do thermochemistry in the column
    - `do_aloft::Bool`       do saturation cond/evap aloft

    Returns:
    - `state::Int`           fastchem state (0: success, 1: critical_fail, 2: elem_fail, 3: conv_fail, 4: both_fail)
    """
    function calc_composition!(atmos::atmosphere.Atmos_t,
                                    do_surf::Bool, do_chem::Bool, do_aloft::Bool)::Int

        state::Int = 0

        # reset composition
        restore_composition!(atmos)

        # surface saturation
        if do_surf
            _sat_surf!(atmos)
        end

        # aloft gas-phase chemistry
        if do_chem
            state = _chem_gas!(atmos, true)
        end

        # aloft saturation
        if do_aloft
            _sat_aloft!(atmos)

            # sum total condensate as surface + aloft
            for c in atmos.condensates
                atmos.ocean_tot[c] += atmos.cond_accum[c]
            end
        end

        return state
    end


end # end module
