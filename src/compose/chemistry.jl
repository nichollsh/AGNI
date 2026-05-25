# This file is part of AGNI. License is GPL-3.0: https://www.gnu.org/licenses

"""
This module handles chemistry, condensation, and evaporation.

Note the important distinctions between the variables which store atmospheric composition.
 * `gas_ovmr` stores VMRs inputted by the user, which are usually constant in height.
 * `gas_cvmr` stores the VMRs calculated by fastchem, which are used as the starting point for condensation calculations.
 * `gas_vmr` stores the runtime gas volume mixing ratios, after all calculations are performed.
"""
module chemistry

    # System libraries
    using Printf
    using Logging

    # Local files
    import ..consts: SMALLFLOAT
    import ..phys
    import ..atmosphere
    import ..species
    import ..ocean
    import ..fastchem


    """
    **Normalise gas VMRs, keeping condensates unchanged**

    Only acts on a single model level.

    Parameters:
    - `atmos::atmosphere.Atmos_t`       atmosphere structure
    - `i::Int64`                        level index to normalise
    """
    function normalise_vmrs!(atmos::atmosphere.Atmos_t, i::Int64)
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
            @. atmos.gas_vmr[g]  = atmos.gas_ovmr[g]
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

            dp = 0.0

            # If de-mixed
            if atmos.demixing && (atmos.tmp_surf < species.get_Tdemix(atmos.gas_dat[c], atmos.p_boa, atmos.gas_vmr[c][end]))
                # rainout completely, and skip condensation
                dp = -p_gas[c]
                atmos.p_boa += dp
                @debug @sprintf("            %s de-mixed, partial pressure += %+.3e bar", c, dp/1e5)
                any_changed = true

            # If supercritical, do nothing
            elseif atmos.tmp_surf > atmos.gas_dat[c].T_crit
                continue

            # Else, handle saturation normally
            else

                # Calculate partial pressure and saturation pressure for this condensable
                #     Work out amount of sub(+) or super(-) saturation
                dp = species.get_Psat(atmos.gas_dat[c], atmos.tmp_surf) - p_gas[c]

                # Negligible
                if abs(dp) < SMALLFLOAT*10.0
                    continue

                # Super-saturated at the surface...
                elseif dp < 0
                    @debug @sprintf("            %s super-saturated, partial pressure += %+.3e bar", c, dp/1e5)

                # Sub-saturated at the surface...
                #     work out change in partial pressure based on initial reservoir amount
                else
                    dp = min(dp, atmos.ocean_ini[c] * atmos.grav_surf / (atmos.gas_dat[c].mmw/atmos.layer_μ[end]))
                    @debug @sprintf("            %s sub-saturated, partial pressure += %+.3e bar", c, dp/1e5)
                end

                # Record that at least one component has as changed
                any_changed = true
            end

            # Reduce or increase total pressure and partial pressure
            atmos.p_boa += dp
            p_gas[c] += dp

            # Change in surface reservoir, with opposite sign to change in pressure
            atmos.ocean_tot[c] -= (dp / atmos.grav_surf) * (atmos.gas_dat[c].mmw/atmos.layer_μ[end])

        end # end condensable

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

    Does not update aerosols or clouds, which should be handled separately.

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

        # Handle condensation (loop from bottom up)
        for i in range(start=atmos.nlev_c, stop=1, step=-1)

            # For each condensate
            for c in atmos.condensates

                # Handle de-mixing
                if atmos.demixing && (atmos.tmp[i] < species.get_Tdemix(atmos.gas_dat[c], atmos.p[i], atmos.gas_vmr[c][i]))
                    # rainout completely at this layer, and skip condensation
                    atmos.gas_vmr[c][i] = 0.0
                    atmos.gas_sat[c][i] = true
                end

                # check criticality
                supcrit = atmos.tmp[i] > atmos.gas_dat[c].T_crit+1.0e-5

                # saturation mixing ratio
                x_sat = species.get_Psat(atmos.gas_dat[c], atmos.tmp[i]) / atmos.p[i]

                # Apply cold trap, implicitly accounting for condensable not making it
                #   this high up by vertical dynamical-transport processes.
                if (atmos.gas_vmr[c][i] > maxvmr[c]) && atmos.coldtrap
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
                                            (atmos.g[i] * atmos.layer_μ[i])

                    # set new vmr to saturated value
                    #   this will always be <= to the current value
                    atmos.gas_vmr[c][i] = x_sat

                    # store vmr for cold trapping at levels above this one
                    maxvmr[c] = x_sat

                    # flag condensate as actively condensing at this level
                    atmos.gas_sat[c][i] = true

                end # end saturation check

                # recalculate layer properties after raining-out this species
                atmosphere.calc_layer_props!(atmos)

            end # end condensate

            # Re-normalise VMRs, keeping condensate VMRs fixed
            normalise_vmrs!(atmos, i)

        end # end i levels

        # Ensure that all yields are positive, at this point in the code
        for c in atmos.condensates
            clamp!(atmos.cond_yield[c], 0.0, Inf)
        end

        # Work out total condensate yield  and do evaporation in lower layers
        for c in atmos.condensates

            # set to zero at TOA
            atmos.cond_accum[c] = 0.0

            # skip here if no evaporation
            if atmos.evap_efficiency < SMALLFLOAT*10
                continue
            end

            # no rain? go to next condensable
            if sum(atmos.cond_yield[c]) < SMALLFLOAT*10
                continue
            end

            # loop from top down (rain always goes downwards)
            for j in 1:atmos.nlev_c

                # raining in this layer...
                #     don't evaporate
                #     add condensate to total budget
                if atmos.cond_yield[c][j] > 0.0
                    atmos.cond_accum[c] += atmos.cond_yield[c][j]
                    continue
                end

                # in a dry layer...

                # skip if no rain entering from above
                if atmos.cond_accum[c] < SMALLFLOAT*10
                    continue
                end

                # exit loop if supercritical, because then condensate mixes miscibly
                if atmos.tmp[j] >= atmos.gas_dat[c].T_crit
                    atmos.cond_accum[c] = 0.0
                    break
                end

                # change in partial pressure that would saturate
                dp_sat = species.get_Psat(atmos.gas_dat[c], atmos.tmp[j]) -
                                                atmos.gas_vmr[c][j]*atmos.p[j]

                # production of gas mass (kg/m2) that would saturate
                dm_sat = atmos.gas_dat[c].mmw * dp_sat/
                                        (atmos.g[j] * atmos.layer_μ[j])

                # Evaporation efficiency factor
                #   This is how close the layer can be brought to saturation by evap.
                #   In reality, this would depend on the microphysical processes.
                dm_sat *= atmos.evap_efficiency

                # don't evaporate more rain than the total available
                dm_sat = min(dm_sat, atmos.cond_accum[c])

                # offset condensate yield at this level by the evaporation
                atmos.cond_yield[c][j] -= dm_sat

                # convert evaporated mass back to partial pressure
                dp_sat = dm_sat * atmos.g[j] * atmos.layer_μ[j] / atmos.gas_dat[c].mmw

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

            # Layer properties
            atmosphere.calc_layer_props!(atmos)

        end # end loop over condensates

        # Layer properties
        for i in 1:atmos.nlev_c
            normalise_vmrs!(atmos, i)
        end
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
    - `state::Int64`                    fastchem state (0: success, 1: critical_fail, 2: elem_fail, 3: conv_fail, 4: both_fail)
    """
    function _chem_gas!(atmos::atmosphere.Atmos_t, write_cfg::Bool)::Int

        @debug "Running equilibrium chemistry"

        # Return code
        state::Int64 = 0

        # Check fastchem enabled
        if !atmos.flag_fastchem
            @warn "Fastchem is not enabled but chemistry was called. Have you set FC_DIR?"
            return 1
        end

        # Check minimum temperature
        if maximum(atmos.tmpl) < atmos.fastchem_floor
            @debug "Whole atmosphere too cold for FC (fc_floor=$(atmos.fastchem_floor))"
        end

        # Call fastchem wrapper
        state = fastchem.run_fastchem!(atmos, write_cfg)

        # Also record this result in gas_cvmr dictionary
        for g in atmos.gas_names
            @. atmos.gas_cvmr[g] = atmos.gas_vmr[g]
        end

        # recalculate layer properties
        atmosphere.calc_layer_props!(atmos)

        # Warn on failure
        if state > 0
            @warn "FastChem internal failure; elements may not be conserved (state=$state)"
        end

        # See docstring for return codes
        return state
    end

    """
    **Reset gas VMRs to post-chemistry values, overwriting condensate effect on VMR**

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    """
    function reset_to_chem!(atmos::atmosphere.Atmos_t)
        for g in atmos.gas_names
            @. atmos.gas_vmr[g] = atmos.gas_cvmr[g]
        end
        return nothing
    end

    """
    **Run condensation and chemistry schemes as required.**

    This function is designed as a wrapper for appropriately handling these three
    schemes together in the correct order, so that variables are appropriately updated.
    Steps:
    1. call `_sat_surf!` to handle saturation at the surface, and adjust the pressure grid
    2. call `_chem_gas!` to handle gas-phase chemistry in the column
    3. call `_sat_aloft!` to handle saturation aloft, above the surface.
    4. update aerosols and clouds

    Arguments:
    - `atmos::Atmos_t`       the atmosphere struct instance to be used.
    - `do_surf::Bool`        do saturation cond/evap at surface
    - `do_chem::Bool`        do thermochemistry in the column
    - `do_aloft::Bool`       do saturation cond/evap aloft

    Returns:
    - `state::Int64`         fastchem state (0: success, 1: critical_fail, 2: elem_fail, 3: conv_fail, 4: both_fail)
    """
    function calc_composition!(atmos::atmosphere.Atmos_t,
                                    do_surf::Bool, do_chem::Bool, do_aloft::Bool)::Int

        state::Int64 = 0

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

        # keep cold trapping rainout effect on abundances?
        # if !atmos.coldtrap
        #     reset_to_chem!(atmos)
        # end

        return state
    end


end # end module
