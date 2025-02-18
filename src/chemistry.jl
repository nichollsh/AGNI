# Contains things relating to atmospheric chemistry

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module chemistry

    # System libraries
    using Printf
    using Logging
    import DelimitedFiles:readdlm

    # Local files
    import ..atmosphere
    import ..phys

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
    **Adjust gas VMRs according to saturation and cold-trap requirements**

    Volatiles which are allowed to condense are rained-out at condensing levels
    until the gas is exactly saturated, not supersaturated. If evaporation is enabled here,
    it will lead to enhanced mixing ratios at deeper levels as rain is converted back
    into gas from a liquid state.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    """
    function handle_saturation!(atmos::atmosphere.Atmos_t)

        # Single gas case does not apply here
        if atmos.gas_num == 1
            return
        end

        # Parameters
        evap_enabled::Bool =        false   # Enable re-vaporation of rain
        evap_efficiency::Float64 =  0.5     # Evaporation efficiency

        # Work arrays
        maxvmr::Dict{String, Float64} = Dict{String, Float64}() # max running VMR for each condensable
        cond_kg::Dict{String,Float64} = Dict{String, Float64}() # condensed kg/m2 for each condensable
        x_sat::Float64 = 0.0
        supcrit::Bool =  false

        # Set maximum value (for cold trapping)
        for c in atmos.condensates
            maxvmr[c] = atmos.gas_vmr[c][end]
        end

        # Reset mixing ratios to surface values
        # Reset phase change flags
        # Reset condensation yield values
        for g in atmos.gas_names
           fill!(atmos.gas_vmr[g][1:end-1], atmos.gas_vmr[g][end])
           fill!(atmos.gas_sat[g],          false)
           fill!(atmos.gas_yield[g],        0.0)
        end

        # Reset water cloud
        fill!(atmos.cloud_arr_r, 0.0)
        fill!(atmos.cloud_arr_l, 0.0)
        fill!(atmos.cloud_arr_f, 0.0)

        # Handle condensation
        for i in range(start=atmos.nlev_c-1, stop=1, step=-1)

            # For each condensate
            for c in atmos.condensates

                # Reset condensation and rain
                cond_kg[c] = 0.0   # kg/m2 of 'c' condensate produced at this level

                # check criticality
                supcrit = atmos.tmp[i] > atmos.gas_dat[c].T_crit+1.0e-5

                # saturation mixing ratio
                x_sat = phys.get_Psat(atmos.gas_dat[c], atmos.tmp[i]) / atmos.p[i]

                # cold trap
                if atmos.gas_vmr[c][i] > maxvmr[c]
                    atmos.gas_vmr[c][i] = maxvmr[c]
                    atmos.gas_sat[c][i] = true
                end

                # condense if supersaturated
                if (atmos.gas_vmr[c][i] > x_sat) && !supcrit

                    # set rainout kg/m2
                    cond_kg[c] = atmos.gas_dat[c].mmw*atmos.p[i]*
                                            (atmos.gas_vmr[c][i] - x_sat)/
                                            (atmos.layer_grav[i] * atmos.layer_μ[i])

                    # condensation yield at this level
                    atmos.gas_yield[c][i] += cond_kg[c]

                    # set new vmr
                    atmos.gas_vmr[c][i] = x_sat

                    # store vmr for cold trapping at levels above this one
                    maxvmr[c] = x_sat

                    # flag condensate as actively condensing at this level
                    atmos.gas_sat[c][i] = true
                    # @printf("%d: %s rain generated %.3e \n", i, c, cond_kg[c])

                    # Set water cloud at this level
                    if c == "H2O"
                        # mass mixing ratio (take ratio of mass surface densities [kg/m^2])
                        atmos.cloud_arr_l[i] = (cond_kg["H2O"] * atmos.cond_alpha) /
                                                    atmos.layer_mass[i]

                        if atmos.cloud_arr_l[i] > 1.0
                            @warn "Water cloud mass mixing ratio is greater than unity (level $i)"
                        end

                        # droplet radius and area fraction (fixed values)
                        atmos.cloud_arr_r[i] = atmos.cloud_val_r
                        atmos.cloud_arr_f[i] = atmos.cloud_val_f
                    end
                end # end saturation check

            end # end condensate

            normalise_vmrs!(atmos, i)
        end # end i levels

        # recalculate layer properties
        atmosphere.calc_layer_props!(atmos)

        # Evaporate rain within unsaturated layers
        if evap_enabled

            total_rain::Float64 = 0.0
            i_top_dry::Int = 1
            i_bot_dry::Int = atmos.nlev_c

            # For each condensable
            for c in atmos.condensates

                # reset dry region
                i_top_dry = 1
                i_bot_dry = atmos.nlev_c

                # accumulate rain
                total_rain = sum(atmos.gas_yield[c])

                # no rain? go to next condensable
                if total_rain < 1.0e-10
                    continue
                end

                # locate top of dry region
                for j in 1:atmos.nlev_c-1
                    if atmos.gas_sat[c][j]
                        # this layer is saturated => dry region must be below it
                        i_top_dry = j+1
                    end
                end

                # locate bottom of dry region
                for j in range(start=atmos.nlev_c, stop=i_top_dry+1, step=-1)
                    if atmos.tmp[j] > atmos.gas_dat[c].T_crit
                        # this layer is supercritical => dry region must be above it
                        i_bot_dry = j-1
                    end
                end
                if i_bot_dry < i_top_dry
                    i_bot_dry = i_top_dry
                end

                # evaporate rain in dry region
                for j in i_top_dry+1:i_bot_dry

                    # evaporate up to saturation
                    dp_sat = phys.get_Psat(atmos.gas_dat[c], atmos.tmp[j]) -
                                                    atmos.gas_vmr[c][j]*atmos.p[j]

                    # Calculate kg/m2 of gas that would saturate layer j
                    dm_sat = atmos.gas_dat[c].mmw * dp_sat/
                                            (atmos.layer_grav[j] * atmos.layer_μ[j])

                    # can we evaporate all rain within this layer?
                    if total_rain < dm_sat
                        # yes, so don't evaporate more rain than the total
                        dm_sat = total_rain
                    end
                    atmos.gas_sat[c][j] = true

                    # Evaporation efficiency factor
                    #   This fraction of the rain that *could* be evaporated
                    #   at this layer *is* converted to vapour in this layer.
                    dm_sat *= evap_efficiency

                    # offset condensate yield at this level by the evaporation
                    atmos.gas_yield[c][j] -= dm_sat

                    # add partial pressure from evaporation to pp at this level
                    dp_sat = dm_sat * atmos.layer_grav[j] *
                                              atmos.layer_μ[j] / atmos.gas_dat[c].mmw

                    # convert extra pp to extra vmr
                    atmos.gas_vmr[c][j] += dp_sat / atmos.p[j]

                    # reduce total rain correspondingly
                    total_rain -= dm_sat

                    # Recalculate layer mmw
                    atmos.layer_μ[j] = 0.0
                    for g in atmos.gas_names
                        atmos.layer_μ[j] += atmos.gas_vmr[g][j] * atmos.gas_dat[g].mmw
                    end

                end # go to next j level (below)

            end # end loop over condensates

        end # end evaporation scheme

        # Recalculate layer properties
        atmosphere.calc_layer_props!(atmos)

        return nothing
    end

    """
    **Calculate composition assuming chemical equilibrium at each level.**

    Uses FastChem to calculate the gas composition at each level of the atmosphere.
    Volatiles are converted to bulk elemental abundances, which are then provided to
    FastChem alongside the temperature/pressure profile. FastChem is currently called as an
    executable, which is not optimal.

    This function DOES NOT automatically recalculate layer properties (e.g. mmw, density).

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `chem_type::Int`                  chemistry type (see wiki)
    - `write_cfg::Bool`                 write config and elements

    Returns:
    - `state::Int`                      fastchem state (0: success, 1: critical_fail, 2: elem_fail, 3: conv_fail, 4: both_fail)
    """
    function fastchem_eqm!(atmos::atmosphere.Atmos_t, chem_type::Int, write_cfg::Bool)::Int

        @debug "Running equilibrium chemistry"

        # Return code
        state::Int = 0

        # Check fastchem enabled
        if !atmos.fastchem_flag
            @warn "Fastchem is not enabled but `fastchem_eqm!` was called. Have you set FC_DIR?"
            return 1
        end

        # Check minimum temperature
        if maximum(atmos.tmpl) < atmos.fastchem_floor
            @warn "Temperature profile is entirely too cold for FastChem. Not doing chemistry."
            return 1
        end

        count_elem_nonzero::Int = 0

        # Paths
        execpath::String = joinpath(atmos.FC_DIR,       "fastchem")             # Executable file
        confpath::String = joinpath(atmos.fastchem_work,"config.input")         # Configuration by AGNI
        elempath::String = joinpath(atmos.fastchem_work,"elements.dat")         # Elements by AGNI
        chempath::String = joinpath(atmos.fastchem_work,"chemistry.dat")        # Chemistry by FastChem

        # Check file exists
        write_cfg = write_cfg || !isfile(confpath) || !isfile(elempath)

        # Write config, elements
        if write_cfg
            # Write config (fastchem is quite particular about the format)
            open(confpath,"w") do f
                write(f,"#Atmospheric profile input file \n")
                write(f,joinpath(atmos.fastchem_work,"pt.dat")*" \n\n")

                type_char = ["g","ce","cr"]
                write(f,"#Chemistry calculation type (gas phase only = g, equilibrium condensation = ce, rainout condensation = cr) \n")
                write(f,"$(type_char[chem_type]) \n\n")

                write(f,"#Chemistry output file \n")
                write(f,joinpath(atmos.fastchem_work,"chemistry.dat")*" "*joinpath(atmos.fastchem_work,"condensates.dat")*" \n\n")

                write(f,"#Monitor output file \n")
                write(f,joinpath(atmos.fastchem_work,"monitor.dat")*" \n\n")

                write(f,"#FastChem console verbose level (1 - 4); 1 = almost silent, 4 = detailed console output \n")
                write(f,"1 \n\n")

                write(f,"#Output mixing ratios (MR) or particle number densities (ND, default) \n")
                write(f,"ND \n\n")

                write(f,"#Element abundance file  \n")
                write(f,joinpath(atmos.fastchem_work,"elements.dat")*" \n\n")

                write(f,"#Species data files    \n")
                logK = joinpath(atmos.FC_DIR, "input/","logK/")
                write(f,joinpath(logK,"logK.dat")*" "*joinpath(logK,"logK_condensates.dat")*" \n\n")

                write(f,"#Accuracy of chemistry iteration \n")
                write(f,@sprintf("%.3e \n\n", atmos.fastchem_xtol))

                write(f,"#Accuracy of element conservation \n")
                write(f,@sprintf("%.3e \n\n", atmos.fastchem_xtol))

                write(f,"#Max number of chemistry iterations  \n")
                write(f,@sprintf("%d \n\n", atmos.fastchem_maxiter*2.5))

                write(f,"#Max number internal solver iterations  \n")
                write(f,@sprintf("%d \n\n", atmos.fastchem_maxiter))
            end

            # Calculate elemental abundances
            # number densities normalised relative to hydrogen
            # for each element X, value = log10(N_X/N_H) + 12
            # N = X(P/(K*T) , where X is the VMR and K is boltz-const
            N_t = zeros(Float64, length(phys.elems_standard))      # total atoms in all gases
            N_g = zeros(Float64, length(phys.elems_standard))      # atoms in current gas
            for gas in atmos.gas_names
                d = phys.count_atoms(gas)
                fill!(N_g, 0.0)
                for (i,e) in enumerate(phys.elems_standard)
                    if e in keys(d)
                        N_g[i] += d[e]
                    end
                end
                # Get gas abundance from original VMR value, since the running
                #    one will be updated using FastChem's output. These will
                #    be normalised later in this function.
                N_g *= atmos.gas_ovmr[gas][atmos.nlev_c] * atmos.p[end] / (phys.k_B * atmos.tmp[end])  # gas contribution
                N_t += N_g  # add atoms in this gas to total atoms
            end

            # Write elemental abundances
            open(elempath,"w") do f
                write(f,"# Elemental abundances derived from AGNI volatiles \n")
                for (i,e) in enumerate(phys.elems_standard)
                    if N_t[i] > 1.0e-30
                        # skip this element if its abundance is too small
                        # normalise relative to hydrogen
                        write(f, @sprintf("%s    %.3f \n",e,log10(N_t[i]/N_t[1]) + 12.0))
                        count_elem_nonzero += 1
                    end
                end
            end
        end

        # Write PT profile
        open(joinpath(atmos.fastchem_work,"pt.dat"),"w") do f
            write(f,"# AGNI temperature structure \n")
            write(f,"# bar, kelvin \n")
            for i in 1:atmos.nlev_c
                write(  f,
                        @sprintf("%.6e    %.6e \n",
                            atmos.p[i]*1e-5,
                            max(atmos.fastchem_floor,atmos.tmp[i])
                            )
                     )
            end
        end

        # Run fastchem
        run(pipeline(`$execpath $confpath`, stdout=devnull))

        # Check monitor output
        monitorpath::String = joinpath(atmos.fastchem_work,"monitor.dat")
        data = readdlm(monitorpath, '\t', String)
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
        if !isfile(chempath)
            @error "Could not find fastchem output"
            return 1
        end
        (data,head) = readdlm(chempath, '\t', Float64, header=true)
        data = transpose(data)  # convert to: gas, level

        # Clear VMRs
        for g in atmos.gas_names
            fill!(atmos.gas_vmr[g], 0.0)
        end

        # Parse gas chemistry
        g_fc::String = "_unset"
        d_fc::Dict = Dict{String, Int}()
        g_in::String = "_unset"
        match::Bool = false
        N_t = data[4,:] # at each level: sum of gas number densities

        for (i,h) in enumerate(head)  # for each column (gas)

            # skip T and P
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

        # Find where we truncated the temperature profile,
        #      and make sure that regions above that use the same gas_vmr values
        for i in range(start=atmos.nlev_c, stop=1, step=-1)
            if atmos.tmp[i] < atmos.fastchem_floor
                for g in atmos.gas_names
                    atmos.gas_vmr[g][1:i] .= atmos.gas_vmr[g][i+1]
                end
                break
            end
        end

        # See docstring for return codes
        return state
    end


end # end module
