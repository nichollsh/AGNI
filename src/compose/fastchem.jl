module fastchem

    # Import packages
    using Printf
    using Logging
    import DelimitedFiles:readdlm

    # Include local modules
    using ..atmosphere
    using ..phys

    # Constants
    const SMOOTH_SCALE::Float64 = 12.0      # smoothing scale for fastchem floor

    """
    **Smoothly transform temperature around floor**

    Approximately returns _x for _x > _floor.
    Approximately returns _floor for _x < _floor.
    Transitions discontinuously for a->Inf
    """
    function _transform_floor(_x::Float64, _floor::Float64)::Float64
        d = 1.0 / ( 1.0 + exp(-SMOOTH_SCALE*(_x-_floor)/_floor) )
        return _x*d + _floor*(1-d)
    end

    """
    **Run FastChem chemistry solver**

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `write_cfg::Bool`                 write config and elements

    Returns:
    - `state::Int64`                    fastchem state (0: success, 1: critical_fail, 2: elem_fail, 3: conv_fail, 4: both_fail)
    """
    function run_fastchem!(atmos::atmosphere.Atmos_t, write_cfg::Bool)::Int64

        state::Int64 = 0
        count_elem_nonzero::Int64 = 0

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
                N_inp_t = zeros(Float64, length(phys.elems_standard)) # total atoms in all gases
                N_inp_g = zeros(Float64, length(phys.elems_standard)) # atoms in current gas
                #    loop over gases
                for gas in atmos.gas_names
                    fill!(N_inp_g, 0.0)

                    # non-zero abundance?
                    if atmos.gas_vmr[gas][end] < SMALL_FLOAT
                        continue
                    end

                    # count atoms in this gas
                    d = phys.count_atoms(gas)
                    for (i,e) in enumerate(phys.elems_standard)
                        if haskey(d, e)
                            N_inp_g[i] += d[e] # N_inp_g stores num of atoms in this gas
                        end
                    end

                    # Get gas abundance from original VMR value
                    #    scale number of atoms by the abundance of the gas (p = Ng kB T)
                    N_inp_g *= atmos.gas_vmr[gas][end] * atmos.p[end] / (phys.k_B * atmos.tmp[end])

                    # Add atoms from this gas to total atoms in the mixture
                    N_inp_t += N_inp_g
                end

                # Check that we have some hydrogen...
                if N_inp_t[1] < SMALL_FLOAT
                    @warn "Cannot calculate metallicity of hydrogen-free mixture!"
                    state = 1
                end

                # Convert elemental abundances to metallicity number ratios, rel to hydrogen
                for (i,e) in enumerate(phys.elems_standard)
                    atmos.metal_calc[e] = N_inp_t[i]/N_inp_t[1]
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
                    if isfinite(atmos.metal_calc[e])
                        write(f, @sprintf("%s    %.3f \n",e,log10(atmos.metal_calc[e]) + 12.0))
                    else
                        @warn "Got non-finite metallicity for $e - adopting solar value"
                        write(f, @sprintf("%s    %.3f \n",e,phys.consts._solar_metallicity[e]))
                    end

                    count_elem_nonzero += 1
                end
            end

        end # end write metallicities

        # Work out which indices are visited by fc
        fc_levels::Array{Int,1} = Float64[atmos.nlev_c]
        if !atmos.fastchem_wellmixed
            fc_levels = collect(Float64,1:atmos.nlev_c)
        end

        # Write PT profile every time
        open(atmos.fastchem_prof,"w") do f
            write(f,"# AGNI temperature structure \n")
            write(f,"# bar, kelvin \n")

            for i in fc_levels
                write(f, @sprintf("%.6e    %.6e \n",
                            atmos.p[i]*1e-5,
                            _transform_floor(atmos.tmp[i], atmos.fastchem_floor) )
                        )
            end # /levels
        end # /file

        # Run fastchem
        run(pipeline(`$(atmos.fastchem_exec) $(atmos.fastchem_conf)`, stdout=devnull))

        # Check monitor output
        data = readdlm(atmos.fastchem_moni, '\t', String)
        fail_elem::String = ""
        fail_conv::String = ""
        for i in fc_levels
            if atmos.fastchem_wellmixed
                if startswith(data[2,6], 'f')
                    fail_elem *= @sprintf("%d ",i)
                end
                if startswith(data[2,5], 'f')
                    fail_conv *= @sprintf("%d ",i)
                end
            else
                if startswith(data[i+1,6],'f')
                    fail_elem *= @sprintf("%d ",i)
                end
                if startswith(data[i+1,5],'f')
                    fail_conv *= @sprintf("%d ",i)
                end
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
            @warn "Could not find FastChem output file '$(atmos.fastchem_chem)'"
            return 1
        end
        (data,head) = readdlm(atmos.fastchem_chem, '\t', Float64, header=true)
        data = transpose(data)  # convert to: gas, level

        # Clear VMRs now that surf metallicity has been recorded
        for g in atmos.gas_names
            fill!(atmos.gas_vmr[g],   0.0)
            fill!(atmos.gas_cvmr[g],  0.0)
        end

        # Parse gas chemistry
        g_fc::String = atmosphere.UNSET_STR  # gas name in fastchem
        d_fc::Dict = Dict{String, Int64}()     # ^ broken in to atoms
        g_in::String = atmosphere.UNSET_STR  # gas name in AGNI, matched by atom count
        match::Bool = false


        # for each column (gas)
        #   i = gas index in fc file
        #   h = gas name in header of fc file
        for (i,h) in enumerate(head)

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
                # convert number densities to VMR, and store
                if atmos.fastchem_wellmixed
                    # just 1 value
                    fill!(atmos.gas_vmr[g_in], data[i,1]/data[4,1])
                else
                    # whole profile (+= because of cis/trans being combined)
                    @. atmos.gas_vmr[g_in] += data[i,:] / data[4,:]
                end
            end # /match
        end # /gas

        return state
    end # /function

end
