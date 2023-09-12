# Contains functions for setting-up the temperature profile analytically

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module setpt

    import phys
    import atmosphere

    using PCHIPInterpolation

    # Read atmosphere T(p) from a CSV file (overwites existing P_boa and P_toa values).
    function fromcsv!(atmos::atmosphere.Atmos_t, fpath::String)

        # Check file exists
        if !isfile(fpath)
            error("The file '$fpath' does not exist")
        end 

        # Read file 
        content = readlines(fpath)

        # Parse file once, to get number of level edges 
        nlev_l = 0  
        for l in content
            if isempty(l) || (l[1] == '#') || (l[1] == '\n')
                continue 
            end
            nlev_l += 1
        end

        # Validate 
        if nlev_l < 3
            error("Csv file contains too few levels (contains $nlev_l edge values)")
        end
        nlev_c = nlev_l - 1

        # Allocate temporary T and P arrays 
        tmpl = zeros(Float64,nlev_l)
        pl   = zeros(Float64,nlev_l)

        # Parse file again, storing data this time
        idx = 1
        for l in content
            # Skip
            if isempty(l) || (l[1] == '#') || (l[1] == '\n')
                continue 
            end

            # Read
            lsplit = split(strip(l,[' ','#','\t','\n']),',')
            p_val = parse(Float64,lsplit[1])  # Pressure [Pa]
            t_val = parse(Float64,lsplit[2])  # Temperature [K]

            # Validate 
            if p_val <= 0.0
                error("Negative pressure(s) in csv file")
            end 
            if t_val <= 0.0
                error("Negative temperature(s) in csv file")
            end 

            # Store
            pl[idx] = p_val 
            tmpl[idx] = t_val

            # Iterate 
            idx += 1
        end 

        # Check if arrays are flipped 
        if pl[1] > pl[2]  
            pl = reverse(pl)
            tmpl = reverse(tmpl)
        end 

        # Check that pressure is monotonic
        for i in 1:nlev_c
            if pl[i] > pl[i+1]
                error("Pressure is not monotonic in csv file")
            end
        end

        # Set pressure grid
        atmos.p_boa = pl[end]
        atmos.p_toa = pl[1]
        atmosphere.generate_pgrid!(atmos)

        # Interpolate from the provided grid to the generated one
        itp = Interpolator(pl, tmpl) # Cell edges 
        atmos.tmpl[:] .= itp.(atmos.pl[:])

        itp = Interpolator(pl, tmpl) # Cell centres 
        atmos.tmp[:] .= itp.(atmos.p[:])
    end

    # Set atmosphere to be isothermal
    function isothermal!(atmos::atmosphere.Atmos_t, set_tmp::Float64)
        if !atmos.is_param
            error("Atmosphere parameters not set")
        end 
        atmos.tmpl[:] .= set_tmp 
        atmos.tmp[:] .= set_tmp
    end 

    # Set atmosphere to dry adiabat
    function dry_adiabat!(atmos::atmosphere.Atmos_t)
        # Validate input
        if !(atmos.is_alloc && atmos.is_param) 
            error("Atmosphere is not setup")
        end 

        # Calculate cell-centre values
        for idx in 1:atmos.nlev_c
            atmos.tmp[idx] = atmos.tmpl[end] * ( atmos.p[idx] / atmos.pl[end] ) ^ ( phys.R_gas / atmos.layer_cp[idx] )
        end
        
        # Calculate cell-edge values
        for idx in 2:atmos.nlev_l-1
            cp = 0.5 * ( atmos.layer_cp[idx-1] + atmos.layer_cp[idx])
            atmos.tmpl[idx] = atmos.tmpl[end] * ( atmos.pl[idx] / atmos.pl[end] ) ^ ( phys.R_gas / cp )
        end

        # Calculate top boundary
        dt = atmos.tmp[1]-atmos.tmpl[2]
        dp = atmos.p[1]-atmos.pl[2]
        atmos.tmpl[1] = atmos.tmp[1] + dt/dp * (atmos.pl[1] - atmos.p[1])
    end 

    # Set atmosphere to have an isothermal stratosphere
    function stratosphere!(atmos::atmosphere.Atmos_t, strat_tmp::Float64)

        # Keep stratosphere below tstar value
        strat_tmp = min(strat_tmp, atmos.tstar)

        # Loop upwards from bottom of model
        strat = false
        for i in range(atmos.nlev_c,1,step=-1)
            # Find tropopause 
            strat = strat || (atmos.tmp[i] < strat_tmp) || (atmos.tmpl[i+1] < strat_tmp) 
            
            # Apply stratosphere to this level if required
            if strat
                atmos.tmp[i]    = strat_tmp 
                atmos.tmpl[i+1] = strat_tmp 
            end
        end

        # Handle topmost level 
        if strat
            atmos.tmpl[1] = strat_tmp 
        end

    end

    # Ensure that the surface isn't supersaturated
    function prevent_surfsupersat!(atmos::atmosphere.Atmos_t)
        if !(atmos.is_alloc && atmos.is_param) 
            error("Atmosphere is not setup")
        end 

        # For each condensible volatile
        for (i_gas, gas) in enumerate(atmos.gases)
            # Get mixing ratio at surface
            mr = atmosphere.get_mr(atmos, gas, lvl=atmos.nlev_c)
            if mr < 1.0e-10
                continue 
            end

            # Check surface pressure (should not be supersaturated)
            psat = phys.calc_Psat(gas, atmos.tmpl[end])
            Tsat = phys.calc_Tdew(gas, atmos.p_boa)
            if (atmos.tmpl[end] < Tsat) && (atmos.tmpl[end] < phys.lookup_T_crit[gas])
                # Reduce amount of volatile until reaches phase curve 
                atmos.p_boa = atmos.pl[end]*(1.0-mr) + psat
                atmos.layer_mr[atmos.nlev_c, i_gas] = psat / atmos.p_boa  

                # Renormalise mixing ratios
                norm_factor = sum(atmos.layer_mr[atmos.nlev_c, :])
                for j_gas in 1:length(atmos.gases)
                    atmos.layer_mr[atmos.nlev_c, j_gas] /= norm_factor
                end

                # Generate new pressure grid 
                atmosphere.generate_pgrid!(atmos)

                # Inform user
                println("Atmosphere: Supersaturated surface - setting partial pressure to saturation pressure")
            end
        end 
    end


    # Set atmosphere to phase curve of 'con' when it enters the condensible region
    function condensing!(atmos::atmosphere.Atmos_t, gas::String)

        if !(atmos.is_alloc && atmos.is_param) 
            error("Atmosphere is not setup")
        end 
        if !(gas in keys(phys.lookup_L_vap))
            error("Invalid condensible $gas")
        end

        # Check if each level is condensing. If it is, place on phase curve
        for i in 1:atmos.nlev_c
            # Cell centre
            Tsat = phys.calc_Tdew(gas,atmos.p[i])
            if atmos.tmp[i] < Tsat
                atmos.tmp[i] = Tsat
            end
                
            # Cell edge
            Tsat = phys.calc_Tdew(gas,atmos.pl[i])
            if atmos.tmpl[i] < Tsat
                atmos.tmpl[i] = Tsat
            end
        end
    end

end 
