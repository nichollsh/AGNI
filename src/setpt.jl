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

    # Read atmosphere T(p) from a CSV file (does not overwrite p_boa and p_toa)
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

        # Extend loaded grid isothermally to lower pressures (prevent domain error)
        if (atmos.pl[1] < pl[1])
            pushfirst!(pl,   atmos.pl[1]/1.1)
            pushfirst!(tmpl, tmpl[1]  )
        end

        # Extend loaded grid isothermally to higher pressures 
        if (atmos.pl[end] > pl[end])
            push!(pl,   atmos.pl[end]*1.1)
            push!(tmpl, tmpl[end])
        end

        # Interpolate from the loaded grid to existing one
        itp = Interpolator(pl, tmpl) 
        atmos.tmpl[:] .= itp.(atmos.pl[:])  # Cell edges 
        atmos.tmp[:]  .= itp.(atmos.p[:])   # Cell centres 

        atmosphere.calc_layer_props!(atmos)
        return nothing
    end

    # Set atmosphere to be isothermal at the given temperature
    function isothermal!(atmos::atmosphere.Atmos_t, set_tmp::Float64)
        if !atmos.is_param
            error("Atmosphere parameters not set")
        end 
        atmos.tmpl[:] .= set_tmp 
        atmos.tmp[:]  .= set_tmp

        atmosphere.calc_layer_props!(atmos)
        return nothing
    end 

    # Set atmosphere to dry adiabat
    function dry_adiabat!(atmos::atmosphere.Atmos_t)
        # Validate input
        if !(atmos.is_alloc && atmos.is_param) 
            error("Atmosphere is not setup or allocated")
        end 

        # Calculate cell-centre values
        for idx in 1:atmos.nlev_c
            cp = atmos.layer_cp[idx] * atmos.layer_mmw[idx]
            atmos.tmp[idx] = atmos.tmpl[end] * ( atmos.p[idx] / atmos.pl[end] ) ^ ( phys.R_gas / cp )
        end
        
        # Calculate cell-edge values
        for idx in 2:atmos.nlev_l-1
            cp = 0.5 * ( atmos.layer_cp[idx-1] * atmos.layer_mmw[idx-1] + atmos.layer_cp[idx] * atmos.layer_mmw[idx])
            atmos.tmpl[idx] = atmos.tmpl[end] * ( atmos.pl[idx] / atmos.pl[end] ) ^ ( phys.R_gas / cp )
        end

        # Calculate top boundary
        dt = atmos.tmp[1]-atmos.tmpl[2]
        dp = atmos.p[1]-atmos.pl[2]
        atmos.tmpl[1] = atmos.tmp[1] + dt/dp * (atmos.pl[1] - atmos.p[1])

        atmosphere.calc_layer_props!(atmos)
        return nothing
    end 

    # Set atmosphere to have an isothermal stratosphere
    function stratosphere!(atmos::atmosphere.Atmos_t, strat_tmp::Float64)

        # Keep stratosphere below tmp_surf value
        strat_tmp = min(strat_tmp, atmos.tmp_surf)

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

        atmosphere.calc_layer_props!(atmos)
        return nothing
    end

    # Set atmosphere to have a log-linear T(p) profile
    function loglinear!(atmos::atmosphere.Atmos_t, top_tmp::Float64)

        # Keep top_tmp below tmp_surf value
        top_tmp = min(top_tmp, atmos.tmp_surf)

        # Set surface and near-surface
        atmos.tmpl[end] = atmos.tmp_surf
        atmos.tmpl[end-1] = atmos.tmp_surf 

        # Loop upwards from bottom of model, assuming temperatures are log-spaced
        dtdi = (top_tmp - atmos.tmp_surf)/(atmos.nlev_l-1)
        for i in range(atmos.nlev_l-2,1,step=-1)
            atmos.tmpl[i] = atmos.tmpl[i+1] + dtdi
        end

        # Set cell-centres 
        atmos.tmp[1:end] .= 0.5 .* (atmos.tmpl[1:end-1] + atmos.tmpl[2:end])

        atmosphere.calc_layer_props!(atmos)
        return nothing
    end


    # Ensure that the surface isn't supersaturated
    function prevent_surfsupersat!(atmos::atmosphere.Atmos_t)
        if !(atmos.is_alloc && atmos.is_param) 
            error("Atmosphere is not setup or allocated")
        end 

        # For each condensible volatile
        for (i_gas, gas) in enumerate(atmos.gases)
            # Get mole fraction at surface
            x = atmosphere.get_x(atmos, gas, atmos.nlev_c)
            if x < 1.0e-10
                continue 
            end

            # Check surface pressure (should not be supersaturated)
            psat = phys.calc_Psat(gas, atmos.tmpl[end])
            Tsat = phys.calc_Tdew(gas, atmos.p_boa*x)
            if (atmos.tmpl[end] < Tsat) && (atmos.tmpl[end] < phys.lookup_safe("T_crit",gas))
                # Reduce amount of volatile until reaches phase curve 
                atmos.p_boa = atmos.pl[end]*(1.0-x) + psat
                atmos.layer_x[atmos.nlev_c, i_gas] = psat / atmos.p_boa  

                # Check that p_boa is still reasonable 
                if atmos.p_boa <= 10.0 * atmos.p_toa 
                    error("Supersaturation check ($gas) resulted in an unreasonably small surface pressure")
                end 

                # Renormalise mole fractions
                norm_factor = sum(atmos.layer_x[atmos.nlev_c, :])
                for j_gas in 1:length(atmos.gases)
                    atmos.layer_x[atmos.nlev_c, j_gas] /= norm_factor
                end

                # Generate new pressure grid 
                atmosphere.generate_pgrid!(atmos)
            end
        end 

        atmosphere.calc_layer_props!(atmos)
        return nothing
    end


    # Set atmosphere to phase curve of gas when it enters the condensible region (does not modify tmp_surf)
    function condensing!(atmos::atmosphere.Atmos_t, gas::String)

        if !(atmos.is_alloc && atmos.is_param) 
            error("Atmosphere is not setup or allocated")
        end 

        # gas is present?
        if !(gas in atmos.gases)
            return nothing
        end

        # gas has thermodynamics setup?
        if phys.lookup_safe("L_vap",gas) < 1e-20
            return nothing
        end

        # apply condensation curve 
        atmosphere.apply_vlcc!(atmos, gas)

        # apply cloud 
        atmosphere.water_cloud!(atmos)

        # calculate properties 
        atmosphere.calc_layer_props!(atmos)

        return nothing
    end

end # end module
