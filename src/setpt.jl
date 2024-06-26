# Contains functions for setting-up the temperature profile analytically

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module setpt

    import ..phys
    import ..atmosphere

    using PCHIPInterpolation
    using NCDatasets

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

        # Extrapolate loaded grid to lower pressures (prevent domain error)
        if (atmos.pl[1] < pl[1])
            p_ext = atmos.pl[1]/1.1
            t_ext = (tmpl[1] - tmpl[3])/(pl[1] - pl[3]) * (p_ext - pl[1])
            pushfirst!(pl,   p_ext)
            pushfirst!(tmpl, t_ext)
        end

        # Extrapolate loaded grid to higher pressures 
        if (atmos.pl[end] > pl[end])
            p_ext = atmos.pl[end]*1.1
            t_ext = (tmpl[end] - tmpl[end-2])/(pl[end] - pl[end-2]) * (p_ext - pl[end])
            push!(pl,   p_ext)
            push!(tmpl, t_ext)
        end

        # Interpolate from the loaded grid to the required one
        #   This uses log-pressures in order to make the interpolation behave 
        #   reasonably across the entire grid.
        itp = Interpolator(log10.(pl), tmpl) 
        atmos.tmpl[:] .= itp.(log10.(atmos.pl[:]))  # Cell edges 
        atmos.tmp[:]  .= itp.(log10.(atmos.p[:]))   # Cell centres 

        atmosphere.calc_layer_props!(atmos)
        return nothing
    end

    """
    Load atmosphere data from NetCDF file (must have same number of levels)
    """
    function fromncdf!(atmos::atmosphere.Atmos_t, fpath::String)

        # Check file exists
        if !isfile(fpath)
            error("The file '$fpath' does not exist")
        end 

        # Open file 
        fpath = abspath(fpath)
        ds = Dataset(fpath,"r")

        # Allocate interleaved temperature profile 
        nlev_c::Int = length(ds["p"][:])
        arr_n::Int = nlev_c + nlev_c + 1
        arr_T::Array{Float64, 1} = zeros(Float64, arr_n)
        arr_P::Array{Float64, 1} = zeros(Float64, arr_n)

        # top
        arr_T[1] = ds["tmpl"][1]
        arr_P[1] = min(ds["pl"][1],atmos.pl[1])         # extend to lower pressures if required

        # middle 
        idx::Int = 0
        for i in 1:nlev_c
            idx = (i-1)*2
            arr_T[idx+1] = ds["tmpl"][i]
            arr_T[idx+2] = ds["tmp"][i]
            arr_P[idx+1] = ds["pl"][i]
            arr_P[idx+2] = ds["p"][i]
        end 

        # bottom
        arr_T[end] = ds["tmpl"][end]
        arr_P[end] = max(ds["pl"][end], atmos.pl[end])      # extend to higher pressures if required 

        # interpolate 
        itp = Interpolator(log10.(arr_P), arr_T) 
        atmos.tmpl[:] .= itp.(log10.(atmos.pl[:]))  # Cell edges 
        atmos.tmp[:]  .= itp.(log10.(atmos.p[:]))   # Cell centres 

        # properties 
        atmos.tmp_surf = ds["tmp_surf"][1]

        # Close file 
        close(ds)

        atmosphere.calc_layer_props!(atmos)
        return nothing
    end # end load_ncdf

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

    # Set atmosphere to be isothermal at the given temperature
    function add!(atmos::atmosphere.Atmos_t, delta::Float64)
        if !atmos.is_param
            error("Atmosphere parameters not set")
        end 
        atmos.tmpl[:] .+= delta 
        atmos.tmp[:]  .+= delta

        atmosphere.calc_layer_props!(atmos)
        return nothing
    end 

    # Set atmosphere to dry adiabat
    function dry_adiabat!(atmos::atmosphere.Atmos_t)
        # Validate input
        if !(atmos.is_alloc && atmos.is_param) 
            error("Atmosphere is not setup or allocated")
        end 

        # Set surface 
        atmos.tmpl[end] = atmos.tmp_surf

        # Thermodynamics etc 
        atmosphere.calc_layer_props!(atmos)

        # Lapse rate dT/dp 
        grad::Float64 = 0.0

        # Calculate values 
        for i in range(start=atmos.nlev_c, stop=1, step=-1)

            # Set cp based on temperature at the level below this one
            tmp_eval = atmos.tmp_surf 
            if i < atmos.nlev_c 
                tmp_eval = atmos.tmp[i+1]
            end
            atmos.layer_cp[i] = 0.0
            for gas in atmos.gas_names
                atmos.layer_cp[i] += atmos.gas_vmr[gas][i] * atmos.gas_dat[gas].mmw * phys.get_Cp(atmos.gas_dat[gas], tmp_eval) / atmos.layer_mmw[i]
            end

            # Cell-edge to cell-centre 
            grad = phys.R_gas * atmos.tmpl[i+1] / (atmos.pl[i+1] * atmos.layer_mmw[i] * atmos.layer_cp[i])
            atmos.tmp[i] = atmos.tmpl[i+1] + grad * (atmos.p[i]-atmos.pl[i+1])
            atmos.tmp[i] = max(atmos.tmp[i], atmos.tmp_floor)

            # Cell-centre to cell-edge 
            grad = phys.R_gas * atmos.tmp[i] / (atmos.p[i] * atmos.layer_mmw[i] * atmos.layer_cp[i])
            atmos.tmpl[i] = atmos.tmp[i] + grad * (atmos.pl[i]-atmos.p[i])
            atmos.tmpl[i] = max(atmos.tmpl[i], atmos.tmp_floor)
        end 

        # Thermodynamics at new temperature profile 
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

        x::Float64 = 0.0
        psat::Float64 = 0.0

        # For each condensable volatile
        for gas in atmos.gas_names
            # Get mole fraction at surface
            x = atmos.gas_vmr[gas][atmos.nlev_c]
            if x < 1.0e-10
                continue 
            end

            # Check criticality 
            if (atmos.tmp_surf > atmos.gas_dat[gas].T_crit)
                continue 
            end 

            # Check surface pressure (should not be supersaturated)
            psat = phys.get_Psat(atmos.gas_dat[gas], atmos.tmp_surf)
            if x*atmos.pl[end] > psat
                # Reduce amount of volatile until reaches phase curve 
                atmos.p_boa = atmos.pl[end]*(1.0-x) + psat
                atmos.gas_vmr[gas][atmos.nlev_c] = psat / atmos.p_boa  

                # Check that p_boa is still reasonable 
                if atmos.p_boa <= 10.0 * atmos.p_toa 
                    error("Supersaturation check ($gas) resulted in an unreasonably small surface pressure")
                end 
            end
        end 

        # Renormalise mole fractions
        tot_vmr = 0.0
        for g in atmos.gas_names
            tot_vmr += atmos.gas_vmr[g][atmos.nlev_c]
        end 
        for g in atmos.gas_names
            atmos.gas_vmr[g][atmos.nlev_c] /= tot_vmr
        end 

        # Generate new pressure grid 
        atmosphere.generate_pgrid!(atmos)
        atmosphere.calc_layer_props!(atmos)
        
        return nothing
    end


    """
    Set T = max(T,T_dew) for a specified gas.

    Does not modify VMRs or surface temperature.
    """
    function saturation!(atmos::atmosphere.Atmos_t, gas::String)

        if !(atmos.is_alloc && atmos.is_param) 
            error("Atmosphere is not setup or allocated")
        end 

        # gas is present?
        if !(gas in atmos.gas_names)
            return nothing
        end

        x::Float64 = 0.0
        Tdew::Float64 = 0.0

        # Check if each level is condensing. If it is, place on phase curve
        for i in 1:atmos.nlev_c

            x = atmos.gas_vmr[gas][i]
            if x < 1.0e-10 
                continue
            end

            # Set cell-centre temperatures
            Tdew = phys.get_Tdew(atmos.gas_dat[gas], atmos.p[i])
            if atmos.tmp[i] < Tdew 
                atmos.tmp[i] = Tdew 
                atmos.gas_sat[gas][i] = true
            end 
        end
        
        # Set cell-edge temperatures
        atmosphere.set_tmpl_from_tmp!(atmos)

        # calculate properties 
        atmosphere.calc_layer_props!(atmos)

        return nothing
    end

end # end module
