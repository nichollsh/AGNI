# Contains functions for setting-up the temperature profile analytically

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module setpt

    import ..phys
    import ..atmosphere
    import ..chemistry

    include("guillot.jl")
    import .guillot

    using NCDatasets
    using Printf
    using LoggingExtras
    import Interpolations: interpolate, Gridded, Linear, Flat, extrapolate, Extrapolation

    # Process a series of requests describing T(p)
    function request!(atmos::atmosphere.Atmos_t, request::Array{String,1})::Bool
        num_req::Int = length(request)          # Number of requests
        idx_req::Int = 1                        # Index of current request
        str_req::String = atmosphere.UNSET_STR  # String of current request
        prt_req::String = "Setting T(p): "
        while idx_req <= num_req
            # get command
            str_req = strip(lowercase(request[idx_req]))
            prt_req *= str_req*", "

            # handle requests
            if str_req == "dry"
                # dry adiabat from surface
                setpt.dry_adiabat!(atmos)

            elseif str_req == "str"
                # isothermal stratosphere
                idx_req += 1
                setpt.stratosphere!(atmos, parse(Float64, request[idx_req]))

            elseif str_req == "loglin"
                # log-linear profile between T_surf and T_top
                idx_req += 1
                setpt.loglinear!(atmos, parse(Float64, request[idx_req]))

            elseif str_req == "iso"
                # isothermal profile
                idx_req += 1
                setpt.isothermal!(atmos, parse(Float64, request[idx_req]))

            elseif str_req == "csv"
                # set from csv file
                idx_req += 1
                setpt.fromcsv!(atmos,request[idx_req])

            elseif str_req == "ncdf"
                # set from NetCDF file
                idx_req += 1
                setpt.fromncdf!(atmos,request[idx_req])

            elseif str_req == "add"
                # add X kelvin from the currently stored T(p)
                idx_req += 1
                setpt.add!(atmos,parse(Float64, request[idx_req]))

            elseif str_req == "surfsat"
                # ensure surface is not super-saturated
                chemistry.restore_composition!(atmos)
                chemistry._sat_surf!(atmos)

            elseif str_req == "sat"
                # condensing a volatile
                idx_req += 1
                setpt.saturation!(atmos, request[idx_req])

            elseif str_req == "ana"
                # analytic solution
                setpt.analytic!(atmos)

            else
                @error "Invalid initial state '$str_req'"
                return false
            end

            atmosphere.calc_layer_props!(atmos)

            # iterate
            idx_req += 1
        end
        @info prt_req[1:end-2]
        return true
    end

    # Set by interpolating from T,P arrays
    function fromarrays!(atmos::atmosphere.Atmos_t, pl::Array, tmpl::Array)

        if !atmos.is_param
            @error "setpt: Atmosphere parameters not set"
            return
        end

        # Check if arrays are flipped
        #      Pressure must be increasing with index
        if pl[1] > pl[2]
            pl = reverse(pl)
            tmpl = reverse(tmpl)
        end

        # Check that pressure is monotonic
        for i in 1:length(pl)-1
            if pl[i] > pl[i+1]
                @error "setpt: input array of pressures is not monotonically increasing"
                return
            end
        end

        # Extrapolate loaded grid to lower pressures (prevent domain error)
        if (atmos.pl[1] < pl[1])
            pushfirst!(pl,   atmos.pl[1]/1.1)
            pushfirst!(tmpl, tmpl[1])
        end

        # Extrapolate loaded grid to higher pressures
        if (atmos.pl[end] > pl[end])
            push!(pl,   atmos.pl[end]*1.1)
            push!(tmpl, tmpl[end])
        end

        # Interpolate from the loaded grid to the required one
        #   This uses log-pressures in order to make the interpolation behave
        #   reasonably across the entire grid.
        itp::Extrapolation = extrapolate(interpolate((log10.(pl),),tmpl, Gridded(Linear())), Flat())
        @. atmos.tmpl = itp(log10(atmos.pl))  # Cell edges
        @. atmos.tmp  = itp(log10(atmos.p ))   # Cell centres

        return
    end

    # Read atmosphere T(p) from a CSV file (does not overwrite p_boa and p_toa)
    function fromcsv!(atmos::atmosphere.Atmos_t, fpath::String)

        # Check file exists
        if !isfile(fpath)
            @error "setpt: file '$fpath' does not exist"
            return
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
            @error "setpt: file contains too few levels (contains $nlev_l edge values)"
            return
        end

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
                @error "setpt: Negative pressure(s) in csv file"
                return
            end
            if t_val <= 0.0
                @error "setpt: Negative temperature(s) in csv file"
                return
            end

            # Store
            pl[idx] = p_val
            tmpl[idx] = t_val

            # Iterate
            idx += 1
        end

        # Use fromarrays function to do the rest
        fromarrays!(atmos, pl, tmpl)

        return
    end

    """
    Load atmosphere data from NetCDF file (must have same number of levels)
    """
    function fromncdf!(atmos::atmosphere.Atmos_t, fpath::String)

        # Check file exists
        if !isfile(fpath)
            @error "setpt: The file '$fpath' does not exist"
            return
        end

        # Open file
        fpath = abspath(fpath)
        @debug "Setting PT from NetCDF file "
        @debug "ALL DEBUG SUPPRESSED"
        with_logger(MinLevelLogger(current_logger(), Logging.Info-200)) do
            ds = Dataset(fpath,"r")

            # Allocate interleaved temperature profile
            nlev_c::Int = length(ds["p"][:])
            arr_n::Int = nlev_c + nlev_c + 1
            arr_T::Array{Float64, 1} = zeros(Float64, arr_n)
            arr_P::Array{Float64, 1} = zeros(Float64, arr_n)

            # top
            arr_T[1] = ds["tmpl"][1]
            arr_P[1] = ds["pl"][1]

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
            arr_P[end] = ds["pl"][end]

            # extend with constant values to avoid issues with interpolator
            newtop = min(atmos.pl[1], arr_P[1])/2.0
            pushfirst!(arr_P,newtop)
            pushfirst!(arr_T, ds["tmpl"][1])

            newbot = max(atmos.pl[end], arr_P[end])*2.0
            push!(arr_P, newbot)
            push!(arr_T, ds["tmpl"][end])

            # properties
            atmos.tmp_surf = ds["tmp_surf"][1]

            # Close file
            close(ds)

            # interpolate
            itp = extrapolate(interpolate((log10.(arr_P),),arr_T,Gridded(Linear())),Flat())
            @. atmos.tmpl = itp(log10(atmos.pl))  # Cell edges
            @. atmos.tmp  = itp(log10(atmos.p ))   # Cell centres

        end
        @debug "ALL DEBUG RESTORED"

        return
    end # end load_ncdf

    # Set atmosphere to be isothermal at the given temperature
    function isothermal!(atmos::atmosphere.Atmos_t, set_tmp::Float64)
        if !atmos.is_param
            @error "setpt: Atmosphere parameters not set"
            return
        end
        fill!(atmos.tmpl, set_tmp)
        fill!(atmos.tmp , set_tmp)

        return
    end

    # Set atmosphere to be isothermal at the given temperature
    function add!(atmos::atmosphere.Atmos_t, delta::Float64)
        if !atmos.is_param
            @error "setpt: Atmosphere parameters not set"
            return
        end
        @. atmos.tmpl += delta
        @. atmos.tmp  += delta

        return nothing
    end

    # Set atmosphere to dry adiabat
    function dry_adiabat!(atmos::atmosphere.Atmos_t)
        # Validate input
        if !(atmos.is_alloc && atmos.is_param)
            @error "setpt: Atmosphere is not setup or allocated"
            return
        end

        # Set surface
        atmos.tmpl[end] = atmos.tmp_surf

        # Set mmw
        atmosphere.calc_profile_mmw!(atmos)

        # Lapse rate dT/dp
        grad::Float64 = 0.0

        # Calculate values
        for i in range(start=atmos.nlev_c, stop=1, step=-1)

            # Set cp and rho based on temperature of the level below this one
            atmos.tmp[i] = atmos.tmp_surf
            if i < atmos.nlev_c
                atmos.tmp[i] = atmos.tmp[i+1]
            end
            atmosphere.calc_single_cpkc!(atmos, i)
            atmosphere.calc_single_density!(atmos, i)

            # Evaluate lapse rate dT/dp
            grad = 1 / (atmos.layer_ρ[i] * atmos.layer_cp[i])

            # Cell-edge to cell-centre
            atmos.tmp[i] = atmos.tmpl[i+1] + grad * (atmos.p[i]-atmos.pl[i+1])
            atmos.tmp[i] = max(atmos.tmp[i], atmos.tmp_floor)

            # Cell-centre to cell-edge
            atmos.tmpl[i] = atmos.tmp[i] + grad * (atmos.pl[i]-atmos.p[i])
            atmos.tmpl[i] = max(atmos.tmpl[i], atmos.tmp_floor)
        end


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

        return nothing
    end

    # Set atmosphere to have a log-linear T(p) profile
    function loglinear!(atmos::atmosphere.Atmos_t, top_tmp::Float64)

        # Keep top_tmp below tmp_surf value
        top_tmp = min(top_tmp, atmos.tmp_surf)

        # Set surface and near-surface
        atmos.tmpl[end] = atmos.tmp_surf

        # Loop upwards from bottom of atmosphere
        dTdP::Float64 = (top_tmp - atmos.tmp_surf)/log10(atmos.pl[end]/atmos.pl[1])
        for i in range(atmos.nlev_l-1,1,step=-1)
            atmos.tmpl[i] = atmos.tmpl[i+1] + dTdP * log10(atmos.pl[i+1]/atmos.pl[i])
        end

        # Set cell-centres
        atmos.tmp[1:end] .= 0.5 .* (atmos.tmpl[1:end-1] + atmos.tmpl[2:end])

        return nothing
    end




    """
    **Set T = max(T,T_dew) for a specified gas.**

    Does not modify VMRs. Does update water-cloud locations.
    """
    function saturation!(atmos::atmosphere.Atmos_t, gas::String)

        if !(atmos.is_alloc && atmos.is_param)
            @error "setpt: Atmosphere is not setup or allocated"
        end

        # gas is present?
        if !(gas in atmos.gas_names)
            return nothing
        end

        xgas::Float64 = 0.0

        # Return the new temperature (Float) and whether saturated (Bool)
        function _tdew(tmp::Float64,pgas::Float64)::Tuple
            # Tiny partial pressure
            if pgas < 1e-99
                return (tmp, false)
            end

            # Supercritical
            if tmp >= atmos.gas_dat[gas].T_crit
                return (tmp, false)
            end

            # Otherwise...
            Tdew = phys.get_Tdew(atmos.gas_dat[gas], pgas)
            if tmp <= Tdew + 1e-2
                return (Tdew, true)
            else
                return (tmp, false)
            end
        end

        # Check if each level is saturated: tmp < Tdew. If it is, set to Tdew.
        for i in range(start=atmos.nlev_l, stop=1, step=-1)

            # Composition
            xgas = atmos.gas_vmr[gas][min(i,atmos.nlev_c)]

            # Cell-centres
            if i <= atmos.nlev_c
                (atmos.tmp[i], atmos.gas_sat[gas][i]) = _tdew(atmos.tmp[i], atmos.p[i]*xgas)
            end

            # Cell-edges
            atmos.tmpl[i] = _tdew(atmos.tmpl[i], atmos.pl[i]*xgas)[1]
        end
        atmos.tmp_surf = atmos.tmpl[end]

        # Set cloud
        if (gas == "H2O") && atmos.control.l_cloud
            fill!(atmos.cloud_arr_r, 0.0)
            fill!(atmos.cloud_arr_l, 0.0)
            fill!(atmos.cloud_arr_f, 0.0)

            atmos.cloud_arr_r[atmos.gas_sat["H2O"][:]] .= atmos.cloud_val_r
            atmos.cloud_arr_l[atmos.gas_sat["H2O"][:]] .= atmos.cloud_val_l
            atmos.cloud_arr_f[atmos.gas_sat["H2O"][:]] .= atmos.cloud_val_f
        end

        return nothing
    end

    """
    **Set temperature profile using Guillot (2010) analytic solution.**
    """
    function analytic!(atmos::atmosphere.Atmos_t)

        # Evalulate Tirr from the instellation
        Tirr = (atmos.instellation / phys.σSB)^0.25

        # Evalulate Tint from flux
        Tint = (atmos.flux_int / phys.σSB)^0.25

        # Evalulate cell-centre temperatures
        for i in 1:atmos.nlev_c
            # get LW optical depth
            τ = guillot.eval_tau(atmos.p[i], atmos.layer_grav[i])

            # set temperature
            atmos.tmp[i] = guillot.eval_T4_cos(τ, Tint, Tirr, atmos.zenith_degrees)^0.25
        end

        # Set cell-edge temperatures
        atmosphere.set_tmpl_from_tmp!(atmos)
        return nothing
    end

end # end module
