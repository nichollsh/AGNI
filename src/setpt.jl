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

    """
    **Parse temperature value from keyword or numeric input.**

    Converts string keywords (`"teq"`, `"tsurf"`) or numeric values to Float64 temperature.

    Arguments:
    - `atmos::Atmos_t`              atmosphere struct instance.
    - `tmp::Union{String, Number}`  temperature (keyword string, numeric string, or number).

    Returns:
    - `Union{Float64, Nothing}` temperature [K], or `nothing` if invalid input.
    """
    function _parse_tmp_str(atmos::atmosphere.Atmos_t,
                                tmp::Union{String, Number})::Union{Float64, Nothing}

        if tmp isa String
            if lowercase(tmp) == "teq"
                return phys.calc_Teq(atmos.instellation, atmos.albedo_b)
            elseif lowercase(tmp) == "tsurf"
                return atmos.tmp_surf
            else
                return tryparse(Float64, tmp)
            end
        elseif tmp isa Number
            return Float64(tmp)
        else
            @warn "Invalid temperature choice '$(tmp)'"
            return nothing
        end
    end

    """
    **Parse required argument from request list.**

    Arguments:
    - `request::Array{String,1}`   list of request commands and parameters.
    - `idx::Int64`                index of the required argument to retrieve.

    Returns:
    - `String`                    the required argument at the specified index
    """
    function _verb_arg(request::Array{String,1}, idx::Int64)::String
        if idx > length(request)
            @warn "Request for initial T(p) structure badly formatted"
            @warn "    Got: $(request)"
            @warn "    Check that all verbs have required parameters (e.g: \"iso\", \"1000\")"
            return atmosphere.UNSET_STR # to intentionally trigger an error
        end
        return request[idx]
    end

    """
    **Process T(p) setup requests.**

    Applies an ordered sequence of commands to set up the atmosphere temperature profile.
    Each command may be followed by a parameter value.

    Arguments:
    - `atmos::Atmos_t`      atmosphere struct instance to modify.
    - `request::Array`      ordered list of commands and their parameters.

    Available commands:
    - `"dry"`: set to dry adiabat from surface
    - `"iso"` + value: set isothermal profile at specified temperature
    - `"str"` + value: apply isothermal stratosphere above tropopause
    - `"loglin"` + value: log-linear profile from surface to top temperature
    - `"csv"` + path: load T(p) from CSV file
    - `"ncdf"` + path: load T(p) from NetCDF file
    - `"add"` + value: add temperature offset to entire profile
    - `"sat"` + gas: apply saturation for specified gas
    - `"surfsat"`: ensure surface is not super-saturated
    - `"ana"`: use Guillot (2010) analytic solution

    Returns:
    - `Bool`                success status.
    """
    function request!(atmos::atmosphere.Atmos_t, request::Array)::Bool

        if !(atmos.is_alloc && atmos.is_param)
            @warn "setpt: Atmosphere is not setup or allocated"
            return false
        end

        succ::Bool = true

        num_req::Int64 = length(request)          # Number of requests
        idx_req::Int64 = 1                        # Index of current request
        str_req::String = atmosphere.UNSET_STR  # String of current request
        prt_req::String = "Setting T(p): "
        while idx_req <= num_req
            # get command
            str_req = strip(lowercase(request[idx_req]))
            prt_req *= str_req*", "

            # handle requests
            if str_req == "dry"
                # dry adiabat from surface
                succ &= setpt.dry_adiabat!(atmos)

            elseif str_req == "str"
                # isothermal stratosphere
                idx_req += 1
                succ &= setpt.stratosphere!(atmos, _verb_arg(request, idx_req))

            elseif str_req == "loglin"
                # log-linear profile between T_surf and T_top
                idx_req += 1
                succ &= setpt.loglinear!(atmos, _verb_arg(request, idx_req))

            elseif str_req == "iso"
                # isothermal profile
                idx_req += 1
                succ &= setpt.isothermal!(atmos, _verb_arg(request, idx_req))

            elseif str_req == "csv"
                # set from csv file
                idx_req += 1
                succ &= setpt.fromcsv!(atmos,_verb_arg(request, idx_req))

            elseif str_req == "ncdf"
                # set from NetCDF file
                idx_req += 1
                succ &= setpt.fromncdf!(atmos,_verb_arg(request, idx_req))

            elseif str_req == "add"
                # add X kelvin from the currently stored T(p)
                idx_req += 1
                succ &= setpt.add!(atmos,_verb_arg(request, idx_req))

            elseif str_req == "surfsat"
                # ensure surface is not super-saturated
                chemistry.restore_composition!(atmos)
                succ &= chemistry._sat_surf!(atmos)

            elseif str_req == "sat"
                # condensing a volatile
                idx_req += 1
                succ &= setpt.saturation!(atmos, _verb_arg(request, idx_req))

            elseif str_req == "ana"
                # analytic solution
                succ &= setpt.analytic!(atmos)
            else
                @warn "Invalid initial state '$str_req'"
                return false
            end

            succ &= atmosphere.calc_layer_props!(atmos)

            # iterate
            idx_req += 1
        end
        @info prt_req[1:end-2]
        return succ
    end

    """
    **Set T(p) by log-pressure interpolation.**

    Interpolates temperature profile from given pressure and temperature arrays
    using log-pressure space for smooth behavior across the atmospheric domain.

    Arguments:
    - `atmos::Atmos_t`      atmosphere struct instance to modify.
    - `pl::Array`           pressure array [Pa] at level edges.
    - `tmpl::Array`         temperature array [K] at level edges.
    - `extrap::Bool`        use linear extrapolation beyond array bounds (default: false).

    Returns:
    - `Bool`                success status.
    """
    function fromarrays!(atmos::atmosphere.Atmos_t, pl::Array, tmpl::Array;
                            extrap::Bool=false)::Bool

        if !(atmos.is_alloc && atmos.is_param)
            @warn "setpt: Atmosphere is not setup or allocated"
            return false
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
                @warn "setpt: input array of pressures is not monotonically increasing"
                return false
            end
        end

        # Extrapolate loaded grid to lower pressures (prevent domain error)
        if (atmos.pl[1] < pl[1])
            pushfirst!(pl,   atmos.pl[1]/1.01)
            pushfirst!(tmpl, tmpl[1])
        end

        # Extrapolate loaded grid to higher pressures
        if (atmos.pl[end] > pl[end])
            push!(pl,   atmos.pl[end]*1.01)
            push!(tmpl, tmpl[end])
        end

        # Interpolate from the loaded grid to the required one
        #   This uses log-pressures in order to make the interpolation behave
        #   reasonably across the entire grid.
        itp::Extrapolation = extrapolate(
                                        interpolate((log10.(pl),),tmpl, Gridded(Linear())),
                                        extrap ? Linear() : Flat()
                                        )
        @. atmos.tmpl = itp(log10(atmos.pl))  # Cell edges
        @. atmos.tmp  = itp(log10(atmos.p ))   # Cell centres

        clamp!(atmos.tmpl, atmos.tmp_floor+0.1, atmos.tmp_ceiling-0.1)
        clamp!(atmos.tmp,  atmos.tmp_floor+0.1, atmos.tmp_ceiling-0.1)

        return true
    end

    """
    **Set T(p) from CSV file.**

    Reads pressure [Pa] and temperature [K] columns from a CSV file and interpolates
    to the model grid. Lines starting with '#' and empty lines are ignored.

    Arguments:
    - `atmos::Atmos_t`      atmosphere struct instance to modify.
    - `fpath::String`       path to CSV file.

    Returns:
    - `Bool`                success status.
    """
    function fromcsv!(atmos::atmosphere.Atmos_t, fpath::String)::Bool

        if !(atmos.is_alloc && atmos.is_param)
            @warn "setpt: Atmosphere is not setup or allocated"
            return false
        end

        # Check file exists
        if !isfile(fpath)
            @warn "setpt: file '$fpath' does not exist"
            return false
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
            @warn "setpt: file contains too few levels (contains $nlev_l edge values)"
            return false
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
                @warn "setpt: Negative pressure(s) in csv file"
                return false
            end
            if t_val <= 0.0
                @warn "setpt: Negative temperature(s) in csv file"
                return false
            end

            # Store
            pl[idx] = p_val
            tmpl[idx] = t_val

            # Iterate
            idx += 1
        end

        # Use fromarrays function to do the rest
        return fromarrays!(atmos, pl, tmpl)
    end

    """
    **Load T(p) from NetCDF file.**

    Reads atmosphere data from a NetCDF file and interpolates to the model grid.
    Also updates surface temperature from the file.

    Arguments:
    - `atmos::Atmos_t`      atmosphere struct instance to modify.
    - `fpath::String`       path to NetCDF file.

    Returns:
    - `Bool`                success status.
    """
    function fromncdf!(atmos::atmosphere.Atmos_t, fpath::String)::Bool

        if !(atmos.is_alloc && atmos.is_param)
            @warn "setpt: Atmosphere is not setup or allocated"
            return false
        end

        # Check file exists
        if !isfile(fpath)
            @warn "setpt: The file '$fpath' does not exist"
            return false
        end

        # Open file
        fpath = abspath(fpath)
        @debug "Setting PT from NetCDF file "
        @debug "ALL DEBUG SUPPRESSED"
        with_logger(MinLevelLogger(current_logger(), Logging.Info-200)) do
            ds = Dataset(fpath,"r")

            # Allocate interleaved temperature profile
            nlev_c::Int64 = length(ds["p"][:])
            arr_n::Int64 = nlev_c + nlev_c + 1
            arr_T::Array{Float64, 1} = zeros(Float64, arr_n)
            arr_P::Array{Float64, 1} = zeros(Float64, arr_n)

            # top
            arr_T[1] = ds["tmpl"][1]
            arr_P[1] = ds["pl"][1]

            # middle
            idx::Int64 = 0
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

        return true
    end # end load_ncdf

    """
    **Set isothermal profile.**

    Sets the entire atmosphere to a single constant temperature.

    Arguments:
    - `atmos::Atmos_t`      atmosphere struct instance to modify.
    - `set_tmp`             temperature [K] (numeric, or `"teq"`, `"tsurf"`).

    Returns:
    - `Bool`                success status.
    """
    function isothermal!(atmos::atmosphere.Atmos_t, set_tmp)
        if !(atmos.is_alloc && atmos.is_param)
            @warn "setpt: Atmosphere is not setup or allocated"
            return false
        end


        set_tmp = _parse_tmp_str(atmos, set_tmp)
        isnothing(set_tmp) && (return false)

        fill!(atmos.tmpl, set_tmp)
        fill!(atmos.tmp , set_tmp)

        return true
    end

    """
    **Add temperature offset to entire profile.**

    Adds a constant temperature offset to every level in the atmosphere.

    Arguments:
    - `atmos::Atmos_t`      atmosphere struct instance to modify.
    - `delta`               temperature offset [K] to add.

    Returns:
    - `Bool`                success status.
    """
    function add!(atmos::atmosphere.Atmos_t, delta)::Bool
        if !(atmos.is_alloc && atmos.is_param)
            @warn "setpt: Atmosphere is not setup or allocated"
            return false
        end

        delta = _parse_tmp_str(atmos, delta)
        isnothing(delta) && (return false)

        @. atmos.tmpl += delta
        @. atmos.tmp  += delta

        return true
    end

    """
    **Set T(p) to dry adiabat.**

    Integrates upward from the surface temperature using the dry adiabatic lapse rate.
    Uses real gas properties (cp, density) at each level based on current composition.

    Arguments:
    - `atmos::Atmos_t`      atmosphere struct instance to modify.

    Returns:
    - `Bool`                success status.
    """
    function dry_adiabat!(atmos::atmosphere.Atmos_t)::Bool
        # Validate input
        if !(atmos.is_alloc && atmos.is_param)
            @warn "setpt: Atmosphere is not setup or allocated"
            return false
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


        return true
    end

    """
    **Apply isothermal stratosphere.**

    Clamps all temperatures to be at least `strat_tmp`, effectively creating
    an isothermal stratosphere in regions where T < strat_tmp.

    Arguments:
    - `atmos::Atmos_t`      atmosphere struct instance to modify.
    - `strat_tmp`           stratospheric temperature [K] (numeric, or `"teq"`, `"tsurf"`).

    Returns:
    - `Bool`                success status.
    """
    function stratosphere!(atmos::atmosphere.Atmos_t, strat_tmp)::Bool
        strat_tmp = _parse_tmp_str(atmos, strat_tmp)
        isnothing(strat_tmp) && (return false)

        clamp!(atmos.tmp,  strat_tmp, atmos.tmp_ceiling)
        clamp!(atmos.tmpl, strat_tmp, atmos.tmp_ceiling)
        return true
    end

    """
    **Set log-linear T(p) profile.**

    Creates a temperature profile that varies linearly with log-pressure,
    from the surface temperature down to `top_tmp` at the top of atmosphere.

    Arguments:
    - `atmos::Atmos_t`      atmosphere struct instance to modify.
    - `top_tmp`             temperature at TOA [K] (numeric, or `"teq"`, `"tsurf"`).

    Returns:
    - `Bool`                success status.
    """
    function loglinear!(atmos::atmosphere.Atmos_t, top_tmp)::Bool

        if !(atmos.is_alloc && atmos.is_param)
            @warn "setpt: Atmosphere is not setup or allocated"
            return false
        end

        top_tmp = _parse_tmp_str(atmos, top_tmp)
        isnothing(top_tmp) && (return false)

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

        return true
    end



    """
    **Enforce saturation constraint.**

    Adjusts temperature profile to satisfy T ≥ T_dew for a specified gas.
    Does not modify gas volume mixing ratios, only temperature.

    Arguments:
    - `atmos::Atmos_t`      atmosphere struct instance to modify.
    - `gas::String`         name of gas to check for saturation (e.g., `"H2O"`, `"CO2"`).
    - `dTdew::Float64`      temperature offset below dew point [K] for stability (default: 0.05).

    Returns:
    - `Bool`                success status.
    """
    function saturation!(atmos::atmosphere.Atmos_t, gas::String; dTdew::Float64=0.05)

        if !(atmos.is_alloc && atmos.is_param)
            @warn "setpt: Atmosphere is not setup or allocated"
            return false
        end

        # gas is present?
        if !(gas in atmos.gas_names)
            return true
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
            if tmp <= Tdew + dTdew
                return (Tdew - dTdew, true)
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
        if gas == "H2O"
            atmosphere.set_cloud!(atmos; from_yield=false)
        end

        return true
    end

    """
    **Set T(p) using Guillot (2010) analytic solution.**

    Uses the analytical radiative-equilibrium solution from Guillot (2010)
    for irradiated planetary atmospheres ( 	https://doi.org/10.1051/0004-6361/200913396 ).

    Arguments:
    - `atmos::Atmos_t`      atmosphere struct instance to modify.

    Returns:
    - `Bool`                success status.
    """
    function analytic!(atmos::atmosphere.Atmos_t)::Bool

        if !(atmos.is_alloc && atmos.is_param)
            @warn "setpt: Atmosphere is not setup or allocated"
            return false
        end

        # Evalulate Tirr from the instellation
        Tirr = (atmos.instellation / phys.σSB)^0.25

        # Evalulate Tint from flux
        Tint = (atmos.flux_int / phys.σSB)^0.25

        # Evalulate cell-centre temperatures
        for i in 1:atmos.nlev_c
            # get LW optical depth
            τ = guillot.eval_tau(atmos.p[i], atmos.g[i])

            # set temperature
            atmos.tmp[i] = guillot.eval_T4_cos(τ, Tint, Tirr, atmos.zenith_degrees)^0.25
        end

        # Set cell-edge temperatures
        atmosphere.set_tmpl_from_tmp!(atmos)
        return true
    end

end # end module
