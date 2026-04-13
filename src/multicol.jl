# Contains the multicol module

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

"""
**Main multicol module for handling multiple 1D-atmospheric columns**

The user should create a single instance of `Atmos_t` first, and then use that to
instantiate the `Globe_t` struct.

SOCRATES is not thread-safe, so the original `atmos` will be modified in-place.
Other variables will be copied to/from it as needed, using auxiliary `atmos` variables.

"""
module multicol

    # System libraries
    using Printf
    using Logging

    # Local files
    import ..atmosphere
    import ..phys

    # Globals - TODO: make configurable
    redist_Pmid::Float64 = 1e5      # 1 bar
    redist_Pwid::Float64 = 1.0      # 1 log unit

    # For copying variables between atmospheres
    SHARED_TYPES     = Union{AbstractArray, AbstractVector, Number, AbstractDict, AbstractString}
    PROTECTED_FIELDS = Symbol[:name, :num_rt_eval, :tim_rt_eval]
    shared_fields = Symbol[]
    for (field,type) in zip(fieldnames(atmosphere.Atmos_t), atmosphere.Atmos_t.types)
        if (type <: SHARED_TYPES) && !(field in PROTECTED_FIELDS)
            push!(shared_fields, field)
        end
    end

    # Contains data for the globe
    mutable struct Globe_t

        # Meta parameters
        is_constructed::Bool               # whether the globe is constructed

        # Worker atmosphere
        atmos_wrk::atmosphere.Atmos_t      # worker atmosphere for calculations

        # Arrays of atmosphere columns
        ncol::Int64                                 # number of columns
        atmos_arr::Array{atmosphere.Atmos_t,1}      # auxiliary atmospheric columns
        lons_arr::Array{Float64,1}                  # longitudes for each column
        lats_arr::Array{Float64,1}                  # latitudes for each column

        Globe_t() = new()
    end # end Globe_t


    """
    **De-construct the globe**

    This resets the globe and deallocates memory for all atmosphere columns.

    Arguments:
    - `globe::Globe_t`      the globe to de-construct
    """
    function deconstruct!(globe::Globe_t)
        for atmos in globe.atmos_arr
            atmosphere.deallocate!(atmos)
        end
        globe.atmos_arr = Array{atmosphere.Atmos_t,1}(undef,0)
        globe.is_constructed = false
    end

    """
    **Setup the globe with the given parameters**

    Arguments:
    - `globe::Globe_t`              globe struct to setup
    - `atmos::Atmos_t`              original 1D-atmosphere to use as a basis
    - `lons::Array{Float64,1}`      longitudes for each column (degrees)
    - `lats::Array{Float64,1}`      latitudes for each column (degrees)

    Returns:
    - `Bool`                        globe was successfully constructed
    """
    function construct!(globe::Globe_t, atmos::atmosphere.Atmos_t,
                            lons::Array{Float64,1},
                            lats::Array{Float64,1},
                        )::Bool

        # Fail safe
        globe.is_constructed = false
        if !atmos.is_alloc
            @warn "Cannot construct globe because initial atmosphere is not allocated"
            return false
        end

        # Store reference to worker atmosphere - will be updated in place
        globe.atmos_wrk = atmos

        # Set up the column locations
        globe.ncol     = length(lons)
        globe.lons_arr = lons
        globe.lats_arr = lats

        # Check that the longitudes and latitudes are valid
        if length(globe.lats_arr) != globe.ncol
            @warn "Number of longitudes and latitudes must match number of columns: $(globe.ncol)"
            return false
        end
        for i in 1:globe.ncol
            if !(0.0 <= globe.lons_arr[i] <= 360.0)
                @warn "Invalid longitude $(globe.lons_arr[i]) for column $i"
                return false
            end
            if !( abs(globe.lats_arr[i]) <= 90.0)
                @warn "Invalid latitude $(globe.lats_arr[i]) for column $i"
                return false
            end
        end

        # Instantiate array of auxiliary atmospheric columns
        globe.atmos_arr = atmosphere.Atmos_t[atmosphere.Atmos_t() for _ in 1:globe.ncol]

        # Print shared fields
        # @debug "Globe shared fields: $(shared_fields)"

        # Copy the original atmosphere to each column
        for i in 1:globe.ncol

            # Make copy of atmosphere for the column
            globe.atmos_arr[i] = deepcopy(atmos)

            # Allocate arrays
            for field in shared_fields
                setfield!(globe.atmos_arr[i], field, deepcopy(getfield(atmos, field)))
            end

            # Set name of atmosphere for the column
            globe.atmos_arr[i].name = atmos.name * @sprintf(".col%03d", i)

            # Set the zenith angle for the column based on its latitude and longitude
            globe.atmos_arr[i].zenith_degrees = atmosphere.calc_zenith_angle(globe.lons_arr[i], globe.lats_arr[i])
            globe.atmos_arr[i].s0_fact        = 1.0

            # Set TOA heating
            globe.atmos_arr[i].toa_heating = atmosphere.calc_toa_heating(globe.atmos_arr[i])

            # Set deep heating
            atmosphere.set_deep_heating!(globe.atmos_arr[i],
                                            redist_Pmid,
                                            redist_Pwid,
                                            0.0, # always flux_rel = 0
                                            0.0, # initially zero, but set later
                                            "pressure", "boundary_flux", "abs")

            # Check that arrays have correct length
            if length(globe.atmos_arr[i].tmp) != atmos.nlev_c
                @warn "Globe column $i has unexpected number of vertical levels ($(length(globe.atmos_arr[i].tmp)))"
                return false
            end
        end

        # Set the globe as constructed
        globe.is_constructed = true
        @debug "Globe constructed with $(globe.ncol) columns"

        return true
    end # end construct!

    """
    **Copy arrays and other data from one atmosphere struct to another.**

    Arguments:
    - `dst::Atmos_t`        destination atmosphere to copy to
    - `src::Atmos_t`        source atmosphere to copy from

    Returns:
    - `Bool`                whether the copy was successful
    """
    function copy_atmos_fields!(dst::atmosphere.Atmos_t, src::atmosphere.Atmos_t)::Bool
        for field in shared_fields
            setfield!(dst, field, getfield(src, field))
        end
        return true
    end

    """
    **Wrapper calling a function across multiple columns in a globe**

    Arguments:
    - `globe::Globe_t`      globe struct instance containing the columns
    - `func::Function`      function to call for each column, which must take an `Atmos_t` as its first argument
    - `args...`             additional arguments to pass to the function, after the `Atmos_t` argument
    - `kwargs...`           additional keyword arguments to pass to the function

    Returns:
    - `succ::Bool`         whether the function succeeded for all cases
    """
    function call_for_globe!(globe::multicol.Globe_t,
                                    func::Function, args...; kwargs...)::Bool

        succ::Bool=true
        for i in 1:globe.ncol
            # Copy data from aux to worker
            succ &= copy_atmos_fields!(globe.atmos_wrk, globe.atmos_arr[i])

            # Run the function for this column, using the worker
            succ &= func(globe.atmos_wrk, args...; kwargs...)

            # Copy data from worker back to aux
            succ &= copy_atmos_fields!(globe.atmos_arr[i], globe.atmos_wrk)
        end
        return succ
    end

    """
    **Wrapper for calling a function agnostically of either an atmosphere or globe**

    If `ag` is a 1D-atmosphere type, `func` will be called directly on it. Alternatively,
    if `ag` is a multicol globe type, `func` will be called for each column using `call_for_globe!`.

    Arguments:
    - `ag::Union{Atmos_t,Globe_t}`      either an atmosphere or globe instance
    - `func::Function`                  function to call, which must take an `Atmos_t` as its first argument
    - `args...`                         additional arguments to pass to the function
    - `kwargs...`                       additional keyword arguments to pass to the function

    Returns:
    - `succ::Bool`                      whether the function succeeded for all cases
    """
    function call_agnostic!(ag::Union{atmosphere.Atmos_t,Globe_t}, func::Function, args...; kwargs...)::Bool
        if ag isa atmosphere.Atmos_t
            return func(ag, args...; kwargs...)
        elseif ag isa Globe_t
            return call_for_globe!(ag, func, args...; kwargs...)
        else
            @warn "Invalid type for call_agnostic!: $(typeof(ag))"
            return false
        end
    end

    """
    **Set heat-redistribution flux profiles for globe's columns.**

    Update heat-redistribution fluxes for all columns in globe.

    Arguments:
    - `globe::Globe_t`              globe struct containing the columns to update
    - `fluxes::Array{Float64,1}`    single flux value for each column

    Returns:
    - `succ::Bool                   whether the update was successful
    """
    function set_redist!(globe::Globe_t, fluxes::Array{Float64,1})::Bool

        succ::Bool = true

        for (i,atmos) in enumerate(globe.atmos_arr)
            succ &= atmosphere.set_deep_heating!(
                        atmos,
                        atmos.deepheat_Pmid, # do not change previously-stored value
                        atmos.deepheat_Pwid, # ^
                        0.0,                 # flux_rel = 0
                        fluxes[i],           # flux_abs = redist flux for this column
                        "pressure", "boundary_flux", "abs"
                    )
        end
        return succ
    end

end # end module
