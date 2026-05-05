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
    using Statistics: median

    # Local files
    import ..atmosphere
    import ..phys
    import ..energy: skin_depth

    # For copying variables between atmospheres
    SHARED_TYPES     = Union{AbstractArray, AbstractVector, Number, AbstractDict, AbstractString}
    PROTECTED_FIELDS = Symbol[:num_rt_eval, :tim_rt_eval]
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

        # Redistribution parameters
        redist_Pmid::Array{Float64,1}              # pressure midpoints for heat redistribution [Pa]
        redist_Pwid::Array{Float64,1}              # pressure widths for heat redistribution [log pressure]
        redist_flux::Array{Float64,1}              # heat redistribution flux for each column [W m-2]

        # Planetary boundary conditions (global)
        flux_int::Float64                  # flux loss from planet [W m-2]
        tmp_magma::Float64                 # magma ocean temperature [K]
        tmp_surf::Float64                  # surface temperature [K]

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
    - `redist_flux::Array{Float64,1}`  heat redistribution flux for each column [W m-2]
    - `redist_Pmid::Array{Float64,1}`  pressure midpoints for heat redistribution [Pa]
    - `redist_Pwid::Array{Float64,1}`  pressure widths for

    Returns:
    - `Bool`                        globe was successfully constructed
    """
    function construct!(globe::Globe_t, atmos::atmosphere.Atmos_t,
                            lons::Array{Float64,1},
                            lats::Array{Float64,1},
                            redist_flux::Array{Float64,1},
                            redist_Pmid::Array{Float64,1},
                            redist_Pwid::Array{Float64,1}
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

        # Check that the longitudes and latitudes are valid
        if (length(lats) != globe.ncol) || (length(lons) != globe.ncol)
            @warn "Num of longitudes & latitudes must match num of columns: $(globe.ncol)"
            return false
        end
        for i in 1:globe.ncol
            if !(0.0 <= lons[i] <= 360.0)
                @warn "Invalid longitude $(lons[i]) for column $i"
                return false
            end
            if !( abs(lats[i]) <= 90.0)
                @warn "Invalid latitude $(lats[i]) for column $i"
                return false
            end
        end

        # Set up other variables
        globe.tmp_magma   = atmos.tmp_magma # BC: magma ocean temperature
        globe.tmp_surf    = atmos.tmp_surf  # BC: surface temperature, same for all
        globe.flux_int    = atmos.flux_int  # BC: global flux loss

        # TODO: should be calculated dynamically in set_redist! at some point
        globe.redist_Pmid = redist_Pmid
        globe.redist_Pwid = redist_Pwid
        globe.redist_flux = redist_flux

        # Instantiate array of auxiliary atmospheric columns
        globe.atmos_arr = atmosphere.Atmos_t[atmosphere.Atmos_t() for _ in 1:globe.ncol]

        # Copy the original atmosphere to each column
        for i in 1:globe.ncol

            # Make copy of atmosphere for the column
            globe.atmos_arr[i] = deepcopy(atmos)

            # Allocate arrays and copy other variables
            for field in shared_fields
                try
                    setfield!(globe.atmos_arr[i], field, deepcopy(getfield(atmos, field)))
                catch e
                    @warn "Cannot copy $field to column $i: $e"
                    return false
                end
            end

            # Set name of atmosphere for the column
            globe.atmos_arr[i].name = atmos.name * @sprintf(".col%03d", i)

            # Set the column's longitude and latitude
            globe.atmos_arr[i].col_lon = lons[i]
            globe.atmos_arr[i].col_lat = lats[i]

            # Set the zenith angle for the column based on its latitude and longitude
            globe.atmos_arr[i].zenith_degrees = atmosphere.calc_zenith_angle(lons[i], lats[i])
            globe.atmos_arr[i].s0_fact        = 1.0

            # Set TOA heating
            globe.atmos_arr[i].toa_heating = atmosphere.calc_toa_heating(globe.atmos_arr[i])

            # Set deep heating
            atmosphere.set_deep_heating!(globe.atmos_arr[i],
                                            globe.redist_Pmid[i],
                                            globe.redist_Pwid[i],
                                            0.0,
                                            globe.redist_flux[i],
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

    """Set a quantity to be a given value, for all columns in the globe**

    Arguments:
    - `globe::Globe_t`      globe struct containing the columns to update
    - `quantity::Symbol`    the quantity to set (e.g. :tmp_surf, :flux_int)
    - `value`               the value to set the quantity to for all columns

    Returns:
    - `succ::Bool`         whether the update was successful
    """
    function set_for_globe!(globe::Globe_t, quantity::Symbol, value)::Bool

        # Check quantity exists
        if !hasfield(atmosphere.Atmos_t, quantity)
            @warn "Invalid quantity $quantity for atmosphere"
            return false
        end

        # Set for worker and for all columns
        setfield!(globe.atmos_wrk, quantity, value)
        for i in 1:globe.ncol
            setfield!(globe.atmos_arr[i], quantity, value)
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
    **Set surface boundary condition for columns in globe**

    Arguments:
    - `globe::Globe_t`      globe struct containing the columns to update
    - `sol_type::Int64`     solution type for surface boundary condition (1: Tsurf, 2: Tmagma, 3: flux_int)

    """
    function set_surface_bc!(globe::Globe_t, sol_type::Int64)::Bool

        flux_tot_avg::Float64 = globe.flux_int

        # Set the surface boundary condition for each column based on the solution type
        if sol_type == 1
            # Enforce Tsurf
            for i in 1:globe.ncol
                globe.atmos_arr[i].tmp_surf = globe.tmp_surf
            end

        elseif sol_type == 2
            # Enforce Tmagma, same for all columns, since they couple with same interior.
            #    To close the system, we update skin_d for all columns, which
            #    depends on Tmagma and allows them to take different Tsurf. This requires
            #    the total flux from each column to be consistent (~ flux_int).
            for i in 1:globe.ncol
                # Enforce Tmagma
                globe.atmos_arr[i].tmp_magma = globe.tmp_magma

                # Get median value of total flux across columns
                #    Will eventually tend to consistent value,
                #    describing the total flux from the planet.
                flux_tot_avg = median([atmos.flux_tot[end] for atmos in globe.atmos_arr])

                # Update the skin depth for this column
                globe.atmos_arr[i].skin_d = skin_depth(globe.atmos_arr[i], flux_tot_avg)
            end

        elseif sol_type == 3
            # Enforce flux_int
            for i in 1:globe.ncol
                globe.atmos_arr[i].flux_int = globe.flux_int
            end

        # TODO: elseif sol_type == 4
        # TODO:     Enforce a particular OLR

        else
            @warn "Invalid solution type for surface BC: $sol_type"
            return false
        end

        return true
    end

    """
    **Set heat-redistribution flux profiles for globe's columns.**

    Update heat-redistribution fluxes and profiles for all columns in globe.

    Arguments:
    - `globe::Globe_t`              globe struct containing the columns to update

    Returns:
    - `succ::Bool                   whether the update was successful
    """
    function set_redist!(globe::Globe_t)::Bool

        succ::Bool = true

        # TODO: calculate redist_Pmid, redist_Pwid, and redist_flux
        #       based on the globe's column locations and other parameters.
        # This will probably be determined by a scaling law depending on:
        #   - column temperature
        #   - column composition
        #   - planet radius
        #   - planet rotation rate

        # Loop through columns
        for (i,atmos) in enumerate(globe.atmos_arr)

            # Set heating profile
            succ &= atmosphere.set_deep_heating!(
                        atmos,
                        globe.redist_Pmid[i], # TODO: calculate in this function
                        globe.redist_Pwid[i], # TODO: calculate in this function
                        0.0,                  # flux_rel = 0, always
                        globe.redist_flux[i], # TODO: calculate in this function
                        "pressure", "boundary_flux", "abs"
                    )
        end

        return succ
    end

end # end module
