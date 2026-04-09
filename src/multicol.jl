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


    # Contains data for the globe
    mutable struct Globe_t

        # Meta parameters
        is_constructed::Bool               # whether the globe is constructed

        # Worker atmosphere
        atmos_wrk::atmosphere.Atmos_t      # worker atmosphere for calculations

        # Arrays of atmosphere columns
        num_cols::Int64                             # number of columns
        atmos_arr::Array{atmosphere.Atmos_t,1} # auxiliary atmospheric columns
        lons_arr::Array{Float64,1}                  # longitudes for each column
        lats_arr::Array{Float64,1}                  # latitudes for each column

        Globe_t() = new()
    end # end Globe_t


    """
    **Destruct the globe**

    Arguments:
    - `globe::Globe_t`      the globe to destruct
    """
    function destruct!(globe::Globe_t)
        for atmos in globe.atmos_arr
            atmosphere.deallocate!(atmos)
        end
        globe.atmos_arr = Array{atmosphere.Atmos_t,1}(undef,0)
        globe.lons_arr = Array{Float64,1}(undef,0)
        globe.lats_arr  = Array{Float64,1}(undef,0)
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
                            lons::Array{Float64,1}, lats::Array{Float64,1})::Bool

        # Ensure that the globe is not already constructed
        if globe.is_constructed
            @warn "Cannot construct globe because it is already constructed"
            return false
        end

        # Store reference to worker atmosphere
        globe.atmos_wrk = atmos

        # Ensure that worker is allocated
        if !globe.atmos_wrk.is_alloc
            @warn "Cannot construct globe because initial atmosphere is not allocated"
            return false
        end

        # Set up the column locations
        globe.num_cols  = length(lons)
        globe.lons_arr  = lons
        globe.lats_arr  = lats

        # Check that the longitudes and latitudes are valid
        if length(globe.lats_arr) != globe.num_cols
            @warn "Number of longitudes and latitudes must match number of columns: $(globe.num_cols)"
            return false
        end
        for i in 1:globe.num_cols
            if !(0.0 <= globe.lons_arr[i] <= 360.0)
                @warn "Invalid longitude $(globe.lons_arr[i]) for column $i"
                return false
            end
            if !( abs(globe.lats_arr[i]) <= 90.0)
                @warn "Invalid latitude $(globe.lats_arr[i]) for column $i"
                return false
            end
        end

        # Allocate array of auxiliary atmospheric columns
        globe.atmos_arr = atmosphere.Atmos_t[atmosphere.Atmos_t() for i in 1:globe.num_cols]

        # Copy the original atmosphere to each column
        for i in 1:globe.num_cols

            # Make copy of atmosphere for the column
            globe.atmos_arr[i] = deepcopy(atmos)

            # Set name of atmosphere for the column
            globe.atmos_arr[i].name = atmos.name * @sprintf(".col%03d", i)

            # Set the zenith angle for the column based on its latitude and longitude
            globe.atmos_arr[i].zenith_degrees = atmosphere.calc_zenith_angle(globe.lons_arr[i], globe.lats_arr[i])
            globe.atmos_arr[i].s0_fact        = 1.0

            # Set TOA heating
            globe.atmos_arr[i].toa_heating = atmosphere.calc_toa_heating(globe.atmos_arr[i])

            # Check that arrays have correct length
            if length(globe.atmos_arr[i].tmp) != atmos.nlev_c
                @warn "Globe column $i has unexpected number of vertical levels ($(length(globe.atmos_arr[i].tmp)))"
                return false
            end
        end

        # Set the globe as constructed
        globe.is_constructed = true
        @debug "Globe constructed with $(globe.num_cols) columns"

        return true
    end # end construct!

    """
    **Copy arrays and other data from one atmosphere to another.**

    Arguments:
    - `dst::Atmos_t`        destination atmosphere to copy to
    - `src::Atmos_t`        source atmosphere to copy from

    Returns:
    - `Bool`                whether the copy was successful
    """
    function copy_atmos_fields!(dst::atmosphere.Atmos_t, src::atmosphere.Atmos_t)::Bool
        for field in fieldnames(src)
            # Only copy fields which are of these types
            if getfield(src, field) <: Union{Float64, Int64, Int, Bool, Dict, String, Char}
                setfield!(dst, field, getfield(src, field))
            end
        end
        return true
    end

end # end module
