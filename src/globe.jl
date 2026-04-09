# Contains the globe module

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

"""
**Main globe module for handling multiple 1D-atmospheric columns**

The user should create a single instance of `Atmos_t` first, and then use that to
instantiate the `Globe_t` struct.

SOCRATES is not thread-safe, so the original `atmos` will be modified in-place.
Other variables will be copied to/from it as needed, using auxiliary `atmos` variables.

"""
module globe

    # System libraries
    using Printf
    using Logging

    # Local files
    import ..atmosphere
    import ..phys


    # Contains data for the globe
    mutable struct Globe_t

        # Meta parameters
        is_ready::Bool64                        # whether the globe is ready to run

        # Worker atmosphere
        atmos_wrk::atmosphere.Atmosphere_t      # worker atmosphere for calculations

        # Arrays of atmosphere columns
        num_cols::Int64                             # number of columns
        atmos_arr::Array{atmosphere.Atmosphere_t,1} # auxiliary atmospheric columns
        longs_arr::Array{Float64,1}                 # longitudes for each column
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
        globe.atmos_arr = Array{atmosphere.Atmosphere_t,1}(undef,0)
        globe.longs_arr = Array{Float64,1}(undef,0)
        globe.lats_arr  = Array{Float64,1}(undef,0)
        globe.is_ready  = false
    end

    """
    **Setup the globe with the given parameters**

    Arguments:
    - `globe::Globe_t`              the globe to setup
    - `atmos::Atmosphere_t`         the original 1D-atmosphere to use as a basis
    - `locs::Array{Tuple,1}`        the longitudes and latitudes for each column (degrees)
    """
    function construct!(globe::Globe_t, atmos::atmosphere.Atmosphere_t, locs::Array{Tuple,1})

        # Destruct the globe if it is already ready
        if globe.is_ready
            @warn "Cannot construct globe because it is already ready"
            return false
        end

        # Store reference to worker atmosphere
        globe.atmos_wrk = atmos

        # Set up the column locations
        globe.num_cols  = length(locs)
        globe.longs_arr = Float64[loc[1] for loc in locs]
        globe.lats_arr  = Float64[loc[2] for loc in locs]

        # Check that the longitudes and latitudes are valid
        for i in 1:globe.num_cols
            if !(0.0 < globe.longs_arr[i] < 360.0)
                @error "Invalid longitude $(globe.longs_arr[i]) for column $i"
                return false
            end
            if !( abs(globe.lats_arr[i]) <= 90.0)
                @error "Invalid latitude $(globe.lats_arr[i]) for column $i"
                return false
            end
        end

        # Allocate array of auxiliary atmospheric columns
        globe.atmos_arr = atmosphere.Atmosphere_t[atmosphere.Atmosphere_t() for i in 1:globe.num_cols]

        # Copy the original atmosphere to each column
        for i in 1:globe.num_cols

            # Make copy of atmosphere for the column
            globe.atmos_arr[i] = deepcopy(atmos)

            # Set the zenith angle for the column based on its latitude and longitude
            globe.atmos_arr[i].zenith_angle = atmosphere.calc_zenith_angle(globe.longs_arr[i], globe.lats_arr[i])
            globe.atmos_arr[i].s0_fact      = 1.0

            # Set TOA heating
            globe.atmos_arr[i].toa_heating = atmosphere.calc_toa_heating(globe.atmos_arr[i])
        end

        # Set the globe as ready
        globe.is_ready = true

    end # end construct!



end # end module
