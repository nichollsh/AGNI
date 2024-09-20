# Load atmosphere from NetCDF file

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module load

    import ..atmosphere

    using NCDatasets
    using LoggingExtras
    using Printf

    """
    **Load atmosphere data from a NetCDF file.**

    This must be called after setup! and allocate! have been called on this atmos struct.
    Does not overwrite energy fluxes within the atmospherem, or other wavelength-dependent
    properties. Only composition, temperatures, layer properties, and boundary conditions
    are updated by this function.

    Arguments:
    - `atmos::Atmos_t`      the atmosphere struct instance to be used.
    - `fname::String`       path to the NetCDF file

    Returns:
    - `ok::Bool`            loaded ok?
    """
    function load_ncdf!(atmos::atmosphere.Atmos_t, fname::String)::Bool

        # Create dataset handle
        fname = abspath(fname)
        @debug "Loading atmos from NetCDF: $fname"

        # Absorb output from these calls, because they spam the Debug logger
        @debug "ALL DEBUG SUPPRESSED"
        with_logger(MinLevelLogger(current_logger(), Logging.Info-200)) do

            ds = Dataset(fname,"r")

            # ----------------------
            # Get dimensions
            nlev_c::Int = ds.dim["nlev_c"]
            ngases::Int = ds.dim["ngases"]

            # Check that these are the same as the allocated atmos struct
            if nlev_c != atmos.nlev_c
                @error "Cannot load NetCDF file with mismatching number of levels"
                return false
            end
            if ngases != atmos.gas_num
                @error "Cannot load NetCDF file with mismatching nunmber of gases"
                return false
            end

            # ----------------------
            # Load scalar quantities
            atmos.tmp_surf =        ds["tmp_surf"][]        # Surface brightness temperature [K]
            atmos.flux_int  =       ds["flux_int"][]        # Internal flux [W m-2]
            atmos.instellation =    ds["instellation"][]    # Solar flux
            atmos.s0_fact =         ds["inst_factor"][]     # Scale factor applied to instellation
            atmos.albedo_b =        ds["bond_albedo"][]     # Bond albedo used to scale-down instellation
            atmos.zenith_degrees =  ds["zenith_angle"][]    # Zenith angle of direct stellar radiation
            atmos.toa_heating =     ds["toa_heating"][]     # TOA SW BC - will be overwritten if radtrans is called
            atmos.tmp_magma =       ds["tmagma"][]          # Magma temperature
            atmos.tmp_floor =       ds["tfloor"][]          # Minimum temperature
            atmos.tmp_ceiling =     ds["tceiling"][]        # Maximum temperature
            atmos.rp =              ds["planet_radius"][]   # Value taken for planet radius
            atmos.grav_surf =       ds["surf_gravity"][]    # Surface gravity
            atmos.skin_d =          ds["cond_skin_d"][]     # Conductive skin thickness
            atmos.skin_k =          ds["cond_skin_k"][]     # Conductive skin thermal conductivity

            # ----------------------
            # Load vector quantities
            atmos.p =           ds["p"][:]
            atmos.pl =          ds["pl"][:]
            atmos.tmp =         ds["tmp"][:]
            atmos.tmpl =        ds["tmpl"][:]
            atmos.cloud_arr_l = ds["cloud_mmr"][:]
            atmos.p_toa =       atmos.pl[1]
            atmos.p_boa =       atmos.pl[end]

            #   gas names
            raw_gases::Array{Char,2} = ds["gases"][:,:]
            num_gas::Int = size(raw_gases)[2]
            input_gases::Array{String,1} = []
            for i in 1:num_gas
                push!(input_gases, strip(String(raw_gases[:,i])))
            end

            # gas VMRs
            raw_vmrs::Array{Float64, 2} = ds["x_gas"][:,:]
            for i in 1:num_gas
                g = input_gases[i]
                atmos.gas_vmr[g]      = zeros(Float64, nlev_c)
                atmos.gas_vmr[g][:]  .= raw_vmrs[i, :]
                atmos.gas_ovmr[g][:] .= atmos.gas_vmr[g][:]
            end

            # recalculate remaining layer properties
            atmosphere.calc_layer_props!(atmos)

            close(ds)

        end # suppress output
        @debug "ALL DEBUG RESTORED"

        return true
    end # end read_ncdf

end
