# Contains the atmosphere module

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module atmosphere

    # System libraries
    using Printf
    using DelimitedFiles
    using PCHIPInterpolation
    using LinearAlgebra
    using Statistics
    using Logging

    # Local files
    include(joinpath(ENV["RAD_DIR"],"julia/src/SOCRATES.jl"))
    import ..phys
    import ..spectrum

    # Contains data pertaining to the atmosphere (fluxes, temperature, etc.)
    mutable struct Atmos_t

        # Code version
        AGNI_VERSION::String

        # Track state of atmos struct
        is_param::Bool      # Params have been set
        is_alloc::Bool      # Arrays have been allocated
        is_solved::Bool     # Current atmosphere state represents the result of a solver
        is_converged::Bool  # Current atmosphere state is converged
        is_out_lw::Bool     # Contains data for LW calculation
        is_out_sw::Bool     # Contains data for SW calculation

        # Directories
        ROOT_DIR::String        # path to AGNI root folder (containing agni.jl)
        OUT_DIR::String         # path to output folder
        THERMO_DIR::String      # path to thermo data
        FC_DIR::String          # path to fastchem exec folder

        # SOCRATES objects
        SOCRATES_VERSION::String
        dimen::SOCRATES.StrDim
        control::SOCRATES.StrCtrl
        spectrum::SOCRATES.StrSpecData
        atm::SOCRATES.StrAtm
        cld::SOCRATES.StrCld
        aer::SOCRATES.StrAer
        bound::SOCRATES.StrBound
        radout::SOCRATES.StrOut

        # Radiation parameters
        num_rt_eval::Int                # Total number of RT evaluations
        all_channels::Bool              # Use all bands?
        spectral_file::String           # Path to spectral file
        star_file::String               # Path to star spectrum
        albedo_b::Float64               # Enforced bond albedo
        zenith_degrees::Float64         # Solar zenith angle [deg]
        toa_heating::Float64            # Derived downward shortwave radiation flux at topmost level [W m-2]
        instellation::Float64           # Solar flux at top of atmopshere [W m-2]
        s0_fact::Float64                # Scale factor to instellation (cronin+14)

        surface_material::String        # Surface material file path
        albedo_s::Float64               # Grey surface albedo when surface=blackbody
        surf_r_arr::Array{Float64,1}    # Spectral surface reflectance
        surf_e_arr::Array{Float64,1}    # Spectral surface emissivity

        tmp_surf::Float64               # Surface brightness temperature [K]
        tmp_int::Float64                # Effective temperature of the planet [K]
        grav_surf::Float64              # Surface gravity [m s-2]
        overlap_method::Int             # Absorber overlap method to be used

        # Band edges
        nbands::Int
        bands_min::Array{Float64,1}    # Lower wavelength [m]
        bands_max::Array{Float64,1}    # Upper wavelength [m]

        # Pressure-temperature grid (with i=1 at the top of the model)
        nlev_c::Int         # Cell centre (count)
        nlev_l::Int         # Cell edge (count)

        p_boa::Float64      # Pressure at bottom [Pa]
        p_toa::Float64      # Pressure at top [Pa]
        rp::Float64         # surface radius [m]

        tmp::Array{Float64,1}   # cc temperature [K]
        tmpl::Array{Float64,1}  # ce temperature
        p::Array{Float64,1}     # cc pressure
        pl::Array{Float64,1}    # ce pressure
        z::Array{Float64,1}     # cc height
        zl::Array{Float64,1}    # ce height

        tmp_floor::Float64      # Temperature floor to prevent numerics [K]
        tmp_ceiling::Float64    # Temperature ceiling to prevent numerics [K]

        # Conductive skin
        skin_d::Float64         # skin thickness [m]
        skin_k::Float64         # skin thermal conductivity [W m-1 K-1] (You can find reasonable values here: https://doi.org/10.1016/S1474-7065(03)00069-X)
        tmp_magma::Float64      # Mantle temperature [K]

        # Gas variables (incl gases which are not in spectralfile)
        gas_num::Int                                # Number of gases
        gas_names::Array{String,1}                  # List of gas names
        gas_vmr::Dict{String, Array{Float64,1}}     # Layer volume mixing ratios in dict, (key,value) = (gas_name,array)
        gas_sat::Dict{String, Array{Bool, 1}}       # Layer is saturated or cold-trapped
        gas_dat::Dict{String, phys.Gas_t}           # struct variables containing thermodynamic data
        gas_yield::Dict{String, Array{Float64,1}}   # condensate yield [kg/m^2] at each level (can be negative, representing evaporation)
        gas_ovmr::Dict{String, Array{Float64,1}}    # original VMR values at model initialisation
        condensates::Array{String, 1}               # List of condensing gases (strings)
        condense_any::Bool                          # length(condensates)>0 ?

        # Gases (only those in spectralfile)
        gas_soc_num::Int
        gas_soc_names::Array{String,1}

        # Layers' average properties
        thermo_funct::Bool                  # use temperature-dependent evaluation of thermodynamic properties
        layer_density::Array{Float64,1}     # density [kg m-3]
        layer_mmw::Array{Float64,1}         # mean molecular weight [kg mol-1]
        layer_cp::Array{Float64,1}          # heat capacity at const-p [J K-1 kg-1]
        layer_kc::Array{Float64,1}          # thermal conductivity at const-p [W m-1 K-1]
        layer_grav::Array{Float64,1}        # gravity [m s-2]
        layer_thick::Array{Float64,1}       # geometrical thickness [m]
        layer_mass::Array{Float64,1}        # mass per unit area [kg m-2]

        # Calculated bolometric radiative fluxes (W m-2)
        flux_int::Float64                   # Effective flux  [W m-2] for sol_type=3
        target_olr::Float64                 # Target olr [W m-2] for sol_type=4

        flux_d_lw::Array{Float64,1}         # down component, lw
        flux_u_lw::Array{Float64,1}         # up component, lw
        flux_n_lw::Array{Float64,1}         # net upward, lw

        flux_d_sw::Array{Float64,1}         # down component, sw
        flux_u_sw::Array{Float64,1}         # up component, sw
        flux_n_sw::Array{Float64,1}         # net upward, sw

        flux_d::Array{Float64,1}            # down component, lw+sw
        flux_u::Array{Float64,1}            # up component, lw+sw
        flux_n::Array{Float64,1}            # net upward, lw+sw

        # Calculated per-band radiative fluxes (W m-2)
        band_d_lw::Array{Float64,2}         # down component, lw
        band_u_lw::Array{Float64,2}         # up component, lw
        band_n_lw::Array{Float64,2}         # net upward, lw

        band_d_sw::Array{Float64,2}         # down component, sw
        band_u_sw::Array{Float64,2}         # up component, sw
        band_n_sw::Array{Float64,2}         # net upward, sw

        # Contribution function (to outgoing flux) per-band
        contfunc_norm::Array{Float64,2}     # LW only, and normalised by maximum value

        # Sensible heating
        C_d::Float64                        # Turbulent exchange coefficient [dimensionless]
        U::Float64                          # Wind speed [m s-1]
        flux_sens::Float64                  # Turbulent flux

        # Convection
        mask_c::Array{Bool,1}               # Layers transporting convective flux
        flux_cdry::Array{Float64,1}         # Dry convective fluxes from MLT
        Kzz::Array{Float64,1}               # Eddy diffusion coefficient from MLT

        # Conduction
        flux_cdct::Array{Float64,1}         # Conductive flux [W m-2]

        # Phase change
        phs_tau_mix::Float64                # Time scale (mixed composition)
        phs_tau_sgl::Float64                # Time scale (single gas)
        phs_wrk_df::Array{Float64,1}        # work array: flux difference
        phs_wrk_fl::Array{Float64,1}        # work array: edge fluxes
        flux_l::Array{Float64, 1}           # Latent heat energy flux [W m-2]
        mask_l::Array{Bool,1}               # Layers transporting latent heat

        # Clouds
        cond_alpha::Float64                 # Condensate retention fraction (i.e. how much goes into forming clouds vs rainout)
        cloud_arr_r::Array{Float64,1}       # Characteristic dimensions of condensed species [m].
        cloud_arr_l::Array{Float64,1}       # Mass mixing ratios of condensate. 0 : saturated water vapor does not turn liquid ; 1 : the entire mass of the cell contributes to the cloud
        cloud_arr_f::Array{Float64,1}       # Total cloud area fraction in layers. 0 : clear sky cell ; 1 : the cloud takes over the entire area of the Cell
        cloud_val_r::Float64                # \
        cloud_val_l::Float64                #  |-> Default scalar values to above arrays
        cloud_val_f::Float64                # /

        # Cell-internal heating
        ediv_add::Array{Float64, 1}     # Additional energy dissipation inside each cell [W m-3] (e.g. from advection)

        # Total energy flux
        flux_dif::Array{Float64,1}      # Flux lost at each level [W m-2]
        flux_tot::Array{Float64,1}      # Total upward-directed flux at cell edges [W m-2]

        # Heating rate felt at each level [K/day]
        heating_rate::Array{Float64,1}

        # FastChem stuff
        fastchem_flag::Bool             # Fastchem enabled?
        fastchem_work::String           # Path to fastchem working directory

        # Observing properties
        transspec_p::Float64            # (INPUT PARAMETER) pressure level probed in transmission [Pa]
        transspec_r::Float64            # planet radius probed in transmission [m]
        transspec_m::Float64            # mass [kg] enclosed by transspec_r
        transspec_rho::Float64          # bulk density [kg m-3] implied by r and m
        interior_rho::Float64           # interior density [kg m-3]
        interior_mass::Float64          # interior mass [kg]

        Atmos_t() = new()
    end

    # Deallocate SOCRATES arrays
    function deallocate!(atmos::atmosphere.Atmos_t)
        SOCRATES.deallocate_atm(     atmos.atm)

        SOCRATES.deallocate_cld(     atmos.cld)
        SOCRATES.deallocate_cld_prsc(atmos.cld)

        SOCRATES.deallocate_aer(     atmos.aer)
        SOCRATES.deallocate_aer_prsc(atmos.aer)

        SOCRATES.deallocate_bound(   atmos.bound)
        SOCRATES.deallocate_out(     atmos.radout)

        atmos.is_alloc = false
        return nothing
    end

    """
    **Set parameters of the atmosphere.**

    This is used to set up the struct immediately after creation. It must be
        done before `allocate!()` is called, and before any RT calculcations
        are performed.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `ROOT_DIR::String`                AGNI root directory.
    - `OUT_DIR::String`                 Output directory.
    - `spfile::String`                  path to spectral file.
    - `instellation::Float64`           bolometric solar flux at the top of the atmosphere [W m-2]
    - `s0_fact::Float64`                scale factor to account for planetary rotation (i.e. S_0^*/S_0 in Cronin+14)
    - `albedo_b::Float64`               bond albedo scale factor applied to instellation in order to imitate shortwave reflection
    - `zenith_degrees::Float64`         angle of radiation from the star, relative to the zenith [degrees].
    - `tmp_surf::Float64`               effective surface temperature to provide upward longwave flux at the bottom of the atmosphere [K].
    - `gravity::Float64`                gravitational acceleration at the surface [m s-2].
    - `radius::Float64`                 planet radius at the surface [m].
    - `nlev_centre::Int`                number of model levels.
    - `p_surf::Float64`                 total surface pressure [bar].
    - `p_top::Float64`                  total top of atmosphere pressure [bar].
    - `mf_dict::Dict`                   dictionary of VMRs in the format (key,value)=(gas,mf).
    - `mf_path::String`                 path to file containing VMRs at each level.

    Optional arguments:
    - `condensates`                     list of condensates (gas names).
    - `surface_material::String`        surface material (default is "blackbody", but can point to file instead).
    - `albedo_s::Float64`               grey surface albedo used when `surface_material="blackbody"`.
    - `tmp_floor::Float64`              temperature floor [K].
    - `C_d::Float64`                    turbulent heat exchange coefficient [dimensionless].
    - `U::Float64`                      surface wind speed [m s-1].
    - `tmp_magma::Float64`              mantle temperature [K] for sol_type==2.
    - `skin_d::Float64`                 skin thickness [m].
    - `skin_k::Float64`                 skin thermal conductivity [W m-1 K-1].
    - `overlap_method::Int`             gaseous overlap scheme (2: rand overlap, 4: equiv extinct, 8: ro+resort+rebin).
    - `target_olr::Float64`             target OLR [W m-2] for sol_type==4.
    - `tmp_int::Float64`                planet's effective (or internal) brightness temperature [K] for sol_type==3.
    - `all_channels::Bool`              use all channels available for RT?
    - `flag_rayleigh::Bool`             include rayleigh scattering?
    - `flag_gcontinuum::Bool`           include generalised continuum absorption?
    - `flag_continuum::Bool`            include continuum absorption?
    - `flag_aerosol::Bool`              include aersols?
    - `flag_cloud::Bool`                include clouds?
    - `thermo_functions::Bool`          use temperature-dependent thermodynamic properties
    - `use_all_gases::Bool`             store information on all supported gases, incl those not provided in cfg

    Returns:
        Nothing
    """
    function setup!(atmos::atmosphere.Atmos_t,
                    ROOT_DIR::String, OUT_DIR::String,
                    spfile::String,
                    instellation::Float64, s0_fact::Float64, albedo_b::Float64,
                    zenith_degrees::Float64, tmp_surf::Float64,
                    gravity::Float64, radius::Float64,
                    nlev_centre::Int, p_surf::Float64, p_top::Float64,
                    mf_dict, mf_path::String;

                    condensates =               String[],
                    surface_material::String =  "blackbody",
                    albedo_s::Float64 =         0.0,
                    tmp_floor::Float64 =        2.0,
                    C_d::Float64 =              0.001,
                    U::Float64 =                2.0,
                    tmp_magma::Float64 =        3000.0,
                    skin_d::Float64 =           0.05,
                    skin_k::Float64 =           2.0,
                    overlap_method::Int =       4,
                    target_olr::Float64 =       0.0,
                    tmp_int::Float64 =          0.0,
                    all_channels::Bool  =       true,
                    flag_rayleigh::Bool =       false,
                    flag_gcontinuum::Bool =     false,
                    flag_continuum::Bool =      false,
                    flag_aerosol::Bool =        false,
                    flag_cloud::Bool =          false,
                    thermo_functions::Bool =    true,
                    use_all_gases::Bool =       false,
                    fastchem_work::String =     ""
                    )::Bool

        if !isdir(OUT_DIR) && !isfile(OUT_DIR)
            mkdir(OUT_DIR)
        end

        # Code versions
        atmos.AGNI_VERSION = "0.7.0"
        atmos.SOCRATES_VERSION = readchomp(joinpath(ENV["RAD_DIR"],"version"))
        @debug "AGNI VERSION = $(atmos.AGNI_VERSION)"
        @debug "Using SOCRATES at $(ENV["RAD_DIR"])"
        @debug "SOCRATES VERSION = $(atmos.SOCRATES_VERSION)"

        atmos.num_rt_eval = 0

        atmos.dimen =       SOCRATES.StrDim()
        atmos.control =     SOCRATES.StrCtrl()
        atmos.spectrum =    SOCRATES.StrSpecData()
        atmos.atm =         SOCRATES.StrAtm()
        atmos.cld =         SOCRATES.StrCld()
        atmos.aer =         SOCRATES.StrAer()
        atmos.bound =       SOCRATES.StrBound()
        atmos.radout =      SOCRATES.StrOut()

        # Set the parameters (and make sure that they're reasonable)
        atmos.ROOT_DIR =        abspath(ROOT_DIR)
        atmos.OUT_DIR =         abspath(OUT_DIR)
        atmos.THERMO_DIR  =     joinpath(atmos.ROOT_DIR, "res", "thermodynamics")
        atmos.spectral_file =   abspath(spfile)
        atmos.all_channels =    all_channels
        atmos.overlap_method =  overlap_method

        atmos.thermo_funct  = thermo_functions

        atmos.tmp_floor =       max(0.1,tmp_floor)
        atmos.tmp_ceiling =     2.0e4

        if nlev_centre < 25
            nlev_centre = 25
            @warn "Adjusted number of levels to $nlev_centre"
        end
        atmos.nlev_c         =  nlev_centre
        atmos.nlev_l         =  atmos.nlev_c + 1
        atmos.tmp_surf =        max(tmp_surf, atmos.tmp_floor)
        atmos.tmp_int =         tmp_int
        atmos.grav_surf =       max(1.0e-7, gravity)
        atmos.zenith_degrees =  max(min(zenith_degrees,89.8), 0.2)
        atmos.surface_material= surface_material
        atmos.albedo_s =        max(min(albedo_s, 1.0 ), 0.0)
        atmos.instellation =    max(instellation, 0.0)
        atmos.albedo_b =        max(min(albedo_b,1.0), 0.0)
        atmos.s0_fact =         max(s0_fact,0.0)
        atmos.toa_heating =     atmos.instellation * (1.0 - atmos.albedo_b) *
                                    s0_fact * cosd(atmos.zenith_degrees)

        atmos.flux_int =        phys.sigma * (atmos.tmp_int)^4.0
        atmos.target_olr =      max(1.0e-20,target_olr)

        atmos.C_d =             max(0,C_d)
        atmos.U =               max(0,U)

        atmos.tmp_magma =       max(atmos.tmp_floor, tmp_magma)
        atmos.skin_d =          max(1.0e-6,skin_d)
        atmos.skin_k =          max(1.0e-6,skin_k)

        if p_top > p_surf
            @error "p_top must be less than p_surf"
            return false
        end

        atmos.p_toa =           p_top * 1.0e+5 # Convert bar -> Pa
        atmos.p_boa =           p_surf * 1.0e+5
        atmos.rp =              max(1.0, radius)

        # derived statistics
        atmos.interior_mass =   atmos.grav_surf * atmos.rp^2 / phys.G_grav
        atmos.interior_rho  =   3.0 * atmos.interior_mass / ( 4.0 * pi * atmos.rp^3)
        atmos.transspec_p   =   1e2     # 1 mbar = 100 Pa

        # absorption contributors
        atmos.control.l_gas::Bool =         true
        atmos.control.l_rayleigh::Bool =    flag_rayleigh
        atmos.control.l_continuum::Bool =   flag_continuum
        atmos.control.l_cont_gen::Bool =    flag_gcontinuum
        atmos.control.l_aerosol::Bool =     flag_aerosol
        atmos.control.l_cloud::Bool =       flag_cloud
        atmos.control.l_drop::Bool =        flag_cloud
        atmos.control.l_ice::Bool  =        false

        # Initialise temperature grid
        atmos.tmpl = zeros(Float64, atmos.nlev_l)
        atmos.tmp =  zeros(Float64, atmos.nlev_c)

        # Initialise pressure grid with current p_toa and p_boa
        generate_pgrid!(atmos)
        atmos.z          = zeros(Float64, atmos.nlev_c)
        atmos.zl         = zeros(Float64, atmos.nlev_l)
        atmos.layer_grav = ones(Float64, atmos.nlev_c) * atmos.grav_surf

        # Initialise thermodynamics
        atmos.layer_mmw     = zeros(Float64, atmos.nlev_c)
        atmos.layer_density = zeros(Float64, atmos.nlev_c)
        atmos.layer_cp      = zeros(Float64, atmos.nlev_c)
        atmos.layer_kc      = zeros(Float64, atmos.nlev_c)
        atmos.layer_mass    = zeros(Float64, atmos.nlev_c)
        atmos.layer_thick   = zeros(Float64, atmos.nlev_c)

        # Initialise cloud arrays
        atmos.cloud_arr_r   = zeros(Float64, atmos.nlev_c)
        atmos.cloud_arr_l   = zeros(Float64, atmos.nlev_c)
        atmos.cloud_arr_f   = zeros(Float64, atmos.nlev_c)

        # Phase change timescales [seconds]
        atmos.phs_tau_mix = 1.0e4   # mixed composition case
        atmos.phs_tau_sgl = 1.0e4   # single gas case

        # Hardcoded cloud properties
        atmos.cond_alpha    = 0.0     # 0% of condensate is retained (i.e. complete rainout)
        atmos.cloud_val_r   = 1.0e-5  # 10 micron droplets
        atmos.cloud_val_l   = 0.8     # 80% of the saturated vapor turns into cloud
        atmos.cloud_val_f   = 0.8     # 100% of the cell "area" is cloud

        # Read VMRs
        if !isempty(mf_path) && !isempty(mf_dict)
            @error "VMRs provided twice"
            return false
        end
        if isempty(mf_path) && isempty(mf_dict)
            @error "VMRs not provided"
            return false
        end

        mf_source::Int = 1  # source for mf (0: dict, 1: file)
        if isempty(mf_path)
            mf_source = 0
        end

        # The values will be stored in a dict of arrays
        atmos.gas_names =   Array{String}(undef, 0)           # list of names
        atmos.gas_dat =     Dict{String, phys.Gas_t}()        # dict of data structures
        atmos.gas_vmr  =    Dict{String, Array{Float64,1}}()  # dict of VMR arrays
        atmos.gas_ovmr  =   Dict{String, Array{Float64,1}}()  # ^ backup of initial values
        atmos.gas_sat  =    Dict{String, Array{Bool, 1}}()    # dict for saturation
        atmos.gas_yield =   Dict{String, Array{Float64,1}}()  # dict of condensate yield
        atmos.gas_num   =   0                                 # number of gases
        atmos.condensates   =   Array{String}(undef, 0)       # list of condensates

        # Dict input case
        if mf_source == 0
            @info "Composition set by dict"
            for (key, value) in mf_dict  # store as arrays
                gas_valid = strip(key, [' ','\t','\n'])

                # Check if repeated
                if gas_valid in atmos.gas_names
                    @warn "    skipping duplicate gas $gas_valid"

                # Not repeated...
                else
                    # Store VMR
                    atmos.gas_vmr[gas_valid] = ones(Float64, atmos.nlev_c)*value
                    push!(atmos.gas_names, gas_valid)
                    atmos.gas_num += 1
                end
            end
        end # end read VMR from dict

        # File input case
        if mf_source == 1
            # check file
            if !isfile(mf_path)
                @error "Could not read VMR file '$mf_path'"
                return false
            end
            @info "Composition set by file"

            # get header
            mf_head::String =   readline(abspath(mf_path))

            # remove comment symbol at start
            mf_head =           mf_head[2:end]

            # remove whitespace
            mf_head =            replace(mf_head, " " => "")

            # split by column and drop first three
            heads::Array{String,1} = split(mf_head, ",")[4:end]

            # create arrays
            for h in heads
                gas_valid = strip(h, [' ','\t','\n'])
                # Check if repeated
                if !(gas_valid in atmos.gas_names)
                    # Store zero VMR for now
                    atmos.gas_vmr[gas_valid] = zeros(Float64, atmos.nlev_c)
                    push!(atmos.gas_names, gas_valid)
                    atmos.gas_num += 1
                end
            end

            # get body
            mf_body::Array{Float64,2} = readdlm(abspath(mf_path), ',', Float64;
                                                header=false, skipstart=2)
            mf_body = transpose(mf_body)

            # set composition by interpolating with pressure array
            gidx::Int=0
            for li in 4:lastindex(heads)
                # Gas index
                gidx += 1

                # Arrays from file
                arr_p::Array{Float64,1} = mf_body[1,:]
                arr_x::Array{Float64,1} = mf_body[li,:]

                # Extend loaded profile to lower pressures (prevents domain error)
                if arr_p[1] > atmos.p_toa
                    pushfirst!(arr_p,   atmos.p_toa/1.1)
                    pushfirst!(arr_x,   arr_x[1] )
                end

                # Extend loaded profile to higher pressures (prevents domain error)
                if arr_p[end] < atmos.p_boa
                    push!(arr_p,   atmos.p_boa*1.1)
                    push!(arr_x,   arr_x[end])
                end

                # Set up interpolator using file data
                itp::Interpolator = Interpolator(arr_p, arr_x)

                # Set values in atmos struct
                for i in 1:atmos.nlev_c
                    atmos.gas_vmr[atmos.gas_names[gidx]][i] = itp(atmos.p[i])
                end
            end

        end # end read VMR from file

        # add extra gases if required
        if use_all_gases

            # for each gas in the file
            open(joinpath(atmos.THERMO_DIR, "standard.txt"), "r") do hdl
                for gas in readlines(hdl)
                    # comment line
                    if isempty(gas) || occursin("#",gas)
                        continue
                    end
                    gas = strip(gas)

                    # duplicate
                    if gas in atmos.gas_names
                        continue
                    end

                    # add gas
                    atmos.gas_vmr[gas] = zeros(Float64, atmos.nlev_c)
                    push!(atmos.gas_names, gas)
                    atmos.gas_num += 1
                end
            end
        end

        # backup mixing ratios from current state
        for k in keys(atmos.gas_vmr)
            atmos.gas_ovmr[k] = zeros(Float64, atmos.nlev_c)
            atmos.gas_ovmr[k][:] .= atmos.gas_vmr[k][:]
        end

        # set condensation mask and yield values [kg]
        for g in atmos.gas_names
            atmos.gas_sat[g]   = falses(atmos.nlev_c)
            atmos.gas_yield[g] = zeros(Float64, atmos.nlev_c)
        end

        # Check that we actually stored some values
        if atmos.gas_num == 0
            @error "No gases were stored"
            return false
        end

        # Normalise VMRs at each level
        tot_vmr::Float64 = 0.0
        for i in 1:atmos.nlev_c
            # get total
            tot_vmr = 0.0
            for g in atmos.gas_names
                tot_vmr += atmos.gas_vmr[g][i]
            end
            # normalise to 1
            for g in atmos.gas_names
                atmos.gas_vmr[g][i] /= tot_vmr
            end
        end

        # Load gas thermodynamic data
        for g in atmos.gas_names
            atmos.gas_dat[g] = phys.load_gas(atmos.THERMO_DIR, g, atmos.thermo_funct)
        end

        # store condensates
        for c in condensates
            if !atmos.gas_dat[c].stub && !atmos.gas_dat[c].no_sat
                push!(atmos.condensates, c)
            end
        end

        # Validate condensate names
        atmos.condense_any = false
        if length(condensates) > 0
            for c in condensates
                if !(c in atmos.gas_names)
                    @error "Invalid condensate '$c'"
                    return false
                end
            end
            atmos.condense_any = true
        end

        # Except for single gas case, must have at least one non-condensable gas
        if (length(condensates) == atmos.gas_num) && (atmos.gas_num > 1)
            @error "There must be at least one non-condensable gas"
            return false
        end

        # Set initial temperature profile to a small value which still keeps
        #   all of the gases supercritical. This should be a safe condition to
        #   default to, although the user must specify a profile in the cfg.
        for g in atmos.gas_names
            atmos.tmpl[end] = max(atmos.tmpl[end], atmos.gas_dat[g].T_crit+5.0)
            fill!(atmos.tmpl, atmos.tmpl[end])
            fill!(atmos.tmp, atmos.tmpl[end])
        end

        # Fastchem
        # enabled?
        atmos.fastchem_flag = false
        # path to fc working directory
        if isempty(fastchem_work)
            atmos.fastchem_work = joinpath(atmos.OUT_DIR, "fastchem/")
        else
            atmos.fastchem_work = abspath(fastchem_work)
        end
        if ("FC_DIR" in keys(ENV))

            @debug "FastChem env has been set"

            # check folder
            atmos.FC_DIR = abspath(ENV["FC_DIR"])
            if !isdir(atmos.FC_DIR)
                @error "Could not find fastchem folder at '$(atmos.FC_DIR)'"
                return false
            end

            # check executable
            atmos.fastchem_flag = isfile(joinpath(atmos.FC_DIR,"fastchem"))
            if !atmos.fastchem_flag
                @error "Could not find fastchem executable inside '$(atmos.FC_DIR)' "
                return false
            else
                @info "Found FastChem executable"
            end
        else
            @debug "FastChem env variable not set"
        end

        # Record that the parameters are set
        atmos.is_param = true
        atmos.is_solved = false
        atmos.is_converged = false

        @debug "Setup complete"
        return true
    end # end function setup

    """
    **Calculate observed radius and bulk density.**

    This is done at the layer probed in transmission.

    Arguments:
        - `atmos::Atmos_t`          the atmosphere struct instance to be used.

    Returns:
        Nothing
    """
    function calc_observed_rho!(atmos::atmosphere.Atmos_t)

        # transspec_p::Float64            # (INPUT) level probed in transmission [Pa]
        # transspec_r::Float64            # planet radius probed in transmission [m]
        # transspec_m::Float64            # mass [kg] of atmosphere + interior
        # transspec_rho::Float64          # bulk density [kg m-3] implied by r and m

        # get the observed height
        idx::Int = findmin(abs.(atmos.p .- atmos.transspec_p))[2]
        atmos.transspec_r = atmos.z[idx] + atmos.rp

        # get mass of whole atmosphere, assuming hydrostatic
        atmos.transspec_m = atmos.p_boa * 4 * pi * atmos.rp^2 / atmos.grav_surf

        # add mass of the interior component
        atmos.transspec_m += atmos.interior_mass

        # get density of all enclosed by observed layer
        atmos.transspec_rho = 3.0 * atmos.transspec_m / (4.0 * pi * atmos.transspec_r^3)

        return nothing
    end

    """
    **Calculate properties within each layer of the atmosphere (e.g. density, mmw).**

    Assumes that the atmosphere may be treated as an ideal gas.

    Arguments:
        - `atmos::Atmos_t`          the atmosphere struct instance to be used.
        - `ignore_errors::Bool`     do not generate errors from hydrostatic integrator.

    Returns:
        - `ok::Bool`                function result is ok
    """
    function calc_layer_props!(atmos::atmosphere.Atmos_t; ignore_errors::Bool=false)::Bool
        if !atmos.is_param
            error("atmosphere parameters have not been set")
        end

        ok::Bool = true
        dz_max::Float64 = 1e9

        # Set pressure arrays in SOCRATES
        atmos.atm.p[1, :] .= atmos.p[:]
        atmos.atm.p_level[1, 0:end] .= atmos.pl[:]

        # Set MMW at each level
        fill!(atmos.layer_mmw, 0.0)
        for i in 1:atmos.nlev_c
            for gas in atmos.gas_names
                atmos.layer_mmw[i] += atmos.gas_vmr[gas][i] * atmos.gas_dat[gas].mmw
            end
        end

        # Temporary values
        mmr::Float64 = 0.0
        g1::Float64 = 0.0; p1::Float64 = 0.0; t1::Float64 = 0.0
        g2::Float64 = 0.0; p2::Float64 = 0.0; t2::Float64 = 0.0
        dzc::Float64= 0.0; dzl::Float64 = 0.0
        GMpl::Float64 = atmos.grav_surf * (atmos.rp^2.0)

        # Reset arrays
        fill!(atmos.z         ,  0.0)
        fill!(atmos.zl        ,  0.0)
        fill!(atmos.layer_grav,  0.0)
        fill!(atmos.layer_thick, 0.0)
        fill!(atmos.layer_density,0.0)
        fill!(atmos.layer_cp     ,0.0)
        fill!(atmos.layer_kc     ,0.0)
        fill!(atmos.layer_mass   ,0.0)

        # Integrate from bottom upwards
        for i in range(start=atmos.nlev_c, stop=1, step=-1)

            # Set cp, kc at this level
            atmos.layer_cp[i] = 0.0
            atmos.layer_kc[i] = 0.0
            for gas in atmos.gas_names
                mmr = atmos.gas_vmr[gas][i] * atmos.gas_dat[gas].mmw/atmos.layer_mmw[i]
                atmos.layer_cp[i] += mmr * phys.get_Cp(atmos.gas_dat[gas], atmos.tmp[i])
                atmos.layer_kc[i] += mmr * phys.get_Kc(atmos.gas_dat[gas], atmos.tmp[i])
            end

            # Temporarily copy this cp, kc to the level above
            #     since they are needed for doing the hydrostatic integration
            if i > 1
                atmos.layer_cp[i-1] = atmos.layer_cp[i]
                atmos.layer_kc[i-1] = atmos.layer_kc[i]
            end

            # Technically, g and z should be integrated as coupled equations,
            # but here they are not. This loose integration has been found to
            # be reasonable in all of my tests.

            # Integrate from lower edge to centre
            g1 = GMpl / ((atmos.rp + atmos.zl[i+1])^2.0)
            p1 = 0.5 * (atmos.p[i] + atmos.pl[i+1])
            t1 = 0.5 * (atmos.tmp[i] + atmos.tmpl[i+1])
            dzc = phys.R_gas * t1 /
                        (atmos.layer_mmw[i] * g1 * p1) * (atmos.pl[i+1] - atmos.p[i])
            if !ignore_errors
                if (dzc < 1e-20)
                    @error "Height integration resulted in dz <= 0 at level $i (l -> c)"
                    ok = false
                end
                if  (dzc > dz_max)
                    @error "Height integration blew up at level $i (l -> c)"
                    ok = false
                end
            end
            atmos.z[i] = atmos.zl[i+1] + min(dzc,dz_max)

            # Integrate from centre to upper edge
            g2 = GMpl / ((atmos.rp + atmos.z[i])^2.0)
            p2 = 0.5 * (atmos.p[i] + atmos.pl[i])
            t2 = 0.5 * (atmos.tmp[i] + atmos.tmpl[i])
            dzl = phys.R_gas * t2 / (
                        atmos.layer_mmw[i] * g2 * p2) * (atmos.p[i]- atmos.pl[i])
            if !ignore_errors
                if (dzl < 1e-20)
                    @error "Height integration resulted in dz <= 0 at level $i (c -> l)"
                    ok = false
                end
                if (dzl > 1e9)
                    @error "Height integration blew up at level $i (c -> l)"
                    ok = false
                end
            end
            atmos.zl[i] = atmos.z[i] + min(dzl,dz_max)

            # Layer average gravity [m s-2]
            atmos.layer_grav[i] = GMpl / ((atmos.rp + atmos.z[i])^2.0)

            # Layer geometrical thickness [m]
            atmos.layer_thick[i] = atmos.zl[i] - atmos.zl[i+1]
        end

        # Mass (per unit area, kg m-2) and density (kg m-3)
        # Loop from top downards
        for i in 1:atmos.nlev_c
            atmos.layer_mass[i] = (atmos.pl[i+1] - atmos.pl[i])/atmos.layer_grav[i]
            atmos.atm.mass[1, i] = atmos.layer_mass[i]          # pass to SOCRATES

            atmos.layer_density[i] = ( atmos.p[i] * atmos.layer_mmw[i] ) /
                                                        (phys.R_gas * atmos.tmp[i])
            atmos.atm.density[1,i] = atmos.layer_density[i]     # pass to SOCRATES
        end

        return ok
    end

    """
    **Generate pressure grid.**

    Almost-equally log-spaced between p_boa and p_boa. The near-surface layers
    are smaller than they would be on an equally log-spaced grid, to avoid f
    numerical weirdness at the bottom boundary.

    Arguments:
    - `atmos::Atmos_t`              the atmosphere struct instance to be used.
    """
    function generate_pgrid!(atmos::atmosphere.Atmos_t)

        # Allocate arrays
        atmos.p  = zeros(Float64, atmos.nlev_c)
        atmos.pl = zeros(Float64, atmos.nlev_l)

        # First, assign log10'd values...

        # Top and bottom boundaries
        atmos.pl[end] = log10(atmos.p_boa)
        atmos.pl[1]   = log10(atmos.p_toa)

        # Almost-surface layer is VERY small
        atmos.pl[end-1]  = atmos.pl[end]*(1.0-1e-4)

        # Logarithmically-spaced levels above
        atmos.pl[1:end-1] .= collect(Float64, range( start=atmos.pl[1],
                                                      stop=atmos.pl[end-1],
                                                      length=atmos.nlev_l-1))

        # Shrink near-surface layers by stretching all layers above
        p_fact::Float64 = 0.6
        p_mid::Float64 = atmos.pl[end-1]*p_fact + atmos.pl[end-2]*(1.0-p_fact)
        atmos.pl[1:end-2] .= collect(Float64, range( start=atmos.pl[1],
                                                     stop=p_mid,
                                                     length=atmos.nlev_l-2))

        # Set cell-centres at midpoint of cell-edges
        atmos.p[1:end] .= 0.5 .* (atmos.pl[1:end-1] .+ atmos.pl[2:end])

        # Finally, convert arrays to 'real' pressure units
        atmos.p[:]  .= 10.0 .^ atmos.p[:]
        atmos.pl[:] .= 10.0 .^ atmos.pl[:]

        return nothing
    end


    """
    **Allocate atmosphere arrays, prepare spectral files, and final steps.**

    Will not modify spectral file if `stellar_spectrum`` is an empty string.

    Arguments:
    - `atmos::Atmos_t`                 the atmosphere struct instance to be used.
    - `stellar_spectrum::String`       path to stellar spectrum csv file
    """
    function allocate!(atmos::atmosphere.Atmos_t, stellar_spectrum::String)::Bool

        @debug "Allocate atmosphere"
        if !atmos.is_param
            @error "Atmosphere parameters have not been set"
            return false
        end

        atmos.atm.n_profile = 1

        #########################################
        # spectral data
        #########################################

        # Validate files
        if !isfile(atmos.spectral_file)
            @error "Spectral file '$(atmos.spectral_file)' does not exist"
            return false
        end

        spectral_file_run::String  = joinpath([atmos.OUT_DIR, "runtime.sf"])
        spectral_file_runk::String = joinpath([atmos.OUT_DIR, "runtime.sf_k"])

        # Setup spectral file
        socstar::String = ""
        if !isempty(stellar_spectrum)
            @debug "Inserting stellar spectrum"

            if !isfile(stellar_spectrum)
                @error "Stellar spectrum file '$(stellar_spectrum)' does not exist"
                return false
            end

            # File to be loaded
            rm(spectral_file_run , force=true)
            rm(spectral_file_runk, force=true)

            atmos.star_file = abspath(stellar_spectrum)

            # Write spectrum in required format
            socstar = joinpath([atmos.OUT_DIR, "socstar.dat"])
            wl::Array{Float64,1}, fl::Array{Float64,1} = spectrum.load_from_file(atmos.star_file)
            spectrum.write_to_socrates_format(wl, fl, socstar) || return false

            # Insert stellar spectrum and rayleigh scattering, if required
            spectrum.insert_stellar_and_rscatter(atmos.spectral_file,
                                                    socstar, spectral_file_run,
                                                    atmos.control.l_rayleigh)

        else
            # Stellar spectrum was not provided, which is taken to mean that
            #       the spectral file includes it already.
            atmos.star_file = "IN_SPECTRAL_FILE"
            spectral_file_run  = atmos.spectral_file
            spectral_file_runk = atmos.spectral_file*"_k"
            @info "Using pre-existing spectral file without modifications"
        end

        # Validate files


        # Read-in spectral file to be used at runtime
        atmos.control.spectral_file = spectral_file_run

        SOCRATES.set_spectrum(spectrum=atmos.spectrum,
                                spectral_file=atmos.control.spectral_file,
                                l_all_gasses=true)

        # Remove temporary star file if it exists
        rm(socstar, force=true)

        #########################################
        # diagnostics
        #########################################
        atmos.control.l_actinic_flux = Bool(atmos.spectrum.Basic.l_present[2])
        atmos.control.l_photolysis_rate = atmos.spectrum.Photol.n_pathway > 0
        atmos.control.l_flux_div = atmos.spectrum.Photol.n_pathway > 0

        #########################################
        # parameters
        #########################################

        if atmos.all_channels
            n_channel = atmos.spectrum.Basic.n_band
        else
            n_channel = 1
        end

        atmos.nbands = atmos.spectrum.Basic.n_band
        atmos.bands_max = zeros(Float64, atmos.spectrum.Basic.n_band)
        atmos.bands_min = zeros(Float64, atmos.spectrum.Basic.n_band)

        for i in 1:atmos.nbands
            atmos.bands_min[i] = atmos.spectrum.Basic.wavelength_short[i]
            atmos.bands_max[i] = atmos.spectrum.Basic.wavelength_long[i]
        end

        # modules_gen/dimensions_field_cdf_ucf.f90
        npd_direction = 1                   # Maximum number of directions for radiances
        npd_layer = atmos.nlev_c            # Number of layers
        npd_column = 24                     # Maximum number of cloudy subcolumns
        npd_profile = 1

        # BRDF reflections
        npd_max_order = 101                 # Maximum order of spherical harmonics used
        npd_brdf_basis_fnc = 2              # Number of BRDF basis functions
        npd_brdf_trunc = 5                  # Order of BRDF truncation

        # prescribed aerosol optical properties
        npd_profile_aerosol_prsc = 9        # Size allocated for profiles
        npd_opt_level_aerosol_prsc = 170    # Size allocated for levels

        # prescribed cloudy optical properties
        npd_profile_cloud_prsc = 9          # Size allocated for profiles
        npd_opt_level_cloud_prsc = 170      # Size allocated for levels

        # modules_gen/dimensioms_fixed_pcf.f90
        npd_cloud_component        =  1     # Number of components of clouds.
        npd_cloud_type             =  1     # Number of permitted types of clouds.
        npd_overlap_coeff          = 18     # Number of overlap coefficients for cloud
        npd_source_coeff           =  2     # Number of coefficients for two-stream sources
        npd_region                 =  1     # Number of regions in a layer

        atmos.dimen.nd_profile                = npd_profile
        atmos.dimen.nd_flux_profile           = npd_profile
        atmos.dimen.nd_2sg_profile            = npd_profile
        atmos.dimen.nd_radiance_profile       = npd_profile
        atmos.dimen.nd_j_profile              = 1
        atmos.dimen.nd_layer                  = npd_layer
        atmos.dimen.nd_layer_clr              = npd_layer
        atmos.dimen.id_cloud_top              = 1
        atmos.dimen.nd_channel                = n_channel
        atmos.dimen.nd_column                 = npd_column
        atmos.dimen.nd_max_order              = npd_max_order
        atmos.dimen.nd_direction              = npd_direction
        atmos.dimen.nd_viewing_level          = npd_layer
        atmos.dimen.nd_brdf_basis_fnc         = npd_brdf_basis_fnc
        atmos.dimen.nd_brdf_trunc             = npd_brdf_trunc
        atmos.dimen.nd_profile_aerosol_prsc   = npd_profile_aerosol_prsc
        atmos.dimen.nd_profile_cloud_prsc     = npd_profile_cloud_prsc
        atmos.dimen.nd_opt_level_aerosol_prsc = npd_opt_level_aerosol_prsc
        atmos.dimen.nd_opt_level_cloud_prsc   = npd_opt_level_cloud_prsc
        atmos.dimen.nd_cloud_component        = npd_cloud_component
        atmos.dimen.nd_cloud_type             = npd_cloud_type
        atmos.dimen.nd_overlap_coeff          = npd_overlap_coeff
        atmos.dimen.nd_source_coeff           = npd_source_coeff
        atmos.dimen.nd_region                 = npd_region
        atmos.dimen.nd_point_tile             = 1
        atmos.dimen.nd_tile                   = 1
        atmos.dimen.nd_subcol_gen             = 1
        atmos.dimen.nd_subcol_req             = 1
        atmos.dimen.nd_aerosol_mode           = 1

        # Set to true to enable custom surface emission through the
        #   variables `planck%flux_ground(l)` and `d_planck_flux_surface`.
        atmos.control.l_flux_ground = false

        # Allocate arrays, etc.
        SOCRATES.allocate_atm(  atmos.atm,   atmos.dimen, atmos.spectrum)
        SOCRATES.allocate_cld(  atmos.cld,   atmos.dimen, atmos.spectrum)
        SOCRATES.allocate_aer(  atmos.aer,   atmos.dimen, atmos.spectrum)
        SOCRATES.allocate_bound(atmos.bound, atmos.dimen, atmos.spectrum)

        # Fill with zeros - will be set inside of radtrans function at call time
        atmos.bound.flux_ground[1, :] .= 0.0


        ###########################################
        # Number of profiles, and profile coordinates
        ###########################################
        atmos.atm.n_layer =     npd_layer
        atmos.atm.n_profile =   1
        atmos.atm.lat[1] =      0.0
        atmos.atm.lon[1] =      0.0

        ###########################################
        # Range of bands
        ###########################################

        atmos.control.last_band = atmos.spectrum.Basic.n_band
        atmos.control.first_band = 1
        n_band_active = atmos.control.last_band - atmos.control.first_band + 1

        # Map spectral bands into output channels
        if ( (n_channel*floor(n_band_active/n_channel) != n_band_active)  &&
            (atmos.spectrum.Var.n_sub_band >= n_channel) )
            # Number of bands not a multiple of channels so use sub-bands
            atmos.control.l_map_sub_bands = true
        end

        SOCRATES.allocate_control(atmos.control, atmos.spectrum)

        if n_channel == 1
            atmos.control.map_channel[1:atmos.spectrum.Basic.n_band] .= 1
        elseif n_channel == atmos.spectrum.Basic.n_band
            atmos.control.map_channel[1:atmos.spectrum.Basic.n_band] .= 1:n_channel
        else
            @error "n_channel $n_channel != 1 and != $n_band_active not supported "
            return false
        end

        # Calculate the weighting for the bands.
        atmos.control.weight_band .= 1.0

        # 'Entre treatment of optical depth for direct solar flux (0/1/2)'
        # '0: no scaling; 1: delta-scaling; 2: circumsolar scaling'
        atmos.control.i_direct_tau = 1

        atmos.control.n_order_forward = 2


        ############################################
        # Check Options
        ############################################

        if atmos.control.l_rayleigh
            if !Bool(atmos.spectrum.Basic.l_present[3])
                @error "The spectral file contains no rayleigh scattering data"
                return false
            end
        end

        if atmos.control.l_aerosol
            if !Bool(atmos.spectrum.Basic.l_present[11])
                @error "The spectral file contains no aerosol data"
                return false
            end
        end

        if atmos.control.l_gas
            if !Bool(atmos.spectrum.Basic.l_present[5])
                @error "The spectral file contains no gaseous absorption data"
                return false
            end
        end

        if atmos.control.l_continuum
            if !Bool(atmos.spectrum.Basic.l_present[9])
                @error "The spectral file contains no continuum absorption data"
                return false
            end
        end

        if atmos.control.l_cont_gen
            if !Bool(atmos.spectrum.Basic.l_present[19])
                @error "The spectral file contains no generalised continuum absorption data"
                return false
            end
        end


        ################################
        # Gaseous absorption
        #################################

        if atmos.overlap_method == 2
            # random overlap
            atmos.control.i_gas_overlap = SOCRATES.rad_pcf.ip_overlap_random

        elseif atmos.overlap_method == 4
            # equivalent extinction with correct scaling
            atmos.control.i_gas_overlap = SOCRATES.rad_pcf.ip_overlap_k_eqv_scl

        elseif atmos.overlap_method == 8
            # random overlap with resorting and rebinning
            atmos.control.i_gas_overlap = SOCRATES.rad_pcf.ip_overlap_random_resort_rebin

        else
            @error "Invalid overlap method"
            return false
        end

        for j in atmos.control.first_band:atmos.control.last_band
            atmos.control.i_gas_overlap_band[j] = atmos.control.i_gas_overlap
        end

        # Check supported gases
        atmos.gas_soc_num       = atmos.spectrum.Gas.n_absorb               # number of gases
        atmos.gas_soc_names     = Array{String}(undef, atmos.gas_soc_num)   # list of names
        for i_gas in 1:atmos.gas_soc_num   # for each supported gas
            atmos.gas_soc_names[i_gas] =
                SOCRATES.input_head_pcf.header_gas[atmos.spectrum.Gas.type_absorb[i_gas]]
        end

        # VMRs are provided to SOCRATES when radtrans is called
        # For now, they are just stored inside the atmos struct

        # Print info on the gases
        @info "Allocated atmosphere with composition:"
        gas_flags::String = ""
        g::String = ""
        for i in 1:atmos.gas_num
            g = atmos.gas_names[i]
            gas_flags = ""
            if !(g in atmos.gas_soc_names) # flag as not included in radtrans
                gas_flags *= "NO_OPACITY "
            end
            if g in atmos.condensates       # flag as condensable
                gas_flags *= "COND"
            end
            if atmos.gas_dat[g].stub        # flag as containing stub thermo data
                gas_flags *= "NO_THERMO "
            end
            if !isempty(gas_flags)
                gas_flags = "($(gas_flags[1:end-1]))"
            end
            @info @sprintf("    %3d %-7s %6.2e %s", i, g, atmos.gas_vmr[g][end], gas_flags)
        end


        # Calc layer properties using initial temperature profile.
        #    Can generate weird issues since the TOA temperature may be large
        #    large but pressure small, which gives it a low density. With the
        #    hydrostatic integrator, this can cause dz to blow up, especially
        #    with a low MMW gas. Should be okay as long as the T(p) provided
        #    by the user is more reasonable. Silence errors *in this case*.
        calc_layer_props!(atmos, ignore_errors=true)

        ################################
        # Aerosol processes
        #################################

        SOCRATES.allocate_aer_prsc(atmos.aer, atmos.dimen, atmos.spectrum)
        if atmos.control.l_aerosol
            @error "Aerosols not implemented"
            return false
        else
            atmos.dimen.nd_profile_aerosol_prsc   = 1
            atmos.dimen.nd_opt_level_aerosol_prsc = 1
            atmos.dimen.nd_phf_term_aerosol_prsc  = 1

            atmos.aer.mr_source .= SOCRATES.rad_pcf.ip_aersrc_classic_roff

            for i = 1:atmos.spectrum.Dim.nd_aerosol_species
                atmos.aer.mr_type_index[i] = i
            end
        end

        #######################################
        # Clouds
        # see src/aux/input_cloud_cdf.f
        #######################################

        atmos.dimen.nd_profile_cloud_prsc   = 1
        atmos.dimen.nd_opt_level_cloud_prsc = 1
        atmos.dimen.nd_phf_term_cloud_prsc  = 1

        if atmos.control.l_cloud
            # Ice and water mixed homogeneously (-K 1) = ip_cloud_homogen
            # Cloud mixing liquid and ice (-K 2) = ip_cloud_ice_water

            atmos.control.i_cloud_representation = SOCRATES.rad_pcf.ip_cloud_homogen
            atmos.control.i_cloud     = SOCRATES.rad_pcf.ip_cloud_mix_max      # Goes with ip_max_rand
            atmos.control.i_overlap   = SOCRATES.rad_pcf.ip_max_rand           # Maximum/random overlap in a mixed column (-C 2)
            atmos.control.i_inhom     = SOCRATES.rad_pcf.ip_homogeneous        # Homogeneous cloud
            atmos.control.i_st_water  = 5                                      # Liquid Water Droplet type 5 (-d 5)
            atmos.control.i_cnv_water = 5                                      # Convective Liquid Water Droplet type 5
            atmos.control.i_st_ice    = 11                                     # Water Ice type 11 (-i 11)
            atmos.control.i_cnv_ice   = 11                                     # Convective Water Ice type 11
        else
            atmos.control.i_cloud = SOCRATES.rad_pcf.ip_cloud_off # 5 (clear sky)
        end

        SOCRATES.allocate_cld_prsc(atmos.cld, atmos.dimen, atmos.spectrum)

        if atmos.control.l_cloud
            atmos.cld.n_condensed       = 1
            atmos.cld.type_condensed[1] = SOCRATES.rad_pcf.ip_clcmp_st_water
            atmos.cld.n_cloud_type      = 1
            atmos.cld.i_cloud_type[1]   = SOCRATES.rad_pcf.ip_cloud_type_homogen
            atmos.cld.i_condensed_param[1] = SOCRATES.rad_pcf.ip_drop_pade_2
        else
            atmos.cld.n_condensed  = 0
            atmos.cld.n_cloud_type = 0
        end

        atmos.control.i_angular_integration = SOCRATES.rad_pcf.ip_two_stream

        ###########################################
        # Surface properties
        ###########################################
        # allocate reflectance and emissivity arrays
        atmos.surf_r_arr = zeros(Float64, atmos.nbands) # directional-hemispheric reflect.
        atmos.surf_e_arr = ones(Float64, atmos.nbands)  # emissivity

        if atmos.surface_material == "blackbody"
            # grey albedo
            fill!(atmos.surf_r_arr, atmos.albedo_s)
            # Kirchoff's law: set emissivity equal to 1-albedo (spectrally)
            fill!(atmos.surf_e_arr, 1.0-atmos.albedo_s)

        else
            # spectral albedo and emissivity
            # Hapke2012: https://doi.org/10.1017/CBO9781139025683
            # Hammond2024: https://arxiv.org/abs/2409.04386

            # try to find a matching file
            atmos.surface_material = abspath(atmos.surface_material)
            if !isfile(atmos.surface_material)
                @error "Could not find surface albedo file '$(atmos.surface_material)'"
                return false
            end

            # read single-scattering albedo data from file
            _alb_data::Array{Float64,2} = readdlm(atmos.surface_material, Float64)
            _alb_w::Array{Float64, 1} = _alb_data[:,1]     # wavelength [nm]
            _alb_s::Array{Float64, 1} = _alb_data[:,2]     # ss albedo [dimensionless]

            # extrapolate to 0 wavelength, with constant value
            pushfirst!(_alb_w, 0.0)
            pushfirst!(_alb_s, _alb_s[1])

            # extrapolate to large wavelength, with constant value
            push!(_alb_w, 1e10)
            push!(_alb_s, _alb_s[end])

            # convert ss albedo to gamma values (eq 14.3 from Hapke 2012)
            _alb_s[:] .= sqrt.(1 .- _alb_s[:]) # operating in place

            # create interpolator on gamma
            _gamma::Interpolator = Interpolator(_alb_w, _alb_s)

            # use interpolator to fill band values
            ga::Float64 = 0.0
            mu::Float64 = cosd(atmos.zenith_degrees)
            for i in 1:atmos.nbands
                # evaluate gamma at band centre, converting from m to nm
                ga = _gamma(0.5 * 1.0e9 * (atmos.bands_min[i]+atmos.bands_max[i]))

                # calculate dh reflectance (eq 3 from Hammond24)
                atmos.surf_r_arr[i] = (1-ga) / (1+2*ga*mu)

                # calculate emissivity (eq 4 from Hammond24 and eq 15.29 from Hapke 2012)
                atmos.surf_e_arr[i] = 1.0 - ((1-ga)/(1+ga)) * (1- ga/(3*(1+ga)))
            end
        end

        #######################################
        # Output arrays
        #######################################
        atmos.flux_d_lw =         zeros(Float64, atmos.nlev_l)
        atmos.flux_u_lw =         zeros(Float64, atmos.nlev_l)
        atmos.flux_n_lw =         zeros(Float64, atmos.nlev_l)

        atmos.flux_d_sw =         zeros(Float64, atmos.nlev_l)
        atmos.flux_u_sw =         zeros(Float64, atmos.nlev_l)
        atmos.flux_n_sw =         zeros(Float64, atmos.nlev_l)

        atmos.flux_d =            zeros(Float64, atmos.nlev_l)
        atmos.flux_u =            zeros(Float64, atmos.nlev_l)
        atmos.flux_n =            zeros(Float64, atmos.nlev_l)

        atmos.band_d_lw =         zeros(Float64, (atmos.nlev_l,atmos.nbands))
        atmos.band_u_lw =         zeros(Float64, (atmos.nlev_l,atmos.nbands))
        atmos.band_n_lw =         zeros(Float64, (atmos.nlev_l,atmos.nbands))

        atmos.band_d_sw =         zeros(Float64, (atmos.nlev_l,atmos.nbands))
        atmos.band_u_sw =         zeros(Float64, (atmos.nlev_l,atmos.nbands))
        atmos.band_n_sw =         zeros(Float64, (atmos.nlev_l,atmos.nbands))

        atmos.contfunc_norm =     zeros(Float64, (atmos.nlev_c,atmos.nbands))

        atmos.flux_sens =         0.0

        atmos.mask_l =            falses(atmos.nlev_l)      # Phase change
        atmos.mask_c =            falses(atmos.nlev_l)      # Dry convection

        atmos.phs_wrk_df =        zeros(Float64, atmos.nlev_c)  # flux difference
        atmos.phs_wrk_fl =        zeros(Float64, atmos.nlev_l)  # edge fluxes
        atmos.flux_l =            zeros(Float64, atmos.nlev_l)  # Latent heat / phase change
        atmos.flux_cdry =         zeros(Float64, atmos.nlev_l)  # Dry convection
        atmos.flux_cdct =         zeros(Float64, atmos.nlev_l)  # Conduction
        atmos.Kzz =               zeros(Float64, atmos.nlev_l)  # eddy diffusion coeff.

        atmos.flux_tot =          zeros(Float64, atmos.nlev_l)
        atmos.flux_dif =          zeros(Float64, atmos.nlev_c)
        atmos.ediv_add =          zeros(Float64, atmos.nlev_c)
        atmos.heating_rate =      zeros(Float64, atmos.nlev_c)

        # Mark as allocated
        atmos.is_alloc = true
        @debug "Allocate complete"

        return true
    end  # end of allocate


    """
    **Calculate composition assuming chemical equilibrium at each level.**

    Uses FastChem to calculate the gas composition at each level of the atmosphere.
    Volatiles are converted to bulk elemental abundances, which are then provided to
    FastChem alongside the temperature/pressure profile. FastChem is currently called as an
    executable, which is not optimal.

    This function DOES NOT automatically recalculate layer properties (e.g. mmw, density).

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `chem_type::Int`                  chemistry type (see wiki)
    - `write_cfg::Bool`                 write config and elements
    - `tmp_floor::Float64`              temperature floor for T(p) provided to FastChem

    Returns:
    - `state::Int`                      fastchem state (0: success, 1: critical_fail, 2: elem_fail, 3: conv_fail, 4: both_fail)
    """
    function chemistry_eqm!(atmos::atmosphere.Atmos_t, chem_type::Int, write_cfg::Bool; tmp_floor::Float64=200.0)::Int

        @debug "Running equilibrium chemistry"

        # Return code
        state::Int = 0

        # Check fastchem enabled
        if !atmos.fastchem_flag
            @warn "Fastchem is not enabled but `chemistry_eqm!` was called. Have you set FC_DIR?"
            return 1
        end

        # Check minimum temperature
        if maximum(atmos.tmpl) < tmp_floor
            @warn "Temperature profile is entirely too cold for FastChem. Not doing chemistry."
            return 1
        end

        count_elem_nonzero::Int = 0

        # Paths
        execpath::String = joinpath(atmos.FC_DIR,       "fastchem")             # Executable file
        confpath::String = joinpath(atmos.fastchem_work,"config.input")         # Configuration by AGNI
        elempath::String = joinpath(atmos.fastchem_work,"elements.dat")         # Elements by AGNI
        chempath::String = joinpath(atmos.fastchem_work,"chemistry.dat")        # Chemistry by FastChem

        # Check file exists
        write_cfg = write_cfg || !isfile(confpath) || !isfile(elempath)

        # Write config, elements
        if write_cfg
            # Write config (fastchem is quite particular about the format)
            open(confpath,"w") do f
                write(f,"#Atmospheric profile input file \n")
                write(f,joinpath(atmos.fastchem_work,"pt.dat")*" \n\n")

                type_char = ["g","ce","cr"]
                write(f,"#Chemistry calculation type (gas phase only = g, equilibrium condensation = ce, rainout condensation = cr) \n")
                write(f,"$(type_char[chem_type]) \n\n")

                write(f,"#Chemistry output file \n")
                write(f,joinpath(atmos.fastchem_work,"chemistry.dat")*" "*joinpath(atmos.fastchem_work,"condensates.dat")*" \n\n")

                write(f,"#Monitor output file \n")
                write(f,joinpath(atmos.fastchem_work,"monitor.dat")*" \n\n")

                write(f,"#FastChem console verbose level (1 - 4); 1 = almost silent, 4 = detailed console output \n")
                write(f,"1 \n\n")

                write(f,"#Output mixing ratios (MR) or particle number densities (ND, default) \n")
                write(f,"ND \n\n")

                write(f,"#Element abundance file  \n")
                write(f,joinpath(atmos.fastchem_work,"elements.dat")*" \n\n")

                write(f,"#Species data files    \n")
                logK = joinpath(atmos.FC_DIR, "input/","logK/")
                write(f,joinpath(logK,"logK.dat")*" "*joinpath(logK,"logK_condensates.dat")*" \n\n")

                write(f,"#Accuracy of chemistry iteration \n")
                write(f,"1.0e-4 \n\n")

                write(f,"#Accuracy of element conservation \n")
                write(f,"1.0e-4 \n\n")

                write(f,"#Max number of chemistry iterations  \n")
                write(f,"50000 \n\n")

                write(f,"#Max number internal solver iterations  \n")
                write(f,"20000 \n\n")
            end

            # Calculate elemental abundances
            # number densities normalised relative to hydrogen
            # for each element X, value = log10(N_X/N_H) + 12
            # N = X(P/(K*T) , where X is the VMR and K is boltz-const
            N_t = zeros(Float64, length(phys.elements_list))      # total atoms in all gases
            N_g = zeros(Float64, length(phys.elements_list))      # atoms in current gas
            for gas in atmos.gas_names
                d = phys.count_atoms(gas)
                fill!(N_g, 0.0)
                for (i,e) in enumerate(phys.elements_list)
                    if e in keys(d)
                        N_g[i] += d[e]
                    end
                end
                # Get gas abundance from original VMR value, since the running
                #    one will be updated using FastChem's output. These will
                #    be normalised later in this function.
                N_g *= atmos.gas_ovmr[gas][atmos.nlev_c] * atmos.p[end] / (phys.k_B * atmos.tmp[end])  # gas contribution
                N_t += N_g  # add atoms in this gas to total atoms
            end

            # Write elemental abundances
            open(elempath,"w") do f
                write(f,"# Elemental abundances derived from AGNI volatiles \n")
                for (i,e) in enumerate(phys.elements_list)
                    if N_t[i] > 1.0e-30
                        # skip this element if its abundance is too small
                        # normalise relative to hydrogen
                        write(f, @sprintf("%s    %.3f \n",e,log10(N_t[i]/N_t[1]) + 12.0))
                        count_elem_nonzero += 1
                    end
                end
            end
        end

        # Write PT profile
        open(joinpath(atmos.fastchem_work,"pt.dat"),"w") do f
            write(f,"# AGNI temperature structure \n")
            write(f,"# bar, kelvin \n")
            for i in 1:atmos.nlev_c
                write(f,@sprintf("%.6e    %.6e \n", atmos.p[i]*1e-5, max(tmp_floor,atmos.tmp[i]) ))
            end
        end

        # Run fastchem
        run(pipeline(`$execpath $confpath`, stdout=devnull))

        # Check monitor output
        monitorpath::String = joinpath(atmos.fastchem_work,"monitor.dat")
        data = readdlm(monitorpath, '\t', String)
        fail_elem::String = ""
        fail_conv::String = ""
        for i in 1:atmos.nlev_c
            if data[i+1,6][1] == 'f'
                fail_elem *= @sprintf("%d ",i)
            end
            if data[i+1,5][1] == 'f'
                fail_conv *= @sprintf("%d ",i)
            end
        end
        if !isempty(fail_elem)
            @debug "Element conservation failed at levels  "*fail_elem
            state = 2
        end
        if !isempty(fail_conv)
            @debug "FastChem solver failed at levels  "*fail_conv
            if state == 2
                state = 4
            else
                state = 3
            end
        end

        # Get gas chemistry output
        if !isfile(chempath)
            @error "Could not find fastchem output"
            return 1
        end
        (data,head) = readdlm(chempath, '\t', Float64, header=true)
        data = transpose(data)  # convert to: gas, level

        # Clear VMRs
        for g in atmos.gas_names
            fill!(atmos.gas_vmr[g], 0.0)
        end

        # Parse gas chemistry
        g_fc::String = "_unset"
        d_fc::Dict = Dict{String, Int}()
        g_in::String = "_unset"
        match::Bool = false
        N_t = data[4,:] # at each level: sum of gas number densities

        for (i,h) in enumerate(head)  # for each column (gas)

            # skip T and P
            if i <= 5+count_elem_nonzero
                continue
            end

            # parse name
            g_fc = rstrip(lstrip(h))
            if occursin("_", g_fc)
                g_fc = split(g_fc, "_")[1]
            end
            g_fc = replace(g_fc, "cis"=>"", "trans"=>"")

            match = false

            # firstly, check if we have the FC name already stored
            for g in atmos.gas_names
                if atmos.gas_dat[g].fastchem_name == g_fc
                    match = true
                    g_in = g
                    break
                end
            end

            # not stored => search based on atoms
            if !match
                d_fc = phys.count_atoms(g_fc)  # get atoms dict from FC name

                for g in atmos.gas_names
                    if phys.same_atoms(d_fc, atmos.gas_dat[g].atoms)
                        match = true
                        g_in = g
                        break
                    end
                end
            end

            # matched?
            if match
                N_g = data[i,:]  # number densities for this gas
                atmos.gas_vmr[g_in][:] .+= N_g[:] ./ N_t[:]    # VMR for this gas
            end
        end

        # Do not renormalise mixing ratios, since this is done by fastchem
        # If we are missing gases then that's okay.

        # Find where we truncated the temperature profile,
        #      and make sure that regions above that use the same x_gas values
        for i in range(start=atmos.nlev_c, stop=1, step=-1)
            if atmos.tmp[i] < tmp_floor
                for g in atmos.gas_names
                    atmos.gas_vmr[g][1:i] .= atmos.gas_vmr[g][i+1]
                end
                break
            end
        end

        # See docstring for return codes
        return state
    end

    """
    **Normalise gas VMRs, keeping condensates unchanged**

    Only acts on a single model level.

    Parameters:
    - `atmos::atmosphere.Atmos_t`       atmosphere structure
    - `i::Int`                          level index to normalise
    """
    function normalise_vmrs!(atmos::atmosphere.Atmos_t, i::Int)
        # Work variables
        x_con::Float64 =    0.0
        x_dry::Float64 =    0.0
        x_dry_old::Float64= 0.0

        # Work out total VMR of all condensing gases
        for c in atmos.condensates
            if atmos.gas_sat[c][i]
                x_con += atmos.gas_vmr[c][i]
            end
        end

        # Calculate current and target dry VMRs
        x_dry = 1.0 - x_con
        x_dry_old = 0.0
        for g in atmos.gas_names
            # skip condensing gases, since their VMR is set by saturation
            if !atmos.gas_sat[g][i]
                x_dry_old += atmos.gas_vmr[g][i]
            end
        end

        # Renormalise VMR to unity, scaling DRY COMPONENTS ONLY
        for g in atmos.gas_names
            if !atmos.gas_sat[g][i]
                atmos.gas_vmr[g][i] *= x_dry / x_dry_old
            end
        end

        # Check total VMR at this level
        x_tot = 0.0
        for g in atmos.gas_names
            x_tot += atmos.gas_vmr[g][i]
        end
        if abs(x_tot - 1.0) > 1.0e-5
            @warn @sprintf("Mixing ratios sum to %.6e (level %d)",x_tot,i)
        end

        return nothing
    end

    """
    **Adjust gas VMRs according to saturation and cold-trap requirements**

    Volatiles which are allowed to condense are rained-out at condensing levels
    until the gas is exactly saturated, not supersaturated. If evaporation is enabled here,
    it will lead to enhanced mixing ratios at deeper levels as rain is converted back
    into gas from a liquid state.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    """
    function handle_saturation!(atmos::atmosphere.Atmos_t)

        # Single gas case does not apply here
        if atmos.gas_num == 1
            return
        end

        # Parameters
        evap_enabled::Bool =        false   # Enable re-vaporation of rain
        evap_efficiency::Float64 =  0.5     # Evaporation efficiency

        # Work arrays
        maxvmr::Dict{String, Float64} = Dict{String, Float64}() # max running VMR for each condensable
        cond_kg::Dict{String,Float64} = Dict{String, Float64}() # condensed kg/m2 for each condensable
        x_sat::Float64 = 0.0
        supcrit::Bool =  false

        # Set maximum value (for cold trapping)
        for c in atmos.condensates
            maxvmr[c] = atmos.gas_vmr[c][end]
        end

        # Reset mixing ratios to surface values
        # Reset phase change flags
        # Reset condensation yield values
        for g in atmos.gas_names
           atmos.gas_vmr[g][1:end-1] .= atmos.gas_vmr[g][end]
           atmos.gas_sat[g][:] .= false
           atmos.gas_yield[g][:] .= 0.0
        end

        # Reset water cloud
        fill!(atmos.cloud_arr_r, 0.0)
        fill!(atmos.cloud_arr_l, 0.0)
        fill!(atmos.cloud_arr_f, 0.0)

        # Handle condensation
        for i in range(start=atmos.nlev_c-1, stop=1, step=-1)

            # For each condensate
            for c in atmos.condensates

                # Reset condensation and rain
                cond_kg[c] = 0.0   # kg/m2 of 'c' condensate produced at this level

                # check criticality
                supcrit = atmos.tmp[i] > atmos.gas_dat[c].T_crit+1.0e-5

                # saturation mixing ratio
                x_sat = phys.get_Psat(atmos.gas_dat[c], atmos.tmp[i]) / atmos.p[i]

                # cold trap
                if atmos.gas_vmr[c][i] > maxvmr[c]
                    atmos.gas_vmr[c][i] = maxvmr[c]
                    atmos.gas_sat[c][i] = true
                end

                # condense if supersaturated
                if (atmos.gas_vmr[c][i] > x_sat) && !supcrit

                    # set rainout kg/m2
                    cond_kg[c] = atmos.gas_dat[c].mmw*atmos.p[i]*
                                            (atmos.gas_vmr[c][i] - x_sat)/
                                            (atmos.layer_grav[i] * atmos.layer_mmw[i])

                    # condensation yield at this level
                    atmos.gas_yield[c][i] += cond_kg[c]

                    # set new vmr
                    atmos.gas_vmr[c][i] = x_sat

                    # store vmr for cold trapping at levels above this one
                    maxvmr[c] = x_sat

                    # flag condensate as actively condensing at this level
                    atmos.gas_sat[c][i] = true
                    # @printf("%d: %s rain generated %.3e \n", i, c, cond_kg[c])

                    # Set water cloud at this level
                    if c == "H2O"
                        # mass mixing ratio (take ratio of mass surface densities [kg/m^2])
                        atmos.cloud_arr_l[i] = (cond_kg["H2O"] * atmos.cond_alpha) /
                                                    atmos.layer_mass[i]

                        if atmos.cloud_arr_l[i] > 1.0
                            @warn "Water cloud mass mixing ratio is greater than unity (level $i)"
                        end

                        # droplet radius and area fraction (fixed values)
                        atmos.cloud_arr_r[i] = atmos.cloud_val_r
                        atmos.cloud_arr_f[i] = atmos.cloud_val_f
                    end
                end # end saturation check

            end # end condensate

            normalise_vmrs!(atmos, i)
        end # end i levels

        # recalculate layer properties
        calc_layer_props!(atmos)

        # Evaporate rain within unsaturated layers
        if evap_enabled

            total_rain::Float64 = 0.0
            i_top_dry::Int = 1
            i_bot_dry::Int = atmos.nlev_c

            # For each condensable
            for c in atmos.condensates

                # reset dry region
                i_top_dry = 1
                i_bot_dry = atmos.nlev_c

                # accumulate rain
                total_rain = sum(atmos.gas_yield[c])

                # no rain? go to next condensable
                if total_rain < 1.0e-10
                    continue
                end

                # locate top of dry region
                for j in 1:atmos.nlev_c-1
                    if atmos.gas_sat[c][j]
                        # this layer is saturated => dry region must be below it
                        i_top_dry = j+1
                    end
                end

                # locate bottom of dry region
                for j in range(start=atmos.nlev_c, stop=i_top_dry+1, step=-1)
                    if atmos.tmp[j] > atmos.gas_dat[c].T_crit
                        # this layer is supercritical => dry region must be above it
                        i_bot_dry = j-1
                    end
                end
                if i_bot_dry < i_top_dry
                    i_bot_dry = i_top_dry
                end

                # evaporate rain in dry region
                for j in i_top_dry+1:i_bot_dry

                    # evaporate up to saturation
                    dp_sat = phys.get_Psat(atmos.gas_dat[c], atmos.tmp[j]) -
                                                    atmos.gas_vmr[c][j]*atmos.p[j]

                    # Calculate kg/m2 of gas that would saturate layer j
                    dm_sat = atmos.gas_dat[c].mmw * dp_sat/
                                            (atmos.layer_grav[j] * atmos.layer_mmw[j])

                    # can we evaporate all rain within this layer?
                    if total_rain < dm_sat
                        # yes, so don't evaporate more rain than the total
                        dm_sat = total_rain
                    end
                    atmos.gas_sat[c][j] = true

                    # Evaporation efficiency factor
                    #   This fraction of the rain that *could* be evaporated
                    #   at this layer *is* converted to vapour in this layer.
                    dm_sat *= evap_efficiency

                    # offset condensate yield at this level by the evaporation
                    atmos.gas_yield[c][j] -= dm_sat

                    # add partial pressure from evaporation to pp at this level
                    dp_sat = dm_sat * atmos.layer_grav[j] *
                                              atmos.layer_mmw[j] / atmos.gas_dat[c].mmw

                    # convert extra pp to extra vmr
                    atmos.gas_vmr[c][j] += dp_sat / atmos.p[j]

                    # reduce total rain correspondingly
                    total_rain -= dm_sat

                    # @printf("    %d: evaporating %.3e \n", j, dm_sat)

                    # ---------------------
                    # WARNING
                    # THIS IS NOT CORRECT. IT WILL LEAD TO sum(x_gas)>1
                    # IN MANY CASES BECAUSE OF THE PREDEFINED PRESSURE GRID.
                    # THIS NEEDS TO RE-ADJUST THE WHOLE PRESSURE GRID IN
                    # ORDER TO CONSERVE MASS!
                    # ----------------------
                    # normalise_vmrs!(atmos, j)
                    #

                    # Recalculate layer mmw
                    atmos.layer_mmw[j] = 0.0
                    for g in atmos.gas_names
                        atmos.layer_mmw[j] += atmos.gas_vmr[g][j] * atmos.gas_dat[g].mmw
                    end

                end # go to next j level (below)

            end # end loop over condensates

        end # end evaporation scheme

        # Recalculate layer properties
        calc_layer_props!(atmos)

        return nothing
    end

    # Set cloud properties within condensing regions.
    function water_cloud!(atmos::atmosphere.Atmos_t)

        # Reset
        fill!(atmos.cloud_arr_r, 0.0)
        fill!(atmos.cloud_arr_l, 0.0)
        fill!(atmos.cloud_arr_f, 0.0)

        # Get index of water
        if !("H2O" in atmos.gas_names)
            return nothing
        end

        # Set level-by-level
        x::Float64      = 0.0
        for i in 1:atmos.nlev_c-1

            # Water VMR
            x = atmos.gas_vmr["H2O"][i]
            if x < 1.0e-10
                continue
            end

            # Use mask from atmos struct
            if atmos.gas_sat["H2O"][i]
                atmos.cloud_arr_r[i] = atmos.cloud_val_r
                atmos.cloud_arr_l[i] = atmos.cloud_val_l
                atmos.cloud_arr_f[i] = atmos.cloud_val_f
            end
        end

        return nothing
    end

    # Smooth temperature at cell-centres
    # function smooth_centres!(atmos::atmosphere.Atmos_t, width::Int)

    #     if width > 2
    #         if mod(width,2) == 0
    #             width += 1 # window width has to be an odd number
    #         end
    #         atmos.tmp = moving_average.hma(atmos.tmp, width)
    #     end
    #     clamp!(atmos.tmp, atmos.tmp_floor, atmos.tmp_ceiling)

    #     return nothing
    # end

    """
    **Set cell-edge temperatures from cell-centre values.**

    Uses interpolation within the bulk of the column and extrapolation for the
    topmost edge.

    Arguments:
    - `atmos::Atmos_t`            the atmosphere struct instance to be used.
    - `back_interp::Bool=false`   interpolate resultant tmpl back to tmp (should be avoided)
    """
    function set_tmpl_from_tmp!(atmos::atmosphere.Atmos_t; back_interp::Bool=false)

        # Interpolate temperature to bulk cell-edge values (log-linear)
        itp = Interpolator(log.(atmos.p), atmos.tmp)
        atmos.tmpl[2:end-1] .= itp.(log.(atmos.pl[2:end-1]))

        # Extrapolate top edge temperature (log-linear)
        grad_dt::Float64 = atmos.tmp[1] - atmos.tmp[2]
        grad_dp::Float64 = log(atmos.p[1]/atmos.p[2])
        atmos.tmpl[1] = atmos.tmp[1] + grad_dt/grad_dp * log(atmos.pl[1]/atmos.p[1])

        # Set bottom edge to bottom cell-centre value
        # This is fine because the bottom cell is very small (in pressure space)
        atmos.tmpl[end]=atmos.tmp[end]

        # Clamp
        clamp!(atmos.tmpl, atmos.tmp_floor, atmos.tmp_ceiling)

        # Second interpolation back to cell-centres.
        # This can help prevent grid-imprinting issues, but in some cases it can
        # lead to unphysical behaviour, because it gives the interpolator too
        # much control over the temperature profile at small scales.
        if back_interp
            clamp!(atmos.tmpl, atmos.tmp_floor, atmos.tmp_ceiling)
            itp = Interpolator(log.(atmos.pl), atmos.tmpl)
            atmos.tmp[:] .= itp.(log.(atmos.p[:]))
        end

        return nothing
    end

    # Get interleaved cell-centre and cell-edge PTZ grid
    function get_interleaved_ptz(atmos::atmosphere.Atmos_t)
        arr_Z::Array{Float64, 1} = zeros(Float64, atmos.nlev_c + atmos.nlev_l)
        arr_T::Array{Float64, 1} = zeros(Float64, atmos.nlev_c + atmos.nlev_l)
        arr_P::Array{Float64, 1} = zeros(Float64, atmos.nlev_c + atmos.nlev_l)
        idx::Int = 1

        # top
        arr_Z[1] = atmos.zl[1]
        arr_T[1] = atmos.tmpl[1]
        arr_P[1] = atmos.pl[1]

        # middle
        for i in 1:atmos.nlev_c
            idx = (i-1)*2
            arr_T[idx+1] = atmos.tmpl[i]
            arr_T[idx+2] = atmos.tmp[i]

            arr_P[idx+1] = atmos.pl[i]
            arr_P[idx+2] = atmos.p[i]

            arr_Z[idx+1] = atmos.zl[i]
            arr_Z[idx+2] = atmos.z[i]
        end

        # bottom
        arr_Z[end] = atmos.zl[end]
        arr_T[end] = atmos.tmpl[end]
        arr_P[end] = atmos.pl[end]

        return arr_P, arr_T, arr_Z
    end

end
