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
    import ..moving_average
    import ..phys
    import ..spectrum

    
    # Solution types
    SOL_TYPES::Array{String, 1} = [
        "steadystate_ext",  # zero divergence at fixed tmp_surf, extrapolated tmpl[end]
        "cond_skin",        # set tmp_surf such that conductive skin conserves energy flux
        "flux_eff",         # total flux at each level equal to flux_eff
        "target_olr",       # OLR is equal to target_olr
    ]

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
        ROOT_DIR::String
        OUT_DIR::String

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
        albedo_s::Float64               # Surface albedo
        zenith_degrees::Float64         # Solar zenith angle [deg]
        toa_heating::Float64            # Derived downward shortwave radiation flux at topmost level [W m-2]
        instellation::Float64           # Solar flux at top of atmopshere [W m-2]
        s0_fact::Float64                # Scale factor to instellation (cronin+14)
        tmp_surf::Float64               # Surface brightness temperature [K]
        tmp_eff::Float64                # Effective temperature of the planet [K]
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
        rp::Float64         # planet radius [m]

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
        gas_vmr::Dict{String, Array{Float64,1}}     # Layer mole fractions in dict format, (key,value) = (gas_name,array)
        gas_phase::Dict{String, Array{Bool, 1}}     # Layer condensation flag in dict format, (key,value) = (gas_name,array)
        gas_dat::Dict{String, phys.Gas_t}           # struct variables containing thermodynamic data
        gas_cprod::Dict{String, Array{Float64,1}}   # condensate yield at each level (can be negative, representing evaporation)
        condensates::Array{String, 1}                   # List of condensing gases (strings)
        single_component::Bool                          # Does a single gas make up 100% of layer at any point in the column?

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
        flux_eff::Float64                   # Effective flux  [W m-2] for sol_type=3
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
        mask_c::Array{Int,1}                # Layers which are (or were recently) convective (value is set to >0)
        mask_decay::Int                     # How long is 'recent'
        flux_cdry::Array{Float64,1}         # Dry convective fluxes from MLT
        Kzz::Array{Float64,1}               # Eddy diffusion coefficient from MLT

        # Conduction 
        flux_cdct::Array{Float64,1}         # Conductive flux [W m-2]

        # Cloud and condensation
        flux_l::Array{Float64, 1}           # Latent heat energy flux [W m-2]
        mask_l::Array{Int,1}                # Layers which are (or were recently) condensing
        #    These arrays give the cloud properties within layers containing cloud
        cloud_arr_r::Array{Float64,1}       # Characteristic dimensions of condensed species [m]. 
        cloud_arr_l::Array{Float64,1}       # Mass mixing ratios of condensate. 0 : saturated water vapor does not turn liquid ; 1 : the entire mass of the cell contributes to the cloud
        cloud_arr_f::Array{Float64,1}       # Total cloud area fraction in layers. 0 : clear sky cell ; 1 : the cloud takes over the entire area of the Cell
        #    These floats give the default values that the above arrays are filled with upon cloud formation
        cloud_val_r::Float64 
        cloud_val_l::Float64 
        cloud_val_f::Float64 

        # Cell-internal heating 
        ediv_add::Array{Float64, 1}     # Additional energy dissipation inside each cell [W m-3] (e.g. from advection)

        # Total energy flux
        flux_dif::Array{Float64,1}      # Flux lost at each level [W m-2]
        flux_tot::Array{Float64,1}      # Total upward-directed flux at cell edges [W m-2]

        # Heating rate felt at each level [K/day]
        heating_rate::Array{Float64,1} 

        # FastChem stuff 
        fastchem_flag::Bool             # Fastchem enabled?
        fastchem_path::String           # Path to fastchem installation folder (incl. executable)
        fastchem_work::String           # Path to fastchem dir inside AGNI output folder
        
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
    - `s0_fact::Float64`                scale factor applied to instellation to account for planetary rotation (i.e. S_0^*/S_0 in Cronin+14)    
    - `albedo_b::Float64`               bond albedo scale factor applied to instellation in order to imitate shortwave reflection     
    - `zenith_degrees::Float64`         angle of radiation from the star, relative to the zenith [degrees].    
    - `tmp_surf::Float64`               effective surface temperature to provide upward longwave flux at the bottom of the atmosphere [K].    
    - `gravity::Float64`                gravitational acceleration at the surface [m s-2].    
    - `radius::Float64`                 planet radius at the surface [m].    
    - `nlev_centre::Int`                number of model levels.    
    - `p_surf::Float64`                 total surface pressure [bar].    
    - `p_top::Float64`                  total top of atmosphere pressure [bar].    
    - `mf_dict=nothing`                 dictionary of mole fractions in the format (key,value)=(gas,mf).    
    - `mf_path=nothing`                 path to file containing mole fractions at each level. 
    - `condensates::Array{String,1}`    list of condensates (names)

    Optional arguments:
    - `albedo_s::Float64`               surface albedo.   
    - `tmp_floor::Float64`              temperature floor [K].   
    - `C_d::Float64`                    turbulent heat exchange coefficient [dimensionless].   
    - `U::Float64`                      surface wind speed [m s-1].   
    - `tmp_magma::Float64`              mantle temperature [K] for sol_type==2.   
    - `skin_d::Float64`                 skin thickness [m].   
    - `skin_k::Float64`                 skin thermal conductivity [W m-1 K-1].   
    - `overlap_method::Int`             gaseous overlap scheme (2: rand overlap, 4: equiv extinct, 8: ro+resort+rebin).   
    - `target_olr::Float64`             target OLR [W m-2] for sol_type==4.   
    - `tmp_eff::Float64`                planet's effective brightness temperature [K] for sol_type==3.   
    - `all_channels::Bool`              use all channels available for RT?   
    - `flag_rayleigh::Bool`             include rayleigh scattering?   
    - `flag_gcontinuum::Bool`           include generalised continuum absorption?   
    - `flag_continuum::Bool`            include continuum absorption?   
    - `flag_aerosol::Bool`              include aersols?   
    - `flag_cloud::Bool`                include clouds?   
    - `thermo_functions::Bool`          use temperature-dependent thermodynamic properties   
    - `fastchem_path::String`           path to FastChem folder (empty string => disabled)   

    Returns:
        Nothing
    """
    function setup!(atmos::atmosphere.Atmos_t, 
                    ROOT_DIR::String, OUT_DIR::String, 
                    spfile::String, 
                    instellation::Float64, s0_fact::Float64, albedo_b::Float64, zenith_degrees::Float64,
                    tmp_surf::Float64,
                    gravity::Float64, radius::Float64,
                    nlev_centre::Int, p_surf::Float64, p_top::Float64;
                    mf_dict=                    nothing,
                    mf_path =                   nothing,
                    condensates::Array{String,1} = String[],
                    albedo_s::Float64 =         0.0,
                    tmp_floor::Float64 =        80.0,
                    C_d::Float64 =              0.001,
                    U::Float64 =                2.0,
                    tmp_magma::Float64 =        3000.0,
                    skin_d::Float64 =           0.05,
                    skin_k::Float64 =           2.0,
                    overlap_method::Int =       4,
                    target_olr::Float64 =       0.0,
                    tmp_eff::Float64 =          0.0,
                    all_channels::Bool  =       true,
                    flag_rayleigh::Bool =       false,
                    flag_gcontinuum::Bool =     false,
                    flag_continuum::Bool =      false,
                    flag_aerosol::Bool =        false,
                    flag_cloud::Bool =          false,
                    thermo_functions::Bool =    true,
                    fastchem_path::String =     ""
                    )

        if !isdir(OUT_DIR) && !isfile(OUT_DIR)
            mkdir(OUT_DIR)
        end

        # Code versions 
        atmos.AGNI_VERSION = "0.4.0"
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
        atmos.tmp_eff =         tmp_eff
        atmos.grav_surf =       max(1.0e-7, gravity)
        atmos.zenith_degrees =  max(min(zenith_degrees,89.8), 0.2)
        atmos.albedo_s =        max(min(albedo_s, 1.0 ), 0.0)
        atmos.instellation =    max(instellation, 0.0)
        atmos.albedo_b =        max(min(albedo_b,1.0), 0.0)
        atmos.s0_fact =         max(s0_fact,0.0)
        atmos.toa_heating =     atmos.instellation * (1.0 - atmos.albedo_b) * s0_fact * cosd(atmos.zenith_degrees)

        atmos.flux_eff =        phys.sigma * (atmos.tmp_eff)^4.0  
        atmos.target_olr =      max(1.0e-20,target_olr)
        
        atmos.mask_decay =      15 
        atmos.C_d =             max(0,C_d)
        atmos.U =               max(0,U)

        atmos.tmp_magma =       max(atmos.tmp_floor, tmp_magma)
        atmos.skin_d =          max(1.0e-6,skin_d)
        atmos.skin_k =          max(1.0e-6,skin_k)

        if p_top > p_surf 
            error("p_top must be less than p_surf")
        end

        atmos.p_toa =           p_top * 1.0e+5 # Convert bar -> Pa
        atmos.p_boa =           p_surf * 1.0e+5
        atmos.rp =              max(1.0, radius)

        # absorption contributors
        atmos.control.l_gas =           true
        atmos.control.l_rayleigh =      flag_rayleigh
        atmos.control.l_continuum =     flag_continuum
        atmos.control.l_cont_gen =      flag_gcontinuum
        atmos.control.l_aerosol =       flag_aerosol
        atmos.control.l_cloud =         flag_cloud
        atmos.control.l_drop =          flag_cloud
        atmos.control.l_ice  =          false

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

        # Hardcoded cloud properties 
        atmos.cloud_val_r   = 1.0e-5  # 10 micron droplets
        atmos.cloud_val_l   = 0.8     # 80% of the saturated vapor turns into cloud
        atmos.cloud_val_f   = 0.8     # 80% of the cell "area" is cloud

        # Read mole fractions
        if isnothing(mf_dict) && isnothing(mf_path)
            error("No mole fractions provided")
        end

        if !isnothing(mf_dict) && !isnothing(mf_path)
            error("Mole fractions provided twice")
        end

        mf_source::Int = 0  # source for mf (0: dict, 1: file)
        if isnothing(mf_dict) && !isnothing(mf_path)
            mf_source = 1
        end
        
        # The values will be stored in a dict of arrays
        atmos.gas_names =   Array{String}(undef, 0)           # list of names (String)
        atmos.gas_dat =     Dict{String, phys.Gas_t}()        # dict of data structures (phys.Gas_t)
        atmos.gas_vmr  =    Dict{String, Array{Float64,1}}()  # dict of VMR arrays (Float)
        atmos.gas_phase =   Dict{String, Array{Bool, 1}}()    # dict for phase change (Bool) 
        atmos.gas_cprod =   Dict{String, Array{Float64,1}}()  # dict of condensate yield values (Float [kg])
        atmos.gas_num   =   0                                 # number of gases 
        atmos.condensates   =   Array{String}(undef, 0)           # list of condensates (String)

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
                error("Could not read file '$mf_path'")
            end
            @info "Composition set by file"

            # get header
            mf_head = readline(abspath(mf_path))
            mf_head = mf_head[2:end]  # remove comment symbol at start
            mf_head = replace(mf_head, " " => "")  # remove whitespace
            heads   = split(mf_head, ",")[4:end] # split by column and drop first three

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
            mf_body = readdlm(abspath(mf_path), ',', Float64; header=false, skipstart=2)
            mf_body = transpose(mf_body)

            # set composition by interpolating with pressure array 
            gidx::Int=0
            for li in 4:lastindex(heads)
                # Gas index 
                gidx += 1

                # Arrays from file 
                arr_p = mf_body[1,:]
                arr_x = mf_body[li,:]

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
                itp = Interpolator(arr_p, arr_x)
                
                # Set values in atmos struct 
                for i in 1:atmos.nlev_c
                    atmos.gas_vmr[atmos.gas_names[gidx]][i] = itp(atmos.p[i])
                end 
            end 

        end # end read VMR from file 

        # store condensates 
        for c in condensates
            push!(atmos.condensates, c)
        end 

        # Cannot have n_gas==n_cond AND n_cond>1, because it will overspecify
        #   the total pressure within condensing regions
        if (length(atmos.condensates) == atmos.gas_soc_num) && (length(atmos.condensates)>1)
            error("There must be at least one non-condensable gas")
            return 
        end 

        # Validate condensate names 
        if length(condensates) > 0
            for c in condensates
                if !(c in atmos.gas_names)
                    @error "Invalid condensate '$c'"
                    return false
                end 
            end
        end 

        # set condensation mask and yield values [kg]
        for g in atmos.gas_names 
            atmos.gas_phase[g] = falses(atmos.nlev_c)
            atmos.gas_cprod[g] = zeros(Float64, atmos.nlev_c)
        end 

        # Check that we actually stored some values
        if atmos.gas_num == 0
            error("No mole fractions were stored")
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

        # Check if atmosphere is (in reality) composed of a single gas at 
        #    any point, since this has implications for the condensation scheme
        atmos.single_component = false
        for i in 1:atmos.nlev_c 
            for g in atmos.gas_names 
                atmos.single_component = atmos.single_component || (atmos.gas_vmr[g][i]>1.0-1.0e-10)
            end 
        end 

        # Load gas thermodynamic data 
        for g in atmos.gas_names 
            atmos.gas_dat[g] = phys.load_gas(g, atmos.thermo_funct)
            if (g in atmos.condensates) && (atmos.gas_dat[g].stub)
                error("No thermodynamic data found for condensable gas $g")
            end
        end 

        # Set initial temperature profile to a small value which still keeps
        #   all of the gases supercritical. This should be a safe condition to 
        #   default to, although the user should specify a profile in the cfg.
        for g in atmos.gas_names 
            atmos.tmpl[end] = max(atmos.tmpl[end], atmos.gas_dat[g].T_crit+20.0)
            fill!(atmos.tmpl, atmos.tmpl[end])
            fill!(atmos.tmp, atmos.tmpl[end])
        end 

        # Fastchem 
        atmos.fastchem_flag = false 
        if !isempty(fastchem_path)
            # check folder
            atmos.fastchem_path = abspath(fastchem_path)
            if !isdir(fastchem_path)
                @error "Could not find fastchem folder at '$(fastchem_path)'"
            end 
            atmos.fastchem_work = joinpath(atmos.OUT_DIR, "fastchem/")
            
            # check executable 
            atmos.fastchem_flag = isfile(joinpath(fastchem_path,"fastchem")) 
            if !atmos.fastchem_flag 
                @error "Could not find fastchem executable inside '$(atmos.fastchem_path)' "
            else 
                @info "Found FastChem executable"
            end 
        end 

        # Record that the parameters are set
        atmos.is_param = true
        atmos.is_solved = false 
        atmos.is_converged = false

        @debug "Setup complete"
        return nothing
    end # end function setup 

    """
    **Get the mole fraction of a gas within the atmosphere.**

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    - `gas::String`             name of the gas (e.g. "H2O").
    - `lvl::Int`                model level to measure mole fraction

    Returns:
    - `x::Float64`              mole fraction of `gas`.
    """
    function get_x(atmos::atmosphere.Atmos_t, gas::String, lvl::Int)::Float64

        if gas in atmos.gas_names
            if (lvl >= 1) && (lvl <= atmos.nlev_c)
                return atmos.gas_vmr[gas][lvl]
            else
                error("Invalid level provided ($lvl)") 
            end 
        else
            @error "Invalid gas $gas queried in get_x()" 
            return 0.0
        end 
    end

    """
    **Calculate properties within each layer of the atmosphere (e.g. density, mmw).**

    Assumes that the atmosphere may be treated as an ideal gas.

    Arguments:
        - `atmos::Atmos_t`          the atmosphere struct instance to be used.

    Returns:
        - `ok::Bool`                function result is ok
    """
    function calc_layer_props!(atmos::atmosphere.Atmos_t)::Bool
        if !atmos.is_param
            error("atmosphere parameters have not been set")
        end

        ok::Bool = true

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
        g1::Float64 = 0.0
        g2::Float64 = 0.0
        dzc::Float64= 0.0
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
            dzc = phys.R_gas * atmos.tmp[i] / (atmos.layer_mmw[i] * g1 * atmos.p[i]) * (atmos.pl[i+1] - atmos.p[i]) 
            if (dzc < 1e-20)
                @error "Height integration resulted in dz <= 0 at level $i (l -> c)"
                ok = false 
            end
            if  (dzc > 1e9)
                @error "Height integration blew up at level $i (l -> c)"
                ok = false 
                dzc = 1e9
            end
            atmos.z[i] = atmos.zl[i+1] + dzc

            # Integrate from centre to upper edge
            g2 = GMpl / ((atmos.rp + atmos.z[i])^2.0)
            dzl = phys.R_gas * atmos.tmp[i] / (atmos.layer_mmw[i] * g2 * atmos.p[i]) * (atmos.p[i]- atmos.pl[i]) 
            if (dzl < 1e-20)
                @error "Height integration resulted in dz <= 0 at level $i (c -> l)"
                ok = false 
            end
            if (dzl > 1e9)
                @error "Height integration blew up at level $i (c -> l)"
                ok = false 
                dzl = 1e9
            end
            atmos.zl[i] = atmos.z[i] + dzl
            
            # Layer average gravity [m s-2]
            atmos.layer_grav[i] = GMpl / ((atmos.rp + atmos.z[i])^2.0)

            # Layer geometrical thickness [m]
            atmos.layer_thick[i] = atmos.zl[i] - atmos.zl[i+1]
        end 

        # Mass (technically area density [kg m-2]) and density [kg m-3]
        for i in 1:atmos.atm.n_layer
            atmos.layer_mass[i] = (atmos.atm.p_level[1, i] - atmos.atm.p_level[1, i-1])/atmos.layer_grav[i]
            atmos.atm.mass[1, i] = atmos.layer_mass[i]          # pass to SOCRATES

            atmos.layer_density[i] = ( atmos.p[i] * atmos.layer_mmw[i] )  / (phys.R_gas * atmos.tmp[i]) 
            atmos.atm.density[1,i] = atmos.layer_density[i]     # pass to SOCRATES
        end

        return ok
    end 

    """
    **Generate pressure grid.**

    Equally log-spaced between p_boa and p_boa. The topmost and bottommost layers 
    have bespoke relative thicknesses of 1±boundary_scale, in order to avoid 
    numerical instabilities.
    
    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `boundary_scale::Float64=1.0e-4`  scale factor for the thickness of the bottom-most cell.
    """
    function generate_pgrid!(atmos::atmosphere.Atmos_t; boundary_scale::Float64=1.0e-4)

        boundary_scale = max(min(boundary_scale, 1.0-1.0e-8), 1.0e-8)

        # Logarithmically-spaced levels
        atmos.pl           = zeros(Float64, atmos.nlev_l)
        atmos.pl[end]      = atmos.p_boa 
        atmos.pl[1:end-1] .= 10 .^ range( log10(atmos.p_toa), stop=log10(atmos.p_boa*(1.0-boundary_scale)), length=atmos.nlev_l-1)

        # Set pressure cell-centre array using geometric mean
        atmos.p = zeros(Float64, atmos.nlev_c)
        atmos.p[1:end] .= sqrt.(atmos.pl[1:end-1].*atmos.pl[2:end])

        return nothing
    end

    
    """
    **Allocate atmosphere arrays, prepare spectral files, and final steps before RT calculations.**

    Arguments:
    - `atmos::Atmos_t`                 the atmosphere struct instance to be used.
    - `stellar_spectrum::String`       path to stellar spectrum csv file (will not modify spectral file if this is left blank)
    - `one_gas::String`                only include this gas in the RT (all gases used if blank)
    """
    function allocate!(atmos::atmosphere.Atmos_t, stellar_spectrum::String; one_gas::String="")

        @debug "Allocate atmosphere"
        if !atmos.is_param
            error("atmosphere parameters have not been set")
        end

        atmos.atm.n_profile = 0

        #########################################
        # spectral data
        #########################################

        # Validate files
        if !isfile(atmos.spectral_file)
            error("Source spectral file '$(atmos.spectral_file)' does not exist")
        end

        spectral_file_run::String  = joinpath([atmos.OUT_DIR, "runtime.sf"])  
        spectral_file_runk::String = joinpath([atmos.OUT_DIR, "runtime.sf_k"])  

        # Setup spectral file 
        if !isempty(stellar_spectrum)

            if !isfile(stellar_spectrum)
                error("Stellar spectrum file '$(stellar_spectrum)' does not exist")
            end

            # File to be loaded
            rm(spectral_file_run , force=true)
            rm(spectral_file_runk, force=true)

            atmos.star_file = abspath(stellar_spectrum)

            # Write spectrum in required format 
            socstar = joinpath([atmos.OUT_DIR, "socstar.dat"])  
            wl, fl = spectrum.load_from_file(atmos.star_file)
            spectrum.write_to_socrates_format(wl, fl, socstar)

            # Insert stellar spectrum and rayleigh scattering, if required
            spectrum.insert_stellar_and_rscatter(atmos.spectral_file, socstar, spectral_file_run, atmos.control.l_rayleigh)

        else 
            # Stellar spectrum was not provided, which is taken to mean that the spectral file includes it already 
            atmos.star_file = "IN_SPECTRAL_FILE"
            spectral_file_run  = atmos.spectral_file
            spectral_file_runk = atmos.spectral_file*"_k"
            @info "Using pre-existing spectral file without modifications"
        end
     
        # Validate files
        if occursin("NaN", readchomp(spectral_file_runk))
            error("Spectral_k file contains NaN values")
        end 
        if occursin("NaN", readchomp(spectral_file_run))
            error("Spectral file contains NaN values")
        end 

        # Read-in spectral file to be used at runtime
        atmos.control.spectral_file = spectral_file_run

        if isempty(one_gas)
            SOCRATES.set_spectrum(spectrum=atmos.spectrum, 
                                spectral_file=atmos.control.spectral_file, 
                                l_all_gasses=true)
        else   
            # set one to true 
            kw_dict = Dict(eval(Meta.parse(":l_"*lowercase(one_gas))) => true)

            # call set_spectrum
            SOCRATES.set_spectrum(;spectrum=atmos.spectrum, 
                                spectral_file=atmos.control.spectral_file, 
                                l_all_gasses=false, kw_dict...)
        end 


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
        npd_direction = 1           # Maximum number of directions for radiances
        npd_layer = atmos.nlev_c    # Number of layers

        npd_column = 24 # Maximum number of cloudy subcolumns
        npd_profile = 1
        npd_max_order = 101 #       Maximum order of spherical harmonics used
        npd_brdf_basis_fnc = 2 #       Number of BRDF basis functions
        npd_brdf_trunc = 5 #       Order of BRDF truncation
        npd_profile_aerosol_prsc = 9 # Size allocated for profiles of prescribed aerosol optical properties
        npd_profile_cloud_prsc = 9 # Size allocated for profiles of prescribed cloudy optical properties
        npd_opt_level_aerosol_prsc = 170 # Size allocated for levels of prescribed aerosol optical properties
        npd_opt_level_cloud_prsc = 170 # Size allocated for levels of prescribed cloudy optical properties

        # modules_gen/dimensioms_fixed_pcf.f90
        npd_cloud_component        =  4 #   Number of components of clouds.
        npd_cloud_type             =  4 #   Number of permitted types of clouds.
        npd_overlap_coeff          = 18 #   Number of overlap coefficients for cloud
        npd_source_coeff           =  2 #   Number of coefficients for two-stream sources
        npd_region                 =  2 # Number of regions in a layer

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
        
        SOCRATES.allocate_atm(  atmos.atm,   atmos.dimen, atmos.spectrum)
        SOCRATES.allocate_cld(  atmos.cld,   atmos.dimen, atmos.spectrum)
        SOCRATES.allocate_aer(  atmos.aer,   atmos.dimen, atmos.spectrum)
        SOCRATES.allocate_bound(atmos.bound, atmos.dimen, atmos.spectrum)

        atmos.bound.rho_alb[:, SOCRATES.rad_pcf.ip_surf_alb_diff, :] .= atmos.albedo_s   #   Set surface albedo.

        # atm sizes and coordinates 
        atmos.atm.n_layer = npd_layer
        atmos.atm.n_profile = 1
        atmos.atm.lat[1] = 0.0
        atmos.atm.lon[1] = 0.0

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
            error("n_channel $n_channel != 1 and != $n_band_active not supported ")
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
            Bool(atmos.spectrum.Basic.l_present[3]) ||
                error("The spectral file contains no rayleigh scattering data.")
        end
        
        if atmos.control.l_aerosol
            Bool(atmos.spectrum.Basic.l_present[11]) ||
                error("The spectral file contains no aerosol data.")
        end

        if atmos.control.l_gas
            Bool(atmos.spectrum.Basic.l_present[5]) ||
                error("The spectral file contains no gaseous absorption data.")
        end

        if atmos.control.l_continuum
            Bool(atmos.spectrum.Basic.l_present[9]) ||
                error("The spectral file contains no continuum absorption data.")
        end

        if atmos.control.l_cont_gen
            Bool(atmos.spectrum.Basic.l_present[19]) ||
                error("The spectral file contains no generalised continuum absorption data.")
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
            error("Invalid overlap method")
        end

        for j in atmos.control.first_band:atmos.control.last_band
            atmos.control.i_gas_overlap_band[j] = atmos.control.i_gas_overlap
        end

        # Check supported gases
        atmos.gas_soc_num       = atmos.spectrum.Gas.n_absorb               # number of gases 
        atmos.gas_soc_names     = Array{String}(undef, atmos.gas_soc_num)   # list of names   
        for i_gas in 1:atmos.gas_soc_num   # for each supported gas
            atmos.gas_soc_names[i_gas] = SOCRATES.input_head_pcf.header_gas[atmos.spectrum.Gas.type_absorb[i_gas]]
        end 

        # VMRs are provided to SOCRATES when radtrans is called
        # For now, they are just stored inside the atmos struct

        # Print info on the gases
        @info "Allocated atmosphere with composition"
        gas_flags::String = ""
        for g in atmos.gas_names 
            gas_flags = ""
            if !(g in atmos.gas_soc_names)  # flag as included in radtrans
                gas_flags *= "NO_OPACITY "
            end
            if g in atmos.condensates       # flag as condensable 
                gas_flags *= "CNDSBLE "
            end 
            if atmos.gas_dat[g].stub        # flag as containing stub thermo data
                gas_flags *= "NO_THERMO "
            end 
            if !isempty(gas_flags)
                gas_flags = "($(gas_flags[1:end-1]))"
            end 
            @info @sprintf("    %6.2e %-5s %s", atmos.gas_vmr[g][end], g, gas_flags)
        end 
 

        # Calc layer properties using initial temperature profile 
        calc_layer_props!(atmos)

        ################################
        # Aerosol processes
        #################################

        SOCRATES.allocate_aer_prsc(atmos.aer, atmos.dimen, atmos.spectrum)
        if atmos.control.l_aerosol
            error("Aerosols not implemented")
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
            atmos.control.i_cloud_representation = SOCRATES.rad_pcf.ip_cloud_homogen # Ice and water mixed homogeneously (-K 1) = ip_cloud_homogen ; Cloud mixing liquid and ice (-K 2) = ip_cloud_ice_water
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

        atmos.mask_l =            zeros(Int, atmos.nlev_c)
        atmos.mask_c =            zeros(Int, atmos.nlev_c)
        atmos.flux_l =            zeros(Float64, atmos.nlev_l)  # Latent heat
        atmos.flux_cdry =         zeros(Float64, atmos.nlev_l)  # Dry convection 
        atmos.flux_cdct =         zeros(Float64, atmos.nlev_l)  # Conduction 
        atmos.Kzz =               zeros(Float64, atmos.nlev_l)

        atmos.flux_tot =          zeros(Float64, atmos.nlev_l)
        atmos.flux_dif =          zeros(Float64, atmos.nlev_c)
        atmos.ediv_add =          zeros(Float64, atmos.nlev_c)
        atmos.heating_rate =      zeros(Float64, atmos.nlev_c)

        # Mark as allocated
        atmos.is_alloc = true
        @debug "Allocate complete"

        return nothing
    end  # end of allocate


    """
    **Calculate composition assuming chemical equilibrium at each level.**

    Uses FastChem to calculate the gas composition at each level of the atmosphere. Volatiles are converted to bulk
    elemental abundances, which are then provided to FastChem alongside the temperature/pressure profile. FastChem is
    currently called as an executable, which is not optimal.

    This function DOES NOT automatically recalculate layer properties (e.g. mmw, density, height).

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `chem_type::Int`                  chemistry type (see wiki)
    - `write_cfg::Bool`                 write config and elements
    - `tmp_floor::Float64=500.0`        temperature floor for T(p) provided to FastChem

    Returns:
    - `state::Int`                      fastchem state (0: success, 1: critical_fail, 2: elem_fail, 3: conv_fail, 4: both_fail)
    """
    function chemistry_eq!(atmos::atmosphere.Atmos_t, chem_type::Int, write_cfg::Bool; tmp_floor::Float64=500.0)::Int

        # Return code 
        state::Int = 0

        # Check fastchem enabled 
        if !atmos.fastchem_flag
            @warn "Fastchem is not enabled but `chemistry_eq!` was called"
            return 1
        end 

        # Paths
        execpath::String = joinpath(atmos.fastchem_path,"fastchem")             # Executable file 
        confpath::String = joinpath(atmos.fastchem_work,"config.input")         # Configuration by AGNI
        chempath::String = joinpath(atmos.fastchem_work,"chemistry.dat")        # Chemistry by FastChem

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
                logK = joinpath(atmos.fastchem_path,"input/","logK/")
                write(f,joinpath(logK,"logK.dat")*" "*joinpath(logK,"logK_condensates.dat")*" \n\n")

                write(f,"#Accuracy of chemistry iteration \n")
                write(f,"1.0e-5 \n\n")
                
                write(f,"#Accuracy of element conservation \n")
                write(f,"1.0e-4 \n\n")

                write(f,"#Max number of chemistry iterations  \n")
                write(f,"50000 \n\n")

                write(f,"#Max number internal solver iterations  \n")
                write(f,"10000 \n\n")
            end

            # Calculate elemental abundances 
            # number densities normalised relative to hydrogen 
            # for each element X, value = log10(N_X/N_H) + 12 
            # N = X(P/(K*T) , where X is the VMR and K is boltz-const
            N_t = zeros(Float64, length(phys.elements_list))                # total
            N_g = zeros(Float64, length(phys.elements_list))                # this gas
            for gas in atmos.gas_names
                d = phys.count_atoms(gas)
                fill!(N_g, 0.0)
                for (i,e) in enumerate(phys.elements_list)
                    if e in keys(d)
                        N_g[i] += d[e]
                    end 
                end 
                N_g *= get_x(atmos, gas, atmos.nlev_c) * atmos.p[end] / (phys.k_B * atmos.tmp[end])  # gas contribution
                N_t += N_g  # number/m^3
            end 

            # Write elemental abundances 
            open(joinpath(atmos.fastchem_work,"elements.dat"),"w") do f
                write(f,"# Elemental abundances derived from AGNI volatiles \n")
                for (i,e) in enumerate(phys.elements_list)
                    if N_t[i] > 1.0e-30
                        write(f, @sprintf("%s    %.3f \n",e,log10(N_t[i]/N_t[1]) + 12.0))
                    end
                end 
            end 
        end

        # Write PT profile 
        open(joinpath(atmos.fastchem_work,"pt.dat"),"w") do f
            write(f,"# AGNI temperature structure \n")
            write(f,"# bar, kelvin \n")
            for i in 1:atmos.nlev_c 
                write(f,@sprintf("%.6e    %.6e \n", atmos.p[i], max(tmp_floor,atmos.tmp[i]) ))
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

        # Parse gas chemistry
        g::String = ""
        N_t = data[4,:] # at each level: sum of gas number densities 
        for (i,h) in enumerate(head)  # for each column (gas)
            g = rstrip(lstrip(h))

            # check if gas is included in the model
            if g in keys(phys.map_fastchem_name) 
                g = phys.map_fastchem_name[g]  # convert name
                if g in atmos.gas_names
                    N_g = data[i,:]  # number densities for this gas 
                    atmos.gas_vmr[g][:] .= N_g[:] ./ N_t[:]    # mole fraction (VMR) for this gas 
                end
            end # not included => skip
        end 

        # Do not renormalise mixing ratios

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
            if atmos.gas_phase[c][i] 
                x_con += atmos.gas_vmr[c][i]
            end 
        end 
            
        # Calculate current and target dry VMRs
        x_dry = 1.0 - x_con 
        x_dry_old = 0.0
        for g in atmos.gas_names 
            # skip condensing gases, since their VMR is set by saturation
            if !atmos.gas_phase[g][i]
                x_dry_old += atmos.gas_vmr[g][i]
            end 
        end 

        # Renormalise VMR to unity, scaling DRY COMPONENTS ONLY
        for g in atmos.gas_names 
            if !atmos.gas_phase[g][i]
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
    until the gas is exactly saturated, not supersaturated. Rain is then 
    re-evaporated at lower levels without supersaturating them.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    """
    function handle_saturation!(atmos::atmosphere.Atmos_t)

        # Parameters 
        evap_enabled::Bool =            true    # Enable re-vaporation of rain
        evap_efficiency::Float64 =      1.0     # Evaporation efficiency
        cond_retention_frac::Float64 =  0.5     # Condensate retention fraction (0 => complete rainout)

        # Work arrays 
        maxvmr::Dict{String, Float64} = Dict{String, Float64}()     # max running VMR for each condensable 
        cond_kg::Dict{String,Float64} = Dict{String, Float64}()     # condensed kg for each condensable 
        rain_kg::Dict{String,Float64} = Dict{String, Float64}()     # rainout kg for each condensable
        x_sat::Float64 =        0.0
        layer_area::Float64 =   0.0
        supcrit::Bool =         false

        # Set maximum value (for cold trapping)
        for c in atmos.condensates 
            maxvmr[c] = atmos.gas_vmr[c][end]
        end 

        # Reset mixing ratios to surface values
        # Reset phase change flags 
        # Reset condensation yield values 
        for g in atmos.gas_names
           atmos.gas_vmr[g][1:end-1] .= atmos.gas_vmr[g][end]
           atmos.gas_phase[g][:] .= false
           atmos.gas_cprod[g][:] .= 0.0
        end 

        # Reset water cloud
        fill!(atmos.cloud_arr_r, 0.0)
        fill!(atmos.cloud_arr_l, 0.0)
        fill!(atmos.cloud_arr_f, 0.0)


        # Loop from bottom to top (does not include bottommost level)
        for i in range(start=atmos.nlev_c-1, stop=1, step=-1)

            layer_area = 4.0*pi*(atmos.z[i]+atmos.rp)^2.0

            # For each condensate 
            for c in atmos.condensates

                # Reset condensation and rain
                cond_kg[c] = 0.0                    # kg of 'c' condensate produced at this level 
                rain_kg[c] = 0.0                    # kg of ^ rained-out from this level

                # check criticality 
                supcrit = atmos.tmp[i] > atmos.gas_dat[c].T_crit+1.0e-10

                # saturation mixing ratio 
                x_sat = phys.get_Psat(atmos.gas_dat[c], atmos.tmp[i]) / atmos.p[i]

                # cold trap 
                if atmos.gas_vmr[c][i] > maxvmr[c]
                    atmos.gas_vmr[c][i] = maxvmr[c]
                    atmos.gas_phase[c][i] = true 
                end 
                
                # condense if supersaturated
                if (atmos.gas_vmr[c][i] > x_sat) && !supcrit

                    # set rainout kg 
                    cond_kg[c] = layer_area * (atmos.gas_dat[c].mmw/atmos.layer_mmw[i])*atmos.p[i]*(atmos.gas_vmr[c][i] - x_sat)/atmos.layer_grav[i]
                    rain_kg[c] = cond_kg[c] * (1.0 - cond_retention_frac)

                    # condensation yield at this level 
                    atmos.gas_cprod[c][i] += cond_kg[c]

                    # set new vmr
                    atmos.gas_vmr[c][i] = x_sat

                    # store vmr for cold trapping at levels above this one
                    maxvmr[c] = x_sat

                    # flag condensate as actively condensing at this level
                    atmos.gas_phase[c][i] = true 
                end 
            end # end condensates

            # for c in atmos.condensates
            #     @printf("At %d: generated %s rain: %.3e kg \n", i, c, rain_kg[c])
            # end

            normalise_vmrs!(atmos, i)

            # Recalculate layer mmw 
            atmos.layer_mmw[i] = 0.0
            for g in atmos.gas_names
                atmos.layer_mmw[i] += atmos.gas_vmr[g][i] * atmos.gas_dat[g].mmw
            end

            # Set water cloud at this level
            if "H2O" in atmos.condensates
                # mass mixing ratio (take ratio of mass surface densities [kg/m^2])
                atmos.cloud_arr_l[i] = (cond_kg["H2O"] * cond_retention_frac / layer_area) / atmos.layer_mass[i]

                if atmos.cloud_arr_l[i] > 1.0
                    @warn "Water cloud mass mixing ratio is greater than unity (level $i)"
                end 

                # droplet radius and area fraction (fixed values)
                atmos.cloud_arr_r[i] = atmos.cloud_val_r
                atmos.cloud_arr_f[i] = atmos.cloud_val_f
            end 

            # Evaporate condensate withon unsaturated layers below this one
            #   integrate downwards from this condensing layer (i)
            if evap_enabled
                for c in atmos.condensates  # loop over condensates 

                    # no rain => skip this condensate 
                    if rain_kg[c] < 1.0e-10
                        continue 
                    end 

                    # loop downwards from layer i to almost surface
                    for j in range(start=i+1, stop=atmos.nlev_c-1, step=1)
                        # wait for bottom of condensing region  
                        if atmos.gas_phase[c][j]
                            continue 
                        end 

                        # layer area
                        layer_area = 4.0*pi*(atmos.z[j]+atmos.rp)^2.0

                        # Calculate partial pressure change required for saturation
                        supcrit = atmos.tmp[j] > atmos.gas_dat[c].T_crit
                        if !supcrit
                            # subcritical (can evaporate up to saturation)
                            dp_sat = phys.get_Psat(atmos.gas_dat[c], atmos.tmp[j]) - atmos.gas_vmr[c][j]*atmos.p[j]
                        else 
                            # supercritical (can take up any amount of evaporated rain - just choose a big number)
                            dp_sat = 1.0e99
                            rain_kg[c] = 0.0
                        end 

                        # Calculate kg of gas that would saturate 
                        dm_sat = layer_area*(atmos.gas_dat[c].mmw/atmos.layer_mmw[j])*dp_sat/atmos.layer_grav[j]
                        
                        # can we evaporate all rain within this layer?
                        if rain_kg[c] < dm_sat 
                            dm_sat = rain_kg[c]
                        end 

                        # Evaporation efficiency factor 
                        #   This fraction of the rain that *could* be evaporated
                        #   at this layer *is* converted to vapour in this layer.
                        dm_sat *= evap_efficiency

                        # offset yield at this level by evaporation 
                        atmos.gas_cprod[c][j] -= dm_sat

                        # additional partial pressure from evaporation 
                        dp_sat = dm_sat * atmos.layer_grav[j] / (layer_area * atmos.layer_mmw[j] / atmos.gas_dat[c].mmw)

                        # convert extra pp to extra vmr
                        atmos.gas_vmr[c][j] += dp_sat / atmos.p[j]
                        
                        # reduce total kg of rain correspondingly 
                        rain_kg[c] -= dm_sat

                        # flag this region as containing a phase change 
                        if !supcrit
                            atmos.gas_phase[c][j] = true 
                            # @printf("%s at %d: evaporating %.2e kg \n", c, j, dm_sat)
                        end

                        normalise_vmrs!(atmos, j)

                        # Recalculate layer mmw 
                        atmos.layer_mmw[j] = 0.0
                        for g in atmos.gas_names
                            atmos.layer_mmw[j] += atmos.gas_vmr[g][j] * atmos.gas_dat[g].mmw
                        end

                        # all rain evaporated?
                        if rain_kg[c] < 1.0e-10 
                            break 
                        end 

                    end # go to next j level (below)
                    
                    # @printf("      Surface: Landed %s rain %.3e kg \n", c, rain_kg[c])

                end # end condensates
            end # end evaporation 

            # Calculate layer properties 
            # calc_layer_props!(atmos)

        end # go to next i level (above)

        return nothing 
    end 

    # Set cloud properties within condensing regions
    function water_cloud!(atmos::atmosphere.Atmos_t)

        # Reset 
        fill!(atmos.cloud_arr_r, 0.0)
        fill!(atmos.cloud_arr_l, 0.0)
        fill!(atmos.cloud_arr_f, 0.0)

        # Get index of water 
        if !("H2O" in atmos.gas_names)
            return nothing
        end

        # Check if each level is condensing. If it is, place on phase curve
        x::Float64      = 0.0  # vmr of water vapour
        Tsat::Float64   = 0.0  # dew point temperature of water vapour
        for i in 1:atmos.nlev_c

            x = atmos.gas_vmr["H2O"][i]
            if x < 1.0e-10 
                continue
            end

            # Cell centre only
            # Check if partial pressure > P_sat
            if atmos.p[i]*atmos.gas_vmr["H2O"][i] > phys.get_Psat(atmos.gas_dat["H2O"], atmos.tmp[i])
                atmos.cloud_arr_r[i] = atmos.cloud_val_r
                atmos.cloud_arr_l[i] = atmos.cloud_val_l
                atmos.cloud_arr_f[i] = atmos.cloud_val_f
            end
        end
        
        return nothing
    end

    # Smooth temperature at cell-centres 
    function smooth_centres!(atmos::atmosphere.Atmos_t, width::Int)

        if width > 2
            if mod(width,2) == 0
                width += 1 # window width has to be an odd number
            end
            atmos.tmp = moving_average.hma(atmos.tmp, width)
        end 
        clamp!(atmos.tmp, atmos.tmp_floor, atmos.tmp_ceiling)

        return nothing
    end

    """
    **Set cell-edge temperatures from cell-centre values.**

    Uses interpolation within the bulk of the column and extrapolation for the 
    topmost edge. 

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `back_interp::Bool=false`         interpolate resultant tmpl back to tmp (should be avoided)
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

    # Get interleaved cell-centre and cell-edge PT grid
    function get_interleaved_pt(atmos::atmosphere.Atmos_t)
        arr_T::Array{Float64, 1} = zeros(Float64, atmos.nlev_c + atmos.nlev_l)
        arr_P::Array{Float64, 1} = zeros(Float64, atmos.nlev_c + atmos.nlev_l)
        idx::Int = 1

        # top
        arr_T[1] = atmos.tmpl[1]
        arr_P[1] = atmos.pl[1]

        # middle 
        for i in 1:atmos.nlev_c
            idx = (i-1)*2
            arr_T[idx+1] = atmos.tmpl[i]
            arr_T[idx+2] = atmos.tmp[i]
            arr_P[idx+1] = atmos.pl[i]
            arr_P[idx+2] = atmos.p[i]
        end 

        # bottom
        arr_T[end] = atmos.tmpl[end]
        arr_P[end] = atmos.pl[end]
        
        return arr_P, arr_T
    end 

end
