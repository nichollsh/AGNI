# Contains the atmosphere module

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

"""
**Main atmosphere module for storing model data**

This module primarily contains the `Atmos_t` struct which does most of the heavy lifting
in AGNI. There are also functions for setting up and configuring the struct.

Also includes hydrostatic integrator and ways to estimate observable quantities.
"""
module atmosphere

    # System libraries
    using Printf
    using LinearAlgebra
    using Logging
    import Statistics
    import Interpolations: interpolate, Gridded, Linear, Flat, Line, extrapolate, Extrapolation
    import DelimitedFiles:readdlm

    # SOCRATES library
    const SOCRATESjl = abspath(ENV["RAD_DIR"], "julia","src","SOCRATES.jl")
    include(SOCRATESjl)

    # Local modules
    import ..phys
    import ..spectrum

    # Code versions
    const AGNI_VERSION::String     = "1.8.6"  # current agni version
    const SOCVER_minimum::Float64  = 2407.2    # minimum required socrates version

    # Hydrostatic+gravity+mass calculation (constants and limits)
    HYDROGRAV_steps::Int64   = 9000      # total number of steps in height integration
    HYDROGRAV_maxdr::Float64 = 1e9       # maximum dz across each layer [m]
    HYDROGRAV_mindr::Float64 = 1e-5      # minimum dz across each layer [m]
    HYDROGRAV_ming::Float64  = 1e-10     # minimum allowed gravity [m/s^2]
    HYDROGRAV_constg::Bool   = false     # constant gravity with height?
    HYDROGRAV_selfg::Bool    = true      # include self-gravity of the atmosphere?

    # Other constants
    const UNSET_STR::String             = "__AGNI_UNSET_STR"
    const NLEV_minimum::Int             = 25        # minimum allowed number of levels
    const PHS_TIMESCALE_MIN::Float64    = 0.01      # minimum phase change timescale [s]
    const SURF_ROUGHNESS_MIN::Float64   = 1e-5      # [m]
    const SURF_WINDSPEED_MIN::Float64   = 1e-5      # [m/s]
    const RP_MIN::Float64               = 1e4       # [m]
    const SKIN_D_MIN::Float64           = 1e-6      # [m]
    const SKIN_K_MIN::Float64           = 1e-6      # [W K-1 m-1]
    const COND_DISALLOWED::Array        = ["H2","He"]

    # Pressure grid
    const PRESSURE_RATIO_MIN::Float64   = 1.0001    # minimum p_boa/p_toa ratio
    const PRESSURE_FACT_BOT::Float64    = 0.6       # Pressure factor at bottom layer

    # Enum of available radiative transfer schemes
    @enum RTSCHEME RT_SOCRATES=1 RT_GREYGAS=2

    """
    **Parameters for deep atmospheric heating.**

    AGNI treats deep heating as an *additional upward energy flux* that is deposited
    within the column following a Gaussian profile in log-pressure space.

    Two power modes are supported:
    - `:efficiency` — heating flux = `efficiency × instellation` (stellar efficiency)
    - `:flux`       — heating flux = `F_total` (fixed radiative flux in W m⁻²)

    Additional design choices:
    - Pressure-normalised deposition (dF/dP profile; legacy behaviour)
    - Mass-normalised deposition (dm-weighted; aligns with ∂L/∂m forms)
    - Handling deposition below the model bottom as a *bottom boundary flux*

    Notes:
    - Pressures are in Pa throughout AGNI.
    - `efficiency` is always clamped to [0,1].

    Fields:
    - `active::Bool`              Enable/disable deep heating.
    - `P_dep::Float64`            Deposition pressure centre [Pa].
    - `sigma_P::Float64`          Width of Gaussian in log-pressure space [dimensionless].
    - `efficiency::Float64`       Heating efficiency (fraction of instellation), used when `power_mode=:efficiency`.
    - `normalization::Symbol`     `:pressure` (legacy) or `:mass` (dm-weighted).
    - `below_domain::Symbol`      `:clamp` or `:boundary_flux`.
    - `power_mode::Symbol`        `:efficiency` | `:flux`.
    - `F_total::Float64`          Total deposited flux [W m-2], used when `power_mode=:flux`.
    """
    struct DeepHeatingParams
        active::Bool
        P_dep::Float64
        sigma_P::Float64
        efficiency::Float64
        normalization::Symbol
        below_domain::Symbol
        power_mode::Symbol
        F_total::Float64
    end

    # Default constructor with heating disabled
    DeepHeatingParams() = DeepHeatingParams(false, 1.0e5, 1.0, 0.0,
                                            :pressure, :clamp,
                                            :efficiency, 0.0)

    # Contains data pertaining to the atmosphere (fluxes, temperature, etc.)
    mutable struct Atmos_t

        # AGNI version used to generate this struct
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
        FC_DIR::String          # path to fastchem install folder
        RFM_DIR::String         # path to RFM install folder
        FRAMES_DIR::String      # path to frames of animation
        IO_DIR::String          # path to temporary directory, for fast I/O

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

        # Radiation scheme
        rt_scheme::RTSCHEME             # RT scheme (1: SOCRATES, 2: Grey gas)
        benchmark::Bool                 # Benchmark RT?
        num_rt_eval::Int                # Total number of RT evaluations
        tim_rt_eval::UInt64             # Total time spent doing RT evaluations [ns]

        # Radiation parameters
        all_channels::Bool              # Use all bands?
        spectral_file::String           # Path to spectral file
        star_file::String               # Path to star spectrum
        albedo_b::Float64               # Enforced bond albedo
        zenith_degrees::Float64         # Solar zenith angle [deg]
        toa_heating::Float64            # Derived downward shortwave radiation flux at topmost level [W m-2]
        instellation::Float64           # Solar flux at top of atmopshere [W m-2]
        s0_fact::Float64                # Scale factor to instellation (see Cronin+14)
        overlap_method::String          # Absorber overlap method to be used
        κ_grey_lw::Float64              # LW opacity used for grey-gas scheme [m2 kg-1]
        κ_grey_sw::Float64              # SW opacity used for grey-gas scheme [m2 kg-1]

        # Spectral bands
        nbands::Int
        bands_min::Array{Float64,1}    # Lower wavelength [m]
        bands_max::Array{Float64,1}    # Upper wavelength [m]
        bands_cen::Array{Float64,1}    # Midpoint [m]
        bands_wid::Array{Float64,1}    # Width [m]

        # Transparent model?
        transparent::Bool             # Atmosphere configured to be transparent?

        # Pressure-temperature grid (with i=1 at the top of the model)
        nlev_c::Int             # Cell centre (count)
        nlev_l::Int             # Cell edge (count)
        p_oboa::Float64         # Pressure at bottom [Pa], original
        p_boa::Float64          # Pressure at bottom [Pa], calculated
        p_toa::Float64          # Pressure at top [Pa]
        tmp::Array{Float64,1}   # cc temperature [K]
        tmpl::Array{Float64,1}  # ce temperature [K]
        p::Array{Float64,1}     # cc pressure [Pa]
        pl::Array{Float64,1}    # ce pressure [Pa]
        r::Array{Float64,1}     # cc radius [m]
        rl::Array{Float64,1}    # ce radius [m]
        g::Array{Float64,1}     # cc gravity [m s-2]
        gl::Array{Float64,1}    # ce gravity [m s-2]
        m::Array{Float64,1}     # cc mass encl [kg]
        ml::Array{Float64,1}    # ce mass encl [kg]

        tmp_floor::Float64      # Temperature floor to prevent numerics [K]
        tmp_ceiling::Float64    # Temperature ceiling to prevent numerics [K]

        # Surface
        surface_material::String        # Surface material file path
        albedo_s::Float64               # Grey surface albedo when surface=greybody
        surf_r_arr::Array{Float64,1}    # Spectral surface spherical reflectance
        surf_e_arr::Array{Float64,1}    # Spectral surface hemispheric emissivity
        tmp_surf::Float64               # Surface brightness temperature [K]
        rp::Float64                     # surface radius [m]
        grav_surf::Float64              # surface gravity [m s-2]
        skin_d::Float64                 # skin thickness [m]
        skin_k::Float64                 # skin thermal conductivity [W m-1 K-1] (You can find reasonable values here: https://doi.org/10.1016/S1474-7065(03)00069-X)
        tmp_magma::Float64              # Mantle temperature [K]

        # Gas tracking variables (incl gases which are not in spectralfile)
        gas_num::Int                                # Number of gases
        gas_names::Array{String,1}                  # List of gas names
        gas_dat::Dict{String, phys.Gas_t}           # Struct variables containing thermodynamic data for each gas

        # Chemistry and composition
        gas_vmr::Dict{String, Array{Float64,1}}     # runtime calculated VMRs in dict, (key,value) = (gas_name,array)
        gas_cvmr::Dict{String, Array{Float64,1}}    # runtime calculated VMR values after chemistry calculation, before rainout
        gas_ovmr::Dict{String, Array{Float64,1}}    # original VMR values at model initialisation
        metal_orig::Dict{String, Float64}           # input elem ratios (elem num density rel to hydrogen)
        metal_calc::Dict{String, Float64}           # calc'd elem ratios, from gas mixing ratios at surface

        # Whether each gas is reliably 'safe', in terms of the physics modelled
        gas_safe::Dict{String, Bool}

        # Condensation and evaporation aloft
        gas_sat::Dict{String, Array{Bool, 1}}       # Gas is saturated or cold-trapped in each layer?
        cond_yield::Dict{String, Array{Float64,1}}  # condensate yield [kg/m^2] at each level (can be negative, representing evaporation)
        cond_accum::Dict{String, Float64}           # condensate accumulation left after evaporation aloft (implicit surface liquid) [kg/m^2]
        condensates::Array{String, 1}               # List of condensing gases (strings)
        condense_any::Bool                          # length(condensates)>0 ?

        # Ocean tracking variables
        ocean_ini::Dict{String, Float64}    # INPUT: ocean reservoir from user [kg/m^2] - does not change
        ocean_tot::Dict{String, Float64}    # OUTPUT: ocean reservoir [kg/m^2] - after cond/evap, both from surf and aloft

        # Ocean layering variables (surface liquid layering)
        ocean_ob_frac::Float64              # INPUT: ocean basin area, as fraction of planet surface
        ocean_cs_height::Float64            # INPUT: continental shelf height [m]
        ocean_layers::Array{Tuple,1}        # OUTPUT: layer structure of surface liquids
        ocean_maxdepth::Float64             # OUTPUT: ocean depth at deepest part [m]
        ocean_areacov::Float64              # OUTPUT: fraction of planet surface covered by oceans
        ocean_topliq::String                # OUTPUT: name of top-most ocean component

        # Gases (only those in SOCRATES spectralfile)
        gas_soc_num::Int                    # number of gases
        gas_soc_names::Array{String,1}      # names of each gas (as a list)

        # Layers' average properties
        real_gas::Bool                      # use real-gas equations of state where possible
        thermo_funct::Bool                  # use temperature-dependent evaluation of thermodynamic properties
        layer_ρ::Array{Float64,1}           # mass density [kg m-3]
        layer_σ::Array{Float64,1}           # mass of each layer, per unit area  [kg m-2]
        layer_μ::Array{Float64,1}           # mean molecular weight [kg mol-1]
        layer_cp::Array{Float64,1}          # heat capacity at const-p [J K-1 kg-1]
        layer_kc::Array{Float64,1}          # thermal conductivity at const-p [W m-1 K-1]
        layer_thick::Array{Float64,1}       # geometrical thickness [m]
        layer_isbound::Array{Bool,1}        # is this layer strongly bound by gravity?

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

        # Surface planck emission (incl. emissivity)
        surf_flux::Array{Float64, 1}

        # Contribution function (to outgoing flux) per-band
        contfunc_band::Array{Float64,2}     # LW only, not normalised

        # RFM line-by-line calculation
        rfm_fl::Array{Float64,1}            # upward flux [erg/(s cm2 cm-1)]
        rfm_wn::Array{Float64,1}            # wavenumber array [cm-1]
        rfm_npts::Int                       # number of points

        # Sensible heating
        C_d::Float64                        # Turbulent exchange coefficient, to be calc'd
        surf_windspeed::Float64             # Surface wind speed [m s-1]
        surf_roughness::Float64             # Surface roughness scale [m]
        flux_sens::Float64                  # Turbulent flux [W m-2]


        # Convection
        mlt_asymptotic::Bool                # INPUT: Mixing length scales asymptotically, but ~0 near ground
        mlt_criterion::Char                 # INPUT: Stability criterion. Options: (s)chwarzschild, (l)edoux
        Kzz_floor::Float64                  # INPUT: Kzz floor [m2 s-1]
        Kzz_ceiling::Float64                # INPUT: Kzz ceiling [m2 s-1]
        Kzz_pbreak::Float64                 # INPUT: Kzz break point pressure [Pa]
        Kzz_kbreak::Float64                 # INPUT: Kzz break point diffusion [m2 s-1]
        mask_c::Array{Bool,1}               # OUT: Layers transporting convective flux
        flux_cdry::Array{Float64,1}         # OUT: Dry convective fluxes from MLT
        Kzz::Array{Float64,1}               # OUT: Eddy diffusion coefficient from MLT
        w_conv::Array{Float64,1}            # OUT: Convective velocities [m s-1]
        λ_conv::Array{Float64,1}            # OUT: Mixing lengths [m s-1]

        # Conduction
        flux_cdct::Array{Float64,1}         # Conductive flux [W m-2]

        # Phase change
        phs_timescale::Float64                # Time scale (mixed composition)
        evap_efficiency::Float64            # Base evaporation efficiency of rain (0 to 1)
        phs_wrk_df::Array{Float64,1}        # work array: flux difference
        phs_wrk_fl::Array{Float64,1}        # work array: edge fluxes
        flux_l::Array{Float64, 1}           # Latent heat energy flux [W m-2]
        mask_l::Array{Bool,1}               # Layers transporting latent heat

        # Clouds
        cloud_alpha::Float64                # Condensate cloud formation fraction
        cloud_arr_r::Array{Float64,1}       # Characteristic dimensions of condensed species [m].
        cloud_arr_l::Array{Float64,1}       # Mass mixing ratios of condensate. 0 : saturated water vapor does not turn liquid ; 1 : the entire mass of the cell contributes to the cloud
        cloud_arr_f::Array{Float64,1}       # Total cloud area fraction in layers. 0 : clear sky cell ; 1 : the cloud takes over the entire area of the Cell
        cloud_val_r::Float64                # \
        cloud_val_l::Float64                #  |-> Default scalar values to above arrays
        cloud_val_f::Float64                # /

        # Cell-internal heating
        flux_advect::Array{Float64, 1}     # Energy flux advected into each cell [W m-2], treated as entering into its bottom edge

        # Deep atmospheric heating
        deep_heating::DeepHeatingParams     # Parameters for deep heating
        flux_deep::Array{Float64,1}         # Deep heating flux at cell edges [W m-2] # should add at lw flux level

        # Total energy flux
        flux_dif::Array{Float64,1}      # Flux lost at each level [W m-2] (positive is heating up)
        flux_tot::Array{Float64,1}      # Total upward-directed flux at cell edges [W m-2]

        # Heating rate felt at each level [K/day]
        heating_rate::Array{Float64,1}

        # Diagnostic profiles
        timescale_conv::Array{Float64,1}    # Convective timescale [s]
        timescale_rad::Array{Float64,1}     # Radiative timescale [s]
        diagnostic_Ra::Array{Float64,1}     # Rayleigh number at each layer

        # FastChem equilibrium chemistry
        flag_fastchem::Bool             # Fastchem enabled?
        fastchem_floor::Float64         # Minimum temperature allowed to be sent to FC
        fastchem_maxiter_chem::Int      # Maximum FC iterations (chemistry)
        fastchem_maxiter_solv::Int      # Maximum FC iterations (internal solver)
        fastchem_xtol_chem::Float64     # FC solver tolerance (chemistry)
        fastchem_xtol_elem::Float64     # FC solver tolerance (elemental)
        fastchem_exec::String           # Path to executable
        fastchem_work::String           # Path to working directory
        fastchem_conf::String           # Path to input config file
        fastchem_elem::String           # Path to input metallicities file
        fastchem_prof::String           # Path to input T(p) profile
        fastchem_chem::String           # Path to output gas composition
        fastchem_cond::String           # Path to output condensate phases
        fastchem_moni::String           # Path to output monitor file
        fastchem_wellmixed::Bool        # True: perform chemistry at Tsurf. False: do 1D profile.

        # RFM radiative transfer
        flag_rfm::Bool                  # RFM enabled?
        rfm_exec::String                # Path to rfm executable
        rfm_work::String                # Path to rfm working directory
        rfm_parfile::String             # Path to rfm parfile. If empty, do not run RFM.

        # Observing properties
        transspec_p::Float64            # Pressure level probed in transmission [Pa]
        transspec_r::Float64            # planet radius probed in transmission [m]
        transspec_μ::Float64            # mmw probed in transmission [kg mol-1]
        transspec_tmp::Float64          # temperature probed in transmission [K]
        transspec_grav::Float64         # gravity at layer probed in transmission [m s-2]
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

    # Check parameter is within valid range
    function _check_range(name, val; min=nothing, max=nothing)::Bool
        if !isnothing(min) && !isnothing(max) && ((val<min) || (val>max))
            @error "$name is out of range"
            @error "    Got: $min < $val < $max"
        elseif !isnothing(min) && (val < min)
            @error "$name is too small"
            @error "    Got: $val < $min"
        elseif !isnothing(max) && (val > max)
            @error "$name is too large"
            @error "    Got: $val > $max"
        else
            return true
        end
        return false
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
    - `spfile::String`                  path to spectral file ('grey' => use grey gas)
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
    - `IO_DIR::String`                  directory used for fast file operations.
    - `condensates`                     list of condensates (gas names).
    - `metallicities::Dict`             dictionary of elemental metallicities (mass ratio rel to hydrogen)
    - `surface_material::String`        surface material (default is "greybody", but can point to file instead).
    - `albedo_s::Float64`               grey surface albedo used when `surface_material="greybody"`.
    - `tmp_floor::Float64`              temperature floor [K].
    - `surf_roughness::Float64`         surface roughness length scale [m]
    - `surf_windspeed::Float64`         surface wind speed [m s-1].
    - `Kzz_floor::Float64`              min eddy diffusion coefficient, cgs units [cm2 s-1]
    - `mlt_asymptotic::Bool`            mixing length scales asymptotically, but ~0 near ground
    - `mlt_criterion::Char`             MLT stability criterion. Options: (s)chwarzschild, (l)edoux.
    - `tmp_magma::Float64`              mantle temperature [K] for sol_type==2.
    - `skin_d::Float64`                 skin thickness [m].
    - `skin_k::Float64`                 skin thermal conductivity [W m-1 K-1].
    - `overlap_method::String`          gaseous overlap scheme (ro: rand overlap, ee: equiv extinct, rorr: ro+resort+rebin).
    - `target_olr::Float64`             target OLR [W m-2] for sol_type==4.
    - `flux_int::Float64`               planet's internal flux for sol_type==3.
    - `all_channels::Bool`              use all channels available for RT?
    - `flag_rayleigh::Bool`             include rayleigh scattering?
    - `flag_gcontinuum::Bool`           include generalised continuum absorption?
    - `flag_continuum::Bool`            include continuum absorption?
    - `flag_aerosol::Bool`              include aersols?
    - `flag_cloud::Bool`                include clouds?
    - `phs_timescale::Float64`          phase change timescale [s]
    - `evap_efficiency::Float64`        re-evaporatione efficiency compared to saturating amount
    - `real_gas::Bool`                  use real gas EOS where possible
    - `thermo_functions::Bool`          use temperature-dependent thermodynamic properties
    - `use_all_gases::Bool`             store information on all supported gases, incl those not provided in cfg
    - `check_integrity::Bool`           confirm integrity of thermo files using their checksum
    - `κ_grey_lw::Float64`              gas opacity when using grey-gas RT scheme, longwave
    - `κ_grey_sw::Float64`              gas opacity when using grey-gas RT scheme, shortwave
    - `fastchem_work::String`           working directory for fastchem
    - `fastchem_floor::Float64`         temperature floor on profile provided to fastchem
    - `fastchem_maxiter_chem::Float64`  maximum chemical iterations allowed by fastchem
    - `fastchem_maxiter_solv:Float64`   maximum solver iterations allowed by fastchem
    - `fastchem_xtol_chem::Float64`     solution tolerance required of fastchem (chemical)
    - `fastchem_xtol_elem::Float64`     solution tolerance required of fastchem (elemental)
    - `rfm_parfile::String`             path to HITRAN-formatted .par file provided to RFM

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

                    IO_DIR::String   =          UNSET_STR,
                    condensates =               String[],
                    metallicities::Dict =       Dict{String,Float64}(),
                    surface_material::String =  "greybody",
                    albedo_s::Float64 =         0.0,
                    tmp_floor::Float64 =        2.0,
                    tmp_ceiling::Float64 =      2e4,
                    surf_roughness::Float64 =   0.001,
                    surf_windspeed::Float64 =   2.0,
                    Kzz_floor::Float64 =        1e5,
                    mlt_asymptotic::Bool =      true,
                    mlt_criterion::Char =       's',
                    tmp_magma::Float64 =        3000.0,
                    skin_d::Float64 =           0.05,
                    skin_k::Float64 =           2.0,
                    overlap_method::String =    "ee",
                    target_olr::Float64 =       0.0,
                    flux_int::Float64 =         0.0,
                    all_channels::Bool  =       true,
                    flag_rayleigh::Bool =       false,
                    flag_gcontinuum::Bool =     false,
                    flag_continuum::Bool =      false,
                    flag_aerosol::Bool =        false,
                    flag_cloud::Bool =          false,

                    phs_timescale::Float64 =    1e6,
                    evap_efficiency::Float64 =  0.05,

                    real_gas::Bool =            true,
                    thermo_functions::Bool =    true,
                    use_all_gases::Bool =       false,
                    check_integrity::Bool =     true,

                    κ_grey_lw::Float64  =       8e-4,
                    κ_grey_sw::Float64  =       2e-4,

                    fastchem_work::String       =  UNSET_STR,
                    fastchem_floor::Float64     =  400.0,
                    fastchem_maxiter_chem::Int  =  80000,
                    fastchem_maxiter_solv::Int  =  40000,
                    fastchem_xtol_chem::Float64 =  1.0e-3,
                    fastchem_xtol_elem::Float64 =  1.0e-3,
                    fastchem_wellmixed::Bool    =  false,

                    rfm_parfile::String =       UNSET_STR,

                    ocean_ob_frac::Float64 =    0.6,
                    ocean_cs_height::Float64 =  3000.0
                    )::Bool

        # Say hello
        @info  "Setting-up a new atmosphere struct"
        atmos.AGNI_VERSION = AGNI_VERSION
        @debug "AGNI VERSION = "*AGNI_VERSION

        # -------------------------
        # Directories
        # -------------------------

        # Set AGNI root directory
        atmos.ROOT_DIR = abspath(ROOT_DIR)
        if !isfile(joinpath(atmos.ROOT_DIR,"agni.jl"))
            @error "Cannot find `agni.jl` in the root directory provided"
            @error "     ROOT_DIR=$(atmos.ROOT_DIR)"
            return false
        end

        # Locate thermodynamics dir
        atmos.THERMO_DIR = joinpath(atmos.ROOT_DIR, "res", "thermodynamics")

        # Make output directory if does not exist
        if isempty(OUT_DIR)
            @error "Output directory cannot be an empty path"
            return false
        end
        atmos.OUT_DIR = abspath(OUT_DIR)
        if samefile(atmos.OUT_DIR, atmos.ROOT_DIR)
            @error "Output directory cannot be the AGNI root directory"
            @error "    Got: $(atmos.OUT_DIR)"
            return false
        end
        if !isdir(atmos.OUT_DIR) && !isfile(atmos.OUT_DIR)
            mkdir(atmos.OUT_DIR)
        end
        @debug "Using OUT_DIR='$(atmos.OUT_DIR)'"

        # Directory used for fast I/O
        if IO_DIR == UNSET_STR
            # not set, so use output directory
            atmos.IO_DIR = OUT_DIR
        else
            # set, user provided
            if isempty(OUT_DIR)
                @error "I/O directory cannot be an empty path"
                return false
            end
            atmos.IO_DIR = abspath(IO_DIR)
            if !isdir(atmos.IO_DIR) && !isfile(atmos.IO_DIR)
                mkdir(atmos.IO_DIR)
            end
        end
        @debug "Using IO_DIR='$(atmos.IO_DIR)'"

        # Directory used for writing animation frames
        atmos.FRAMES_DIR  =  joinpath(atmos.IO_DIR, "frames")

        # -------------------------
        # Set other parameters
        # -------------------------

        # Set parameters for benchmarking
        atmos.benchmark   = false
        atmos.num_rt_eval = 0
        atmos.tim_rt_eval = 0.0

        # Work out which RT scheme is going to be used
        if strip(lowercase(spfile)) == "greygas"
            atmos.rt_scheme = RT_GREYGAS
            atmos.spectral_file = "greygas"
            atmos.SOCRATES_VERSION = "0000"
            @info "Using grey-gas radiative transfer scheme"

            # check options
            if flag_rayleigh || flag_cloud
                @error "Scattering not supported by grey-gas RT scheme!"
                @error "    In this case, disable rayleigh scattering and clouds"
                return false
            end

        else
            atmos.rt_scheme = RT_SOCRATES
            atmos.spectral_file = abspath(spfile)

            @debug "Using SOCRATES at $(ENV["RAD_DIR"])"

            # Get SOCRATES version
            atmos.SOCRATES_VERSION = readchomp(joinpath(ENV["RAD_DIR"],"version"))
            @debug "SOCRATES VERSION = "*atmos.SOCRATES_VERSION

            # Check SOCRATES version is valid
            if parse(Float64, atmos.SOCRATES_VERSION) < SOCVER_minimum
                @error "SOCRATES is out of date and cannot be used!"
                @error "    found at $(ENV["RAD_DIR"])"
                @error "    version is "*atmos.SOCRATES_VERSION
                return false
            end
        end

        # --------------------------------------
        # Allocate arrays
        # --------------------------------------


        # Allocate SOCRATES structs
        atmos.dimen =       SOCRATES.StrDim()
        atmos.control =     SOCRATES.StrCtrl()
        atmos.spectrum =    SOCRATES.StrSpecData()
        atmos.atm =         SOCRATES.StrAtm()
        atmos.cld =         SOCRATES.StrCld()
        atmos.aer =         SOCRATES.StrAer()
        atmos.bound =       SOCRATES.StrBound()
        atmos.radout =      SOCRATES.StrOut()

        atmos.all_channels =    all_channels
        atmos.overlap_method =  overlap_method

        atmos.real_gas      =   real_gas
        atmos.thermo_funct  =   thermo_functions

        atmos.tmp_floor =       max(1,tmp_floor)
        atmos.tmp_ceiling =     tmp_ceiling
        if atmos.tmp_ceiling <= atmos.tmp_floor+100
            @warn "The temperature limits provided are probably invalid"
            @warn "    Got: tmp_floor   = $(atmos.tmp_floor) K"
            @warn "    Got: tmp_ceiling = $(atmos.tmp_ceiling) K"
        end

        atmos.nlev_c         =  nlev_centre
        atmos.nlev_l         =  atmos.nlev_c + 1
        _check_range("Number of levels", atmos.nlev_c; min=NLEV_minimum) || return false

        atmos.tmp_surf =        max(tmp_surf, atmos.tmp_floor)
        atmos.grav_surf =       max(1.0e-7, gravity)
        atmos.surface_material= surface_material
        atmos.albedo_s =        max(min(albedo_s, 1.0 ), 0.0)
        atmos.albedo_b =        max(min(albedo_b,1.0), 0.0)
        atmos.zenith_degrees =  min(zenith_degrees, 89.9)
        _check_range("Zenith angle", atmos.zenith_degrees; min=0.0, max=90.0) || return false
        atmos.instellation =    instellation
        _check_range("Instellation", atmos.instellation; min=0) || return false
        atmos.s0_fact =         s0_fact
        _check_range("Stellar s0 factor", atmos.s0_fact; min=0, max=1) || return false

        atmos.toa_heating =     atmos.instellation * (1.0 - atmos.albedo_b) *
                                    s0_fact * cosd(atmos.zenith_degrees)

        atmos.flux_int =        flux_int
        atmos.target_olr =      max(1.0e-10, target_olr)

        atmos.phs_timescale =   phs_timescale
        atmos.evap_efficiency = max(min(evap_efficiency, 1.0),0.0)

        atmos.κ_grey_lw = κ_grey_lw
        _check_range("Grey LW opacity", atmos.κ_grey_lw; min=0) || return false
        atmos.κ_grey_sw = κ_grey_sw
        _check_range("Grey SW opacity", atmos.κ_grey_sw; min=0) || return false

        atmos.Kzz_floor =       max(0.0, Kzz_floor / 1e4)  # convert to SI units
        atmos.Kzz_ceiling =     1.0e20 / 1e4
        atmos.Kzz_pbreak =      1e5 # 1 bar as default location for break point
        atmos.Kzz_kbreak =      max(0.0, atmos.Kzz_floor)
        atmos.mlt_asymptotic =  mlt_asymptotic
        atmos.mlt_criterion =   mlt_criterion

        atmos.surf_roughness = surf_roughness
        _check_range("Surface roughness", atmos.surf_roughness; min=SURF_ROUGHNESS_MIN) || return false
        atmos.surf_windspeed = surf_windspeed
        _check_range("Surface windspeed", atmos.surf_windspeed; min=SURF_WINDSPEED_MIN) || return false
        atmos.C_d =  0.001  # placeholder, will be overwritten

        if atmos.real_gas && (atmos.mlt_criterion == 'l')
            @warn "Ledoux criterion not self-consistently supported for real gases"
            @warn "    (Will use Ledoux criterion anyway)"
            # @warn "    Switching criterion to Schwarzschild, neglecting MMW gradients"
            # atmos.mlt_criterion = 's'
        end
        if !(atmos.mlt_criterion in ['s','l'])
            @error "Invalid choice for mlt_criterion: $(atmos.mlt_criterion)"
            @error "    Must be: 's' or 'l' only"
            return false
        end

        atmos.tmp_magma =  max(atmos.tmp_floor, tmp_magma)
        atmos.skin_d = skin_d
        _check_range("Surface CBL thickness", atmos.skin_d; min=SKIN_D_MIN) || return false
        atmos.skin_k = skin_k
        _check_range("Surface CBL conductivity", atmos.skin_k; min=SKIN_K_MIN) || return false

        # pressure boundaries
        atmos.p_toa = p_top * 1.0e5 # Convert bar -> Pa
        atmos.p_boa = p_surf * 1.0e5
        if atmos.p_toa > atmos.p_boa
            @error "Top pressure must be less than surface pressure"
            @error "    p_top  = $p_top bar"
            @error "    p_surf = $p_surf bar"
            return false
        end
        atmos.p_oboa = atmos.p_boa

        # interior radius
        atmos.rp = radius
        _check_range("Planet surface radius", atmos.rp; min=RP_MIN) || return false

        # derived statistics
        atmos.interior_mass  =  atmos.grav_surf * atmos.rp^2 / phys.G_grav
        atmos.interior_rho   =  3.0 * atmos.interior_mass / ( 4.0 * pi * atmos.rp^3)
        atmos.transspec_μ    =  0.0
        atmos.transspec_rho  =  0.0
        atmos.transspec_tmp  =  0.0
        atmos.transspec_grav =  0.0
        atmos.transspec_r    =  0.0
        atmos.transspec_p    =  2e3     # 20 mbar = 2000 Pa
        if atmos.p_toa > atmos.transspec_p
            @error "p_top must be less than transspec_p"
            @error "    Got p_top:       $(atmos.p_toa) Pa"
            @error "    and transspec_p: $(atmos.transspec_p) Pa"
            return false
        end

        # absorption contributors
        atmos.control.l_gas::Bool =         true
        atmos.control.l_rayleigh::Bool =    flag_rayleigh
        atmos.control.l_continuum::Bool =   flag_continuum
        atmos.control.l_cont_gen::Bool =    flag_gcontinuum
        atmos.control.l_aerosol::Bool =     flag_aerosol
        atmos.control.l_cloud::Bool =       flag_cloud
        atmos.control.l_drop::Bool =        flag_cloud
        atmos.control.l_ice::Bool  =        false
        atmos.transparent =                 false

        # Initialise temperature grid
        atmos.tmpl = zeros(Float64, atmos.nlev_l)
        atmos.tmp =  zeros(Float64, atmos.nlev_c)

        # Initialise pressure grid with current p_toa and p_boa
        generate_pgrid!(atmos)

        # Initialise grid and other arrays
        #    radii
        atmos.r             = zeros(Float64, atmos.nlev_c) # radii at cell centres [m]
        atmos.rl            = zeros(Float64, atmos.nlev_l) # radii at cell edges [m]
        atmos.layer_thick   = zeros(Float64, atmos.nlev_c)
        #    gravity
        atmos.g             = ones(Float64, atmos.nlev_c) * atmos.grav_surf
        atmos.gl            = ones(Float64, atmos.nlev_l) * atmos.grav_surf
        #    enclosed mass [kg]
        atmos.m             = ones(Float64, atmos.nlev_c) * atmos.interior_mass
        atmos.ml            = ones(Float64, atmos.nlev_l) * atmos.interior_mass
        #    surface density (mass of each layer per unit area [kg m-2])
        atmos.layer_σ       = zeros(Float64, atmos.nlev_c)
        atmos.layer_isbound = trues(atmos.nlev_c)

        # Initialise thermodynamic properties
        atmos.layer_μ       = zeros(Float64, atmos.nlev_c)
        atmos.layer_ρ       = zeros(Float64, atmos.nlev_c)
        atmos.layer_cp      = zeros(Float64, atmos.nlev_c)
        atmos.layer_kc      = zeros(Float64, atmos.nlev_c)

        # Initialise cloud arrays
        atmos.cloud_arr_r   = zeros(Float64, atmos.nlev_c)
        atmos.cloud_arr_l   = zeros(Float64, atmos.nlev_c)
        atmos.cloud_arr_f   = zeros(Float64, atmos.nlev_c)

        # Phase change timescales [seconds]
        atmos.phs_timescale = phs_timescale
        _check_range("Phase change timescale", atmos.phs_timescale; min=PHS_TIMESCALE_MIN) || return false

        # Evaporation efficiency
        atmos.evap_efficiency = evap_efficiency
        _check_range("Evaporation efficiency", atmos.evap_efficiency; min=0, max=1) || return false

        # Hardcoded cloud properties
        atmos.cloud_alpha   = 0.01    # 1% of condensed water forms substantial clouds
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

        # Gas variables
        atmos.gas_num   =   0                                 # number of gases
        atmos.gas_names =   Array{String}(undef, 0)           # list of names
        atmos.gas_dat =     Dict{String, phys.Gas_t}()        # dict of gas data structs
        atmos.gas_vmr  =    Dict{String, Array{Float64,1}}()  # dict of VMR arrays
        atmos.gas_ovmr  =   Dict{String, Array{Float64,1}}()  # ^ backup of initial values
        atmos.gas_cvmr  =   Dict{String, Array{Float64,1}}()  # ^ backup of initial values
        atmos.gas_safe  =   Dict{String, Bool}()              # considered 'safe'?
        atmos.gas_sat  =    Dict{String, Array{Bool, 1}}()    # mask for saturation

        # Metallicities
        atmos.metal_orig =  Dict{String, Float64}()          # input metallicities rel to H
        atmos.metal_calc =  Dict{String, Float64}()          # calculated metallicities (empty for now)
        for k in keys(metallicities)
            # mass -> mole, by scaling factor 1/mu
            atmos.metal_orig[k] = metallicities[k] * phys._get_mmw("H") / phys._get_mmw(k)

            if atmos.metal_orig["H"] < 1e-30
                @error "Cannot define metallicity of hydrogen relative to itself!"
            end
        end

        # Phase change compositional variables
        #    condensation aloft
        atmos.condensates = Array{String}(undef, 0)           # list of condensates
        atmos.cond_yield =  Dict{String, Array{Float64,1}}()  # cond/evap yield at each layer
        atmos.cond_accum =  Dict{String, Float64}()           # sum of each yield in atmosphere
        #    ocean reservoirs
        atmos.ocean_ini =   Dict{String,Float64}()            # initial ocean reservoir [kg/m^2]
        atmos.ocean_tot =   Dict{String,Float64}()            # final ocean reservoir [kg/m^2]
        #    ocean layers
        atmos.ocean_maxdepth  = 0.0
        atmos.ocean_areacov   = 0.0
        atmos.ocean_topliq    = UNSET_STR
        atmos.ocean_layers    = Tuple[(1,UNSET_STR,0.0,0.0),]  # array of tuples


        # Dict input case
        if mf_source == 0
            @debug "Composition set by dict"
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

        # CSV file case
        elseif mf_source == 1
            # check file
            if !isfile(mf_path)
                @error "Could not read VMR file '$mf_path'"
                return false
            end
            @debug "Composition set by file"

            # get header
            mf_head::String =   readline(abspath(mf_path))

            # remove comment symbol at start
            mf_head = mf_head[2:end]

            # remove whitespace
            mf_head = replace(mf_head, " " => "")

            # split by column and drop first three
            heads::Array{String,1} = split(mf_head, ",")[1:end]

            # get body
            mf_body::Array{Float64,2} = readdlm(abspath(mf_path), ',', Float64;
                                                header=false, skipstart=2)
            mf_body = transpose(mf_body)

            # Get pressure array from file
            arr_p::Array{Float64,1} = mf_body[1,:]
            arr_x::Array{Float64,1} = zero(arr_p)

            # set composition by interpolating with pressure array
            # the pressure array must be the first column in the file
            for li in 2:lastindex(heads)
                gas = strip(heads[li], [' ','\t','\n'])

                # Check if this gas is repeated, or contains invalid chars
                if gas in atmos.gas_names
                    @warn "VMR file contains repeated gas '$gas'"
                    continue

                elseif !occursin(r"^[[:alnum:]]*$", gas)
                    @warn "VMR file contains invalid gas '$gas'"
                    continue

                else
                    # Store with zero VMR
                    atmos.gas_vmr[gas] = zeros(Float64, atmos.nlev_c)
                    push!(atmos.gas_names, gas)
                    @debug "    added $gas"
                    atmos.gas_num += 1
                end

                # Get VMR values from file
                arr_x[:] .= mf_body[li,:]

                # Extend loaded profile to lower pressures (prevents domain error)
                if arr_p[1] > atmos.p_toa
                    pushfirst!(arr_p,   atmos.p_toa/1.1)
                    pushfirst!(arr_x,   arr_x[1] )
                end

                # Extend loaded profile to higher pressures (prevents domain error)
                if arr_p[end] < atmos.p_boa
                    push!(arr_p, atmos.p_boa*1.1)
                    push!(arr_x, arr_x[end])
                end

                # Set up interpolator using file data
                itp = interpolate((log10.(arr_p),),arr_x, Gridded(Linear()))

                # Set values in atmos struct
                for i in 1:atmos.nlev_c
                    atmos.gas_vmr[gas][i] = itp(log10(atmos.p[i]))
                end
            end

        else
            @error "Invalid value for setting mole fraction source"
            @error "    Got mf_source=$mf_source"
            return false
        end

        # add extra gases if required
        if use_all_gases
            for gas in phys.gases_standard
                if !(gas in atmos.gas_names)
                    atmos.gas_vmr[gas] = zeros(Float64, atmos.nlev_c)
                    push!(atmos.gas_names, gas)
                    atmos.gas_num += 1
                end
            end
        else
            if !isempty(atmos.metal_orig)
                @error "Must set `use_all_gases=true` if providing metallicities"
                return false
            end
        end

        # backup mixing ratios from current state
        for k in keys(atmos.gas_vmr)
            atmos.gas_ovmr[k] = deepcopy(atmos.gas_vmr[k])
            atmos.gas_cvmr[k] = deepcopy(atmos.gas_vmr[k])
        end

        # set phase change quantities for each gas
        for g in atmos.gas_names
            # layers at which it is saturated
            atmos.gas_sat[g]    = falses(atmos.nlev_c)

            # production/removal of condensate at each layer
            atmos.cond_yield[g] = zeros(Float64, atmos.nlev_c)

            # total yield of condensate from rain/evap within atmosphere
            #     to be calculated as: sum(cond_yield[g])
            atmos.cond_accum[g] = 0.0

            # amount of condensate at the surface
            atmos.ocean_ini[g]       = 0.0  # TODO: make this input from user
            atmos.ocean_tot[g]       = 0.0  # to be calculated
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
        @info "Loading thermodynamic data"
        for g in atmos.gas_names
            atmos.gas_dat[g] = phys.load_gas(atmos.THERMO_DIR, g,
                                                atmos.thermo_funct, atmos.real_gas;
                                                check_integrity=check_integrity)

            if atmos.gas_dat[g].fail
                @error "Problem when loading thermodynamic data"
                @error "    Try downloading them again and/or updating AGNI."
                return false
            end
        end

        # store condensates
        for c in condensates
            if atmos.gas_dat[c].stub || atmos.gas_dat[c].no_sat || (c in COND_DISALLOWED)
                @warn "$c disallowed from being condensable; treated as dry"
            else
                push!(atmos.condensates, c)
            end
        end

        # Validate condensate names
        atmos.condense_any = length(condensates) > 0
        if atmos.condense_any
            for c in condensates
                if !(c in atmos.gas_names)
                    @error "Invalid condensate. '$c' is not in gas list"
                    return false
                end
            end
        end

        # Ocean params
        atmos.ocean_cs_height = max(0.0, ocean_cs_height)
        _check_range("Continent shelf height", atmos.ocean_cs_height; min=0) || return false
        atmos.ocean_ob_frac = ocean_ob_frac
        _check_range("Ocean basin fraction", atmos.ocean_ob_frac; min=0, max=1) || return false

        # Set initial temperature profile to a small value which still keeps
        #   all of the gases supercritical. This should be a safe condition to
        #   default to, although the user must specify a profile in the cfg.
        for g in atmos.gas_names
            atmos.tmpl[end] = max(atmos.tmpl[end], atmos.gas_dat[g].T_crit+5.0)
            fill!(atmos.tmpl, atmos.tmpl[end])
            fill!(atmos.tmp, atmos.tmpl[end])
        end

        # Check T,P range vs EOS limits
        for g in atmos.gas_names
            if atmos.p_boa > atmos.gas_dat[g].prs_max
                @warn "Surface pressure exceeds the valid range ($g EOS)"
            end
            if maximum(atmos.tmp) > atmos.gas_dat[g].tmp_max
                @warn "Temperature profile exceeds the valid range ($g EOS)"
            end
        end

        # Fastchem directory
        atmos.flag_fastchem = false
        if fastchem_work == UNSET_STR
            # default
            atmos.fastchem_work = joinpath(atmos.IO_DIR, "fastchem")  # default path
            @debug "Fastchem working dir defaulting to $(atmos.fastchem_work)"
        else
            # user-provided
            atmos.fastchem_work = abspath(fastchem_work)
            @debug "Fastchem working dir set to $(atmos.fastchem_work)"
        end
        if samefile(atmos.fastchem_work, atmos.ROOT_DIR)
            @error "FastChem working directory cannot be the AGNI root directory"
            return false
        end

        # Fastchem enabled by environment?
        if "FC_DIR" in keys(ENV)

            @debug "FastChem env has been set"

            # check fastchem installation folder
            atmos.FC_DIR = abspath(ENV["FC_DIR"])
            if !isdir(atmos.FC_DIR)
                @error "Could not find FastChem installed at FC_DIR='$(atmos.FC_DIR)'"
                @error "    Install FastChem with `\$ ./src/get_fastchem.sh`"
                return false
            end

            # check executable
            atmos.fastchem_exec = abspath(atmos.FC_DIR,"fastchem")
            atmos.flag_fastchem = isfile(atmos.fastchem_exec)
            if !atmos.flag_fastchem
                @error "Could not find fastchem executable inside '$(atmos.FC_DIR)'"
                @error "    Install FastChem with `\$ ./src/get_fastchem.sh`"
                return false
            else
                @debug "Found FastChem executable at $(atmos.fastchem_exec)"
            end

            # re-make FC working directory
            rm(atmos.fastchem_work,force=true,recursive=true)
            mkdir(atmos.fastchem_work)
        else
            @debug "FastChem env variable not set, so FC won't be available for use"
        end
        # other parameters for FC
        atmos.fastchem_wellmixed    = fastchem_wellmixed
        atmos.fastchem_floor        = fastchem_floor
        atmos.fastchem_maxiter_chem = fastchem_maxiter_chem
        atmos.fastchem_maxiter_solv = fastchem_maxiter_solv
        atmos.fastchem_xtol_chem    = fastchem_xtol_chem
        atmos.fastchem_xtol_elem    = fastchem_xtol_elem
        atmos.fastchem_conf         = joinpath(atmos.fastchem_work,"config.input")
        atmos.fastchem_elem         = joinpath(atmos.fastchem_work,"elements.dat")
        atmos.fastchem_chem         = joinpath(atmos.fastchem_work,"chemistry.dat")
        atmos.fastchem_cond         = joinpath(atmos.fastchem_work,"condensates.dat")
        atmos.fastchem_prof         = joinpath(atmos.fastchem_work,"pt.dat")
        atmos.fastchem_moni         = joinpath(atmos.fastchem_work,"monitor.dat")

        # RFM
        atmos.flag_rfm = !(rfm_parfile == UNSET_STR)
        atmos.rfm_work = joinpath(atmos.IO_DIR, "rfm")
        if samefile(atmos.rfm_work, atmos.ROOT_DIR)
            @error "RFM working directory cannot be the AGNI root directory"
            return false
        end
        if atmos.flag_rfm
            atmos.rfm_parfile = abspath(rfm_parfile)
            @debug "RFM parfile set: $(atmos.rfm_parfile)"
            if !isfile(atmos.rfm_parfile)
                @error "Could not find parfile provided by user, for running RFM"
                @error "    atmos.rfm_parfile=$(atmos.rfm_parfile)"
                return false
            end

            # re-make working directory
            rm(atmos.rfm_work,force=true,recursive=true)
            mkdir(atmos.rfm_work)
        end

        # Deep atmospheric heating (default: disabled)
        atmos.deep_heating = DeepHeatingParams()

        # Record that the parameters are set
        atmos.is_param = true
        atmos.is_solved = false
        atmos.is_converged = false

        @debug "Setup complete"
        return true
    end # end function setup


    """
    **Set deep atmospheric heating parameters.**

    Configures the deep heating source
    which deposits energy as a Gaussian distribution in log-pressure space.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    - `active::Bool`            enable/disable deep heating.
    - `P_dep::Float64`          deposition pressure [Pa].
    - `sigma_P::Float64`        width of Gaussian in log-pressure space [dimensionless].
    - `efficiency::Float64`     heating efficiency as fraction of instellation.

    Returns:
        Nothing
    """
    function set_deep_heating!(atmos::atmosphere.Atmos_t;
                                active::Bool=false,
                                P_dep::Float64=1.0e5,
                                sigma_P::Float64=1.0,
                                efficiency::Float64=0.0,
                                normalization::Symbol=:pressure,
                                below_domain::Symbol=:clamp,
                                power_mode::Symbol=:efficiency,
                                F_total::Float64=0.0)

        # Validate/normalise parameters
        sigma_P = max(sigma_P, 0.01)  # Prevent division issues
        efficiency = clamp(efficiency, 0.0, 1.0)

        # Normalise symbols
        if !(normalization in (:pressure, :mass))
            error("Invalid deep heating normalization: $(normalization)")
        end
        if !(below_domain in (:clamp, :boundary_flux))
            error("Invalid deep heating below_domain: $(below_domain)")
        end
        if !(power_mode in (:efficiency, :flux))
            error("Invalid deep heating power_mode: $(power_mode)")
        end

        # Deposition pressure handling
        if below_domain == :clamp
            P_dep = clamp(P_dep, atmos.p_toa, atmos.p_boa)
        else
            # Allow P_dep outside domain (used as a signal for boundary-flux mode)
            P_dep = max(P_dep, atmos.p_toa)
        end

        # Create and store the parameters
        atmos.deep_heating = DeepHeatingParams(active, P_dep, sigma_P, efficiency,
                               normalization, below_domain,
                               power_mode, F_total)

        if active
            @info @sprintf("Deep heating enabled: norm=%s below=%s mode=%s",
                           String(normalization), String(below_domain), String(power_mode))
            if power_mode == :efficiency
                @info @sprintf("    P_dep=%.2e Pa, σ_P=%.2f, ε=%.4f", P_dep, sigma_P, efficiency)
            else
                @info @sprintf("    P_dep=%.2e Pa, σ_P=%.2f, F_total=%.4e W/m²", P_dep, sigma_P, F_total)
            end
        end

        return nothing
    end


    """
    **Estimate photosphere.**

    Estimates the location of the photosphere by finding the median of the contribution
    function in each band (0.2 um to 150 um), and then finding the pressure level at which
    these median values are maximised.

    Arguments:
        - `atmos::Atmos_t`          the atmosphere struct instance to be used.

    Returns:
        - `p_ref::Float64`          pressure level of photosphere [Pa]
    """
    function estimate_photosphere!(atmos::atmosphere.Atmos_t)::Float64

        # Params
        wl_min::Float64  = 0.2 * 1e-6 # 200 nm
        wl_max::Float64  = 150 * 1e-6 # 150 um
        p_min::Float64   = 10.0 # 1e-4 bar

        # tracking
        cff_max::Float64 = 0.0
        cff_try::Float64 = 0.0
        atmos.transspec_p = p_min

        # get band indices
        wl_imin = findmin(abs.(atmos.bands_cen .- wl_min))[2]
        wl_imax = findmin(abs.(atmos.bands_cen .- wl_max))[2]

        # reversed?
        if wl_imin > wl_imax
            wl_imin, wl_imax = wl_imax, wl_imin
        end

        # loop over levels
        for i in 1:atmos.nlev_c
            if atmos.p[i] < p_min
                continue
            end

            # maximum contfunc in this band
            cff_try = Statistics.median(atmos.contfunc_band[i,wl_imin:wl_imax])

            # is this more than the existing maximum?
            if cff_try > cff_max
                cff_max = cff_try
                atmos.transspec_p = atmos.p[i]
            end
        end

        return atmos.transspec_p
    end

    """
    **Calculate observed radius and bulk density.**

    This is done at the layer probed in transmission, which is set to a fixed pressure.

    Arguments:
        - `atmos::Atmos_t`          the atmosphere struct instance to be used.

    Returns:
        - `transspec_rho::Float64`  the bulk density observed in transmission
    """
    function calc_observed_rho!(atmos::atmosphere.Atmos_t)::Float64

        # transspec_r::Float64            # planet radius probed in transmission [m]
        # transspec_m::Float64            # mass [kg] of atmosphere + interior
        # transspec_rho::Float64          # bulk density [kg m-3] implied by r and m

        # Store reference pressure in atmos struct
        # estimate_photosphere!(atmos)

        # get the observed height
        idx::Int = findmin(abs.(atmos.p .- atmos.transspec_p))[2]
        atmos.transspec_r    = atmos.r[idx]
        atmos.transspec_μ    = atmos.layer_μ[idx]
        atmos.transspec_tmp  = atmos.tmp[idx]
        atmos.transspec_grav = atmos.g[idx]

        # get mass of whole atmosphere, assuming hydrostatic
        atmos.transspec_m = atmos.p_boa * 4 * pi * atmos.rp^2 / atmos.grav_surf

        # add mass of the interior component
        atmos.transspec_m += atmos.interior_mass

        # get density of all enclosed by observed layer
        # store this in the atmosphere struct
        atmos.transspec_rho = 3.0 * atmos.transspec_m / (4.0 * pi * atmos.transspec_r^3)

        # also return the value
        return atmos.transspec_rho
    end

    """
    **Get pressure at top and bottom of convective zone**

    Arguments:
        - `atmos::Atmos_t`          the atmosphere struct instance to be used.

    Returns:
        - `p_top::Float64`          pressure [Pa] at top of convective zone
        - `p_bot::Float64`          pressure [Pa] at bottom of convective zone
    """
    function estimate_convective_zone(atmos::atmosphere.Atmos_t)::Tuple{Float64,Float64}

        # Defaults to zero, if there's no convection
        p_top::Float64 = 0.0
        p_bot::Float64 = 0.0

        # Loop from top-down to find p_top
        for i in 1:atmos.nlev_l
            if atmos.mask_c[i]
                p_top = atmos.pl[i]
                break
            end
        end

        # Loop from bottom-up to find p_bot
        for i in range(start=atmos.nlev_l, stop=1, step=-1)
            if atmos.mask_c[i]
                p_bot = atmos.pl[i]
                break
            end
        end

        # Return top, bot
        return (p_top, p_bot)
    end

    """
    **Calculate properties within each layer of the atmosphere (e.g. density, mmw).**

    Function will return false if hydrostatic calculcation fails. This is usually when
    the atmosphere becomes unbound.

    Arguments:
        - `atmos::Atmos_t`          the atmosphere struct instance to be used.

    Returns:
        - `ok::Bool`                function result is ok
    """
    function calc_layer_props!(atmos::atmosphere.Atmos_t)::Bool
        if !atmos.is_param
            @error("Atmosphere struct has not been setup")
            return false
        end

        # Status
        ok::Bool = true

        # MMW
        calc_profile_mmw!(atmos)

        # Heat capacity and thermal conductivity
        calc_profile_cpkc!(atmos)

        # Set density at each level
        calc_profile_density!(atmos)

        # Perform hydrostatic integration
        ok = ok && calc_profile_radius!(atmos)

        return ok
    end

    """
    **Calculate radii, gravities, and masses for all layers.**

    Performs hydrostatic integration from the ground upwards.
    Requires density, temperature, pressure to have already been set.

    Does not account for surface ocean height.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.

    Returns:
    - `bound::Bool`             atmosphere is strongly bound by gravity
    """
    function calc_profile_radius!(atmos::atmosphere.Atmos_t)::Bool

        # Reset arrays
        fill!(atmos.r         ,    atmos.rp)
        fill!(atmos.rl        ,    atmos.rp)
        fill!(atmos.g         ,    atmos.grav_surf)
        fill!(atmos.gl        ,    atmos.grav_surf)
        fill!(atmos.m         ,    atmos.interior_mass)
        fill!(atmos.ml        ,    atmos.interior_mass)
        fill!(atmos.layer_thick,   1.0)
        fill!(atmos.layer_σ ,      1.0)
        fill!(atmos.layer_isbound, true)

        # Check config...
        if HYDROGRAV_constg && HYDROGRAV_selfg
            @warn "Incompatible gravity parameters have been set:"
            @warn "    constant with height (HYDROGRAV_constg=$HYDROGRAV_constg)"
            @warn "    atmos self-attraction (HYDROGRAV_selfg=$HYDROGRAV_selfg)"
        end

        # Temporary values
        nsub::Int = round(Int, HYDROGRAV_steps/atmos.nlev_c, RoundUp)

        # Integrate from surface upwards
        for i in range(start=atmos.nlev_c, stop=1, step=-1)

            # ------------
            # Integrate from lower edge to centre
            atmos.r[i], atmos.g[i], atmos.m[i] =
                integ_hydrograv(atmos.rl[i+1], atmos.gl[i+1], atmos.ml[i+1], atmos.pl[i+1],
                                    atmos.p[i], atmos.layer_ρ[i], nsub)

            #   apply radius limiter
            atmos.r[i] = max(atmos.r[i], atmos.rl[i+1] + HYDROGRAV_mindr)
            if atmos.r[i] > atmos.rl[i+1] + HYDROGRAV_maxdr/2
                atmos.r[i] = atmos.rl[i+1] + HYDROGRAV_maxdr/2
                atmos.layer_isbound[i] = false
            end

            #   apply gravity limiter
            if HYDROGRAV_constg
                atmos.g[i] = atmos.grav_surf
            end
            if atmos.g[i] < HYDROGRAV_ming
                atmos.g[i] = HYDROGRAV_ming
                atmos.layer_isbound[i] = false
            end

            # ------------
            # Integrate from centre to upper edge
            atmos.rl[i], atmos.gl[i], atmos.ml[i] =
                integ_hydrograv(atmos.r[i], atmos.g[i], atmos.m[i], atmos.p[i],
                                    atmos.pl[i], atmos.layer_ρ[i], nsub)

            #   apply radius limiter
            atmos.rl[i] = max(atmos.rl[i], atmos.r[i] + HYDROGRAV_mindr)
            if atmos.rl[i] > atmos.r[i] + HYDROGRAV_maxdr/2
                atmos.rl[i] = atmos.r[i] + HYDROGRAV_maxdr/2
                atmos.layer_isbound[i] = false
            end

            #   apply gravity limiter
            if HYDROGRAV_constg
                atmos.gl[i] = atmos.grav_surf
            end
            if atmos.gl[i] < HYDROGRAV_ming
                atmos.gl[i] = HYDROGRAV_ming
                atmos.layer_isbound[i] = false
            end

            # Store: Layer geometrical thickness [m]
            atmos.layer_thick[i] = atmos.rl[i] - atmos.rl[i+1]

            # Mass of layer, per unit area at layer-centre [kg m-2]
            atmos.layer_σ[i] = (atmos.ml[i] - atmos.ml[i+1])/(4 * pi * atmos.r[i]^2)
        end

        return all(atmos.layer_isbound)
    end

    """
    **Integrate hydrostatic and gravity equations across a pressure interval.**

    Uses the classic fourth-order Runge-Kutta method.

    Arguments:
    - `r0::Float64`     radius   at start of interval [m]
    - `g0::Float64`     gravity  at start of interval [m s-2]
    - `m0::Float64`     mass enc at start of interval [kg]
    - `p0::Float64`     pressure at start of interval [Pa]
    - `p1::Float64`     pressure at end   of interval [Pa]
    - `rho::Float64`    density throughout interval, constant [kg m-3]
    - `n::Int`          number of steps for integration (n >= 2)

    Returns:
    - `rj::Float64`     radius   at end of interval [m]
    - `gj::Float64`     gravity  at end of interval [kg]
    - `mj::Float64`     mass enc at end of interval [kg]
    """
    function integ_hydrograv(r0::Float64, g0::Float64, m0::Float64, p0::Float64,
                                    p1::Float64, rho::Float64, n::Int)::Tuple{Float64,Float64,Float64}

        # Work variables
        pj::Float64 = p0    # rolling pressure (decreasing)
        rj::Float64 = r0    # rolling radius   (increasing)
        gj::Float64 = g0    # rolling gravity  (incr, decr, or constant)
        mj::Float64 = m0    # rolling mass     (increasing)

        # Get gravity at r
        function _grav(r)
            if HYDROGRAV_constg
                # gravity is constant
                return g0
            else
                if HYDROGRAV_selfg
                    # gravity changes with mass and radius
                    return phys.grav_accel(mj, rj)
                else
                    # gravity changes with radius only
                    return g0 * (r0/rj)^2
                end
            end
        end

        # Derivative to integrate
        #   dr/dp = -1 / (rho * g(r))
        function _drdp(p,r)
            return -1 / (rho * _grav(r))
        end

        # Parameters
        dp::Float64  = (p1-p0)/max(2,n) # this will be negative
        dp2::Float64 = dp/2
        k1::Float64  = 0.0; k2::Float64 = 0.0
        k3::Float64  = 0.0; k4::Float64 = 0.0

        # Loop over sub-levels between p0 and p1
        for _ in range(p0, stop=p1, step=dp)

            # Integrate radius ...
            k1 = _drdp(pj,       rj)
            k2 = _drdp(pj + dp2, rj + k1*dp2)
            k3 = _drdp(pj + dp2, rj + k2*dp2)
            k4 = _drdp(pj + dp,  rj + k3*dp)
            rj += dp/6 * (k1 + 2*k2 + 2*k3 + k4)

            # Integrate mass enclosed ...
            mj += 4 * pi * rj^2 * (-1 * dp) / gj

            # Integrate pressure (negative change )
            pj += dp

            # Update gravity
            gj = _grav(rj)
        end

        return (rj, gj, mj)
    end

    """
    **Calculate mean molecular weight for all layers.**

    MMW stored as kg/mol.

    Arguments:
        - `atmos::Atmos_t`      the atmosphere struct instance to be used.
    """
    function calc_profile_mmw!(atmos::atmosphere.Atmos_t)
        fill!(atmos.layer_μ, 0.0)
        for gas in atmos.gas_names
            @. atmos.layer_μ += atmos.gas_vmr[gas] * atmos.gas_dat[gas].mmw
        end
        return nothing
    end

    """
    **Calculate specific heat capacity and thermal conductivity for all layers.**

    Specific heat per unit mass: J K-1 kg-1.
    Thermal conductivity: W m-1 K-1.

    Arguments:
        - `atmos::Atmos_t`      the atmosphere struct instance to be used.
    """
    function calc_profile_cpkc!(atmos::atmosphere.Atmos_t)
        # Loop over layers
        @inbounds for i in 1:atmos.nlev_c
            calc_single_cpkc!(atmos, i)
        end
        return nothing
    end

    """
    **Calculate specific heat capacity and thermal conductivity of a single layer.**

    Specific heat per unit mass: J K-1 kg-1.
    Thermal conductivity: W m-1 K-1.

    Arguments:
        - `atmos::Atmos_t`      the atmosphere struct instance to be used.
        - `idx::Int`            index of the layer
    """
    function calc_single_cpkc!(atmos::atmosphere.Atmos_t, idx::Int)
        # Reset
        mmr::Float64 = 0.0
        atmos.layer_cp[idx] = 0.0
        atmos.layer_kc[idx] = 0.0
        # Loop over gases
        for gas in atmos.gas_names
            mmr = atmos.gas_vmr[gas][idx] * atmos.gas_dat[gas].mmw/atmos.layer_μ[idx]
            atmos.layer_cp[idx] += mmr * phys.get_Cp(atmos.gas_dat[gas], atmos.tmp[idx])
            atmos.layer_kc[idx] += mmr * phys.get_Kc(atmos.gas_dat[gas], atmos.tmp[idx])
        end
        return nothing
    end

    """
    **Calculate the mass-density for all layers.**

    Requires temperature, pressure, mmw to be already have been set.

    Arguments:
        - `atmos::Atmos_t`      the atmosphere struct instance to be used.
    """
    function calc_profile_density!(atmos::atmosphere.Atmos_t)
        # Loop over levels
        @inbounds for i in 1:atmos.nlev_c
            calc_single_density!(atmos, i)
        end
        return nothing
    end

    """
    **Calculate the mass-density of a single layer.**

    Requires temperature, pressure, mmw to be already have been set.

    Arguments:
    - `atmos::Atmos_t`      the atmosphere struct instance to be used.
    - `idx::Int`            index of the layer
    """
    function calc_single_density!(atmos::atmosphere.Atmos_t, idx::Int)
        # Array of gas properties
        gas_arr = [atmos.gas_dat[gas]      for gas in atmos.gas_names]
        vmr_arr = [atmos.gas_vmr[gas][idx] for gas in atmos.gas_names]

        # Evaluate density
        atmos.layer_ρ[idx] = phys.calc_rho_mix(
                                                    gas_arr, vmr_arr,
                                                    atmos.tmp[idx], atmos.p[idx],
                                                    atmos.layer_μ[idx]
                                                    )
        return nothing
    end

    """
    **Generate pressure grid.**

    Almost equally log-spaced between p_boa and p_boa. The near-boundary layers
    are smaller than they would be on an equally log-spaced grid, to avoid numerics.

    Arguments:
    - `atmos::Atmos_t`              the atmosphere struct instance to be used.
    """
    function generate_pgrid!(atmos::atmosphere.Atmos_t)

        # Ensure pressures are in bounds
        if atmos.p_boa/atmos.p_toa < PRESSURE_RATIO_MIN
            @warn "Bottom/top pressure ratio is too small"
            @debug "   Got pressure ratio bot/top = $(atmos.p_boa/atmos.p_toa)"
            atmos.p_boa = PRESSURE_RATIO_MIN * atmos.p_toa
        end

        # Allocate arrays
        atmos.p  = zeros(Float64, atmos.nlev_c)
        atmos.pl = zeros(Float64, atmos.nlev_l)

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
        p_mid::Float64 = atmos.pl[end-1]*PRESSURE_FACT_BOT + atmos.pl[end-2]*(1.0-PRESSURE_FACT_BOT)
        atmos.pl[1:end-2] .= collect(Float64, range( start=atmos.pl[1],
                                                            stop=p_mid,
                                                            length=atmos.nlev_l-2))

        # Set cell-centres at midpoint of cell-edges
        atmos.p[1:end] .= 0.5 .* (atmos.pl[1:end-1] .+ atmos.pl[2:end])

        # Shrink top-most layer to avoid doing too much extrapolation
        # atmos.p[1] = atmos.pl[1]*PRESSURE_FACT_TOP + atmos.p[1]*(1-PRESSURE_FACT_TOP)

        # Finally, convert arrays to actual pressure units [Pa]
        @. atmos.p  = 10.0 ^ atmos.p
        @. atmos.pl = 10.0 ^ atmos.pl

        # Ensure pressure grid is strictly decreasing (from surface upwards)
        for i in range(start=atmos.nlev_c, stop=1, step=-1)
            atmos.p[i]  = min(atmos.p[i],  atmos.pl[i+1]/PRESSURE_RATIO_MIN)
            atmos.pl[i] = min(atmos.pl[i], atmos.p[i]/PRESSURE_RATIO_MIN)
        end

        return nothing
    end


    """
    **Allocate atmosphere arrays, prepare spectral files, and final steps.**

    Will not modify spectral file if `stellar_spectrum` is an empty string.
    Will treat star as blackbody with photospheric effective temperature of `stellar_Teff`
    if the parameter `stellar_spectrum` has value of `"blackbody"`.

    Arguments:
    - `atmos::Atmos_t`             the atmosphere struct instance to be used.
    - `stellar_spectrum::String`   path to stellar spectrum csv file
    - `stellar_Teff::Float64`      star effective temperature if blackbody
    - `check_safe_gas::Bool`       require that there be at least one 'safe' gas in the mix
    """
    function allocate!(atmos::atmosphere.Atmos_t, stellar_spectrum::String;
                        stellar_Teff::Float64=-1.0,
                        check_safe_gas::Bool=true)::Bool

        @debug "Allocate atmosphere"
        if !atmos.is_param
            @error "Atmosphere struct has not been setup"
            return false
        end

        atmos.atm.n_profile = 1

        #########################################
        # spectral data
        #########################################

        if atmos.rt_scheme == RT_GREYGAS

            atmos.star_file = "_CALC_AS_BLACKBODY=$stellar_Teff"
            atmos.nbands    = 1
            atmos.bands_max = Float64[1.0,]  # m
            atmos.bands_min = Float64[1e-9,] # m
            atmos.bands_wid = atmos.bands_max - atmos.bands_min
            atmos.bands_cen = Float64[1e-6,]

        elseif atmos.rt_scheme == RT_SOCRATES

            # Validate files
            if !isfile(atmos.spectral_file)
                @error "Spectral file not found"
                @error "    Got: '$(atmos.spectral_file)'"
                @error "To download, try using the `get_data.sh` script"
                @error "    e.g. to get CodenameXX: \$ ./src/get_data.sh anyspec Codename XX"
                return false
            end

            spectral_file_run::String  = joinpath([atmos.IO_DIR, "runtime.sf"])
            spectral_file_runk::String = joinpath([atmos.IO_DIR, "runtime.sf_k"])

            # Setup spectral file
            socstar::String = joinpath([atmos.IO_DIR, "socstar.dat"])
            if !isempty(stellar_spectrum)
                @debug "Inserting stellar spectrum into spectral file"

                # Remove if already exists
                rm(spectral_file_run , force=true)
                rm(spectral_file_runk, force=true)

                # Wl and Fl arrays, dummy values at TOA
                wl::Array{Float64,1} = Float64[1.0, 2.0, 3.0]  # nm
                fl::Array{Float64,1} = Float64[0.0, 0.0, 0.0]  # erg s-1 cm-2 nm-1

                # Insert stellar spectrum
                if lowercase(stellar_spectrum) == "blackbody"
                    # generate blackbody spectrum
                    @debug "Inserting spectrum as blackbody with Teff=$stellar_Teff K"

                    if stellar_Teff < 1.0
                        @error "Invalid stellar photospheric temperature: $stellar_Teff K"
                        @error "    Choose a valid temperature when using blackbody spectrum"
                        return false
                    end
                    atmos.star_file = "_CALC_AS_BLACKBODY=$stellar_Teff"

                    wl, fl = spectrum.blackbody_star(stellar_Teff, atmos.instellation)
                else
                    # use spectrum on disk
                    @debug "Inserting stellar spectrum from file $stellar_spectrum"

                    if !isfile(stellar_spectrum)
                        @error "Stellar spectrum file '$(stellar_spectrum)' does not exist"
                        @error "    Try using: \$ ./src/get_data.sh stellar"
                        return false
                    end
                    atmos.star_file = abspath(stellar_spectrum)

                    wl, fl = spectrum.load_from_file(atmos.star_file)
                end

                # Write stellar spectrum to disk in format required by SOCRATES
                spectrum.write_to_socrates_format(wl, fl, socstar) || return false

                # Insert stellar spectrum and rayleigh scattering, if required
                spectrum.insert_stellar_and_rscatter(atmos.spectral_file,
                                                        socstar, spectral_file_run,
                                                        atmos.control.l_rayleigh)

            else
                # Stellar spectrum was not provided, which is taken to mean that
                #       the spectral file includes it already.
                @info "Using pre-existing spectral file without modifications"
                atmos.star_file = "_ALREADY_IN_SPECTRAL_FILE"
                spectral_file_run  = atmos.spectral_file
                spectral_file_runk = atmos.spectral_file*"_k"
            end

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
            atmos.bands_max = zeros(Float64, atmos.nbands)
            atmos.bands_min = zeros(Float64, atmos.nbands)
            atmos.bands_cen = zeros(Float64, atmos.nbands)
            atmos.bands_wid = zeros(Float64, atmos.nbands)

            for i in 1:atmos.nbands
                atmos.bands_min[i] = atmos.spectrum.Basic.wavelength_short[i]
                atmos.bands_max[i] = atmos.spectrum.Basic.wavelength_long[i]
            end

            @. atmos.bands_cen = 0.5 * (atmos.bands_max + atmos.bands_min)
            @. atmos.bands_wid = abs(atmos.bands_max - atmos.bands_min)

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
            npd_cloud_component        =  4     # Number of components of clouds.
            npd_cloud_type             =  4     # Number of permitted types of clouds.
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

            SOCRATES.allocate_atm(  atmos.atm,   atmos.dimen, atmos.spectrum)
            SOCRATES.allocate_cld(  atmos.cld,   atmos.dimen, atmos.spectrum)
            SOCRATES.allocate_aer(  atmos.aer,   atmos.dimen, atmos.spectrum)
            SOCRATES.allocate_bound(atmos.bound, atmos.dimen, atmos.spectrum)

            # Fill with zeros - will be set inside of radtrans function at call time
            fill!(atmos.bound.flux_ground, 0.0)


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
            fill!(atmos.control.weight_band, 1.0)

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

            if atmos.overlap_method == "ro"
                # random overlap
                atmos.control.i_gas_overlap = SOCRATES.rad_pcf.ip_overlap_random

            elseif atmos.overlap_method == "ee"
                # equivalent extinction with correct scaling
                atmos.control.i_gas_overlap = SOCRATES.rad_pcf.ip_overlap_k_eqv_scl

            elseif atmos.overlap_method == "rorr"
                # random overlap with resorting and rebinning
                atmos.control.i_gas_overlap = SOCRATES.rad_pcf.ip_overlap_random_resort_rebin

            else
                @error "Invalid overlap method $(atmos.overlap_method)"
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

            # Warn user if all absorbers are condensable, which risks opacity going to zero
            if issubset(atmos.gas_soc_names, atmos.condensates)
                @warn "All absorbers are marked as condensable!"
                @warn "    Opacity will be zero if they all rainout of the atmosphere."
            end

        end # end socrates-only

        # VMRs are provided to SOCRATES when radtrans is called
        # For now, they are just stored inside the atmos struct

        # Metallicities provided?
        if !isempty(atmos.metal_orig)
            @info "Composition will be set by metallicity. Elem molar ratios:"
            for e in keys(atmos.metal_orig)
                @info @sprintf("    %-4s %.5f",e,atmos.metal_orig[e])
            end
        end

        # Print info on the gases, and check whether they are 'safe' or not
        gas_flags::String = ""
        g::String = ""
        @info "Allocating atmosphere with initial composition:"
        for i in 1:atmos.gas_num
            g = atmos.gas_names[i]
            gas_flags = ""
            atmos.gas_safe[g] = true

            # flag as not included in radtrans
            if (atmos.rt_scheme == RT_SOCRATES) && !(g in atmos.gas_soc_names)
                gas_flags *= "NO_OPACITY "
                atmos.gas_safe[g] = false
            end

            # flag as condensable
            if g in atmos.condensates
                gas_flags *= "COND "
                atmos.gas_safe[g] = false
            end

            # flag as containing stub thermo data
            if atmos.gas_dat[g].stub
                gas_flags *= "NO_THERMO "
                atmos.gas_safe[g] = false
            end

            # print info for this gas
            gas_flags *= String(Symbol(atmos.gas_dat[g].eos))*" "
            if !isempty(gas_flags)
                gas_flags = "($(gas_flags[1:end-1]))"
            end
            @info @sprintf("    %3d %-7s %6.2e %s", i, g, atmos.gas_vmr[g][end], gas_flags)
        end

        # There must be at least one 'safe' gas
        if !any(values(atmos.gas_safe)) && check_safe_gas
            @error "None of the supplied gases are considered 'safe'"
            @error "There must be at least one gas which satisfies criteria:"
            @error "    a) is dry, i.e. non-condensable"
            @error "    b) has opacity"
            @error "    c) has thermodynamic data"
            return false
        end

        # Calc layer properties using initial temperature profile.
        #    Can generate weird issues since the TOA temperature may be large
        #    large but pressure small, which gives it a low density. With the
        #    hydrostatic integrator, this can cause dr to blow up, especially
        #    with a low MMW gas. Should be okay as long as the T(p) provided
        #    by the user is more reasonable.
        calc_layer_props!(atmos)

        ################################
        # Aerosol processes
        #################################

        if atmos.rt_scheme == RT_SOCRATES

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
                atmos.cld.i_cloud_type[1]   = SOCRATES.rad_pcf.ip_cloud_type_water
                atmos.cld.i_condensed_param[1] = SOCRATES.rad_pcf.ip_drop_pade_2

                # reset parameters
                fill!(atmos.cld.condensed_param_list, 0.0)

                # input_cloud_cdf.f90, line 565
                n_cloud_parameter = 16 # for  ip_drop_pade_2
                for j in 1:atmos.nbands
                    for k in 1:n_cloud_parameter
                        atmos.cld.condensed_param_list[k, 1, j] = atmos.spectrum.Drop.parm_list[k, j,
                                                                        atmos.cld.i_condensed_param[1]]
                    end
                end

                # In-cloud fractions of different types of cloud
                fill!(atmos.cld.frac_cloud, 0.0)
                atmos.cld.frac_cloud[1,:,1] .= 1.0

            else
                atmos.cld.n_condensed  = 0
                atmos.cld.n_cloud_type = 0
            end
        end # end socrates only

        ###########################################
        # Surface properties
        ###########################################
        # allocate reflectance and emissivity arrays
        atmos.surf_r_arr = zeros(Float64, atmos.nbands) # spherical reflectivity
        atmos.surf_e_arr = ones(Float64, atmos.nbands)  # hemispherical emissivity

        if atmos.surface_material == "greybody"
            # grey albedo
            # Kirchoff's law: set emissivity equal to 1-albedo (spectrally)
            fill!(atmos.surf_r_arr, atmos.albedo_s)
            fill!(atmos.surf_e_arr, 1.0-atmos.albedo_s)

        else
            # spectral albedo and emissivity
            # Hapke2012: https://doi.org/10.1017/CBO9781139025683
            # Hammond2024: https://arxiv.org/abs/2409.04386

            # AGNI is flexible in accepting multiple forms of two-column surface data
            #   first column is WL in nm
            #   second column describes the optical properties 'value'
            # Options for second column:
            #    r = spherical reflectance (AKA Bond albedo), converted to e in code
            #    e = hemispherical emissivity, converted to r in code
            #    w = single scattering albedo, converted to r and e in code
            # Which of these 3 are provided must be indicated on 1st data line of file

            # reflectivity: r_hh = r_s = r0 * (1 - γ/(3+3γ) )
            # emissivity: e_h = 1- r_s
            # r0 is the diffusive reflectance = (1 − γ)/(1 + γ)
            # γ is the albedo factor = sqrt(1 − w)
            # w is the single-scattering albedo

            # try to find a matching file
            atmos.surface_material = abspath(atmos.surface_material)
            if !isfile(atmos.surface_material)
                @error "Could not find surface albedo file '$(atmos.surface_material)'"
                @error "    Get it by running: \$ ./src/get_data.sh surfaces"
                return false
            end

            # read spectral surface radiative properties from file
            (_srf_data::Array{Float64,2}, _srf_head) =
                        readdlm(atmos.surface_material, Float64; header=true, comments=true)

            _srf_w::Array{Float64, 1} = _srf_data[:,1]     # wavelength [nm]
            _srf_v::Array{Float64, 1} = _srf_data[:,2]     # value [dimensionless]

            # extrapolate to 0 wavelength with constant value
            pushfirst!(_srf_w, 0.0)
            pushfirst!(_srf_v, _srf_v[1])

            # extrapolate to large wavelength with constant value
            push!(_srf_w, 1e10)
            push!(_srf_v, _srf_v[end])

            # convert wl array from [nm] to [m]
            @. _srf_w = _srf_w / 1e9

            # sort data
            _srf_mask = sortperm(_srf_w)
            _srf_w = _srf_w[_srf_mask]
            _srf_v = _srf_v[_srf_mask]

            # create linear interpolator on the data
            _srf_i = extrapolate(interpolate((_srf_w,),_srf_v,Gridded(Linear())),Flat())

            # calculate r and e depending on what 'value' in the file actually represents
            if _srf_head[1] == "r"
                @debug "Reading surface data, spherical reflectance (r_hh, r_s)"
                for i in 1:atmos.nbands
                    atmos.surf_r_arr[i] = _srf_i(atmos.bands_cen[i])
                end
                @. atmos.surf_e_arr = 1.0 - atmos.surf_r_arr

            elseif _srf_head[1] == "e"
                @debug "Reading surface data, hemispherical emissivity (e_h)"
                for i in 1:atmos.nbands
                    atmos.surf_e_arr[i] = _srf_i(atmos.bands_cen[i])
                end
                @. atmos.surf_r_arr = 1.0 - atmos.surf_e_arr


            elseif _srf_head[1] == "w"
                @debug "Reading surface data, single scattering albedo (w)"
                γ::Float64 = 0.0
                for i in 1:atmos.nbands
                    γ = sqrt( 1 - _srf_i(atmos.bands_cen[i]) )
                    atmos.surf_r_arr[i] = (1−γ)/(1+γ) *(1 - γ/(3+3γ) )
                end
                @. atmos.surf_e_arr = 1.0 - atmos.surf_r_arr

            else
                @error "Unexpected format for surface data"
                @error "File path: $(atmos.surface_material)"
                @error "File header: $_srf_head"
                @error "Try downloading the data again using `src/get_data.sh`"
                return false
            end

            # set dummy value for the scalar variable containing the grey albedo
            atmos.albedo_s = Statistics.median(atmos.surf_r_arr)
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

        atmos.surf_flux =         zeros(Float64, atmos.nbands)

        atmos.contfunc_band =     zeros(Float64, (atmos.nlev_c,atmos.nbands))

        atmos.flux_sens =         0.0

        atmos.mask_l =            falses(atmos.nlev_l)      # Phase change
        atmos.mask_c =            falses(atmos.nlev_l)      # Dry convection

        atmos.phs_wrk_df =        zeros(Float64, atmos.nlev_c)  # flux difference
        atmos.phs_wrk_fl =        zeros(Float64, atmos.nlev_l)  # edge fluxes
        atmos.flux_l =            zeros(Float64, atmos.nlev_l)  # Latent heat / phase change
        atmos.flux_cdct =         zeros(Float64, atmos.nlev_l)  # Conduction

        atmos.flux_cdry =         zeros(Float64, atmos.nlev_l)  # Dry convection
        atmos.Kzz =               zeros(Float64, atmos.nlev_l)  # eddy diffusion coeff [m2 s-1]
        atmos.w_conv =            zeros(Float64, atmos.nlev_l)  # convective velocity [m s-1]
        atmos.λ_conv =            zeros(Float64, atmos.nlev_l)  # mixing length [m]

        atmos.flux_advect =       zeros(Float64, atmos.nlev_l)  # advective heat flux

        atmos.flux_tot =          zeros(Float64, atmos.nlev_l)
        atmos.flux_dif =          zeros(Float64, atmos.nlev_c)
        atmos.flux_deep =         zeros(Float64, atmos.nlev_l)  # Deep atmospheric heating flux
        atmos.heating_rate =      zeros(Float64, atmos.nlev_c)
        atmos.timescale_conv =    zeros(Float64, atmos.nlev_c)
        atmos.timescale_rad =     zeros(Float64, atmos.nlev_c)
        atmos.diagnostic_Ra =     zeros(Float64, atmos.nlev_c)

        # RFM values will be overwritten at runtime
        atmos.rfm_npts =          4
        atmos.rfm_wn =            zeros(Float64, atmos.rfm_npts)
        atmos.rfm_fl =            zeros(Float64, atmos.rfm_npts)

        # Mark as allocated
        atmos.is_alloc = true
        @debug "Allocate complete"

        return true
    end  # end of allocate

    """
    **Set atmosphere properties such that it is effectively transparent.**

    This will modify the surface pressure and disable gas opacity in SOCRATES.
    These changes cannot be reversed directly. To undo them, it is best to create and
    allocate a new atmosphere struct.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    """
    function make_transparent!(atmos::atmosphere.Atmos_t)

        # Turn off clouds
        fill!(atmos.cloud_arr_r, 0.0)
        fill!(atmos.cloud_arr_l, 0.0)
        fill!(atmos.cloud_arr_f, 0.0)

        # Set all gases to use ideal gas EOS
        # Avoid issues with partial pressures near zero
        for g in atmos.gas_names
            atmos.gas_dat[g].eos = phys.EOS_IDEAL
        end

        # Set surface pressure to be very small, but still larger than TOA pressure
        atmos.p_boa       = atmos.p_toa*1.10
        atmos.transspec_p = atmos.p_toa*1.05
        generate_pgrid!(atmos)

        # Set temperatures to be small, except the surface
        fill!(atmos.tmp[1:end-1],  atmos.tmp_floor)
        fill!(atmos.tmpl[1:end-1], atmos.tmp_floor)

        # Turn off gas opacity and rayleigh scattering
        atmos.control.l_continuum = false
        atmos.control.l_cont_gen  = false
        atmos.control.l_gas       = false
        atmos.control.l_rayleigh  = false
        atmos.control.l_cloud     = false

        # Flag as transparent
        atmos.transparent = true

        return nothing
    end

    """
    **Set cell-edge temperatures from cell-centre values.**

    Uses interpolation within the bulk of the column and extrapolation for the
    topmost edge.

    Arguments:
    - `atmos::Atmos_t`            the atmosphere struct instance to be used.
    """
    function set_tmpl_from_tmp!(atmos::atmosphere.Atmos_t)

        # Interpolate temperature to bulk cell-edge values (log-linear)
        itp = extrapolate( interpolate( (log10.(atmos.p[:]),),
                                        atmos.tmp,Gridded(Linear())
                                        ),Line())
        @. atmos.tmpl = itp(log10(atmos.pl))

        # Clamp
        clamp!(atmos.tmpl, atmos.tmp_floor, atmos.tmp_ceiling)

        return nothing
    end

    """
    **Get interleaved cell-centre and cell-edge PTR grid.**

    Arguments:
    - `atmos::Atmos_t`      the atmosphere struct instance to be used.

    Returns:
    - `arr_P::Array`        pressure array [Pa]
    - `arr_T::Array`        temperature array [K]
    - `arr_R::Array`        radius array [m]
    """
    function get_interleaved_ptr(atmos::atmosphere.Atmos_t)
        arr_R::Array{Float64, 1} = zeros(Float64, atmos.nlev_c + atmos.nlev_l)
        arr_T::Array{Float64, 1} = zeros(Float64, atmos.nlev_c + atmos.nlev_l)
        arr_P::Array{Float64, 1} = zeros(Float64, atmos.nlev_c + atmos.nlev_l)
        idx::Int = 1

        # top
        arr_R[1] = atmos.rl[1]
        arr_T[1] = atmos.tmpl[1]
        arr_P[1] = atmos.pl[1]

        # middle
        for i in 1:atmos.nlev_c
            idx = (i-1)*2
            arr_T[idx+1] = atmos.tmpl[i]
            arr_T[idx+2] = atmos.tmp[i]

            arr_P[idx+1] = atmos.pl[i]
            arr_P[idx+2] = atmos.p[i]

            arr_R[idx+1] = atmos.rl[i]
            arr_R[idx+2] = atmos.r[i]
        end

        # bottom
        arr_R[end] = atmos.rl[end]
        arr_T[end] = atmos.tmpl[end]
        arr_P[end] = atmos.pl[end]

        return arr_P, arr_T, arr_R
    end

    """
    **Estimate a diagnostic Rayleigh number in each layer.**

    Assuming that the Rayleigh number scales like `Ra ~ (wλ/κ)^(1/β)`
    Where `κ` is the thermal diffusivity and `β` is the convective beta parameter.

    This quantity must be taken lightly.

    Arguments:
    - `atmos::Atmos_t`      the atmosphere struct instance to be used.
    """
    function estimate_Ra!(atmos::atmosphere.Atmos_t)

        # Thermal diffusivity array
        κ::Array{Float64,1} = zero(atmos.layer_cp)
        @. κ = phys.calc_therm_diffus(atmos.layer_kc, atmos.layer_ρ, atmos.layer_cp)

        # One over beta
        ooβ::Float64 = 1.0 / phys.βRa

        # Estimate Rayleigh number
        @inbounds for i in 1:atmos.nlev_c
            atmos.diagnostic_Ra[i] = ( atmos.w_conv[i] * atmos.λ_conv[i] / κ[i]) ^ ooβ
        end

        return nothing
    end

    """
    **Estimate a diagnostic radiative timescale in each layer.**

    This quantity must be taken lightly.

    Arguments:
    - `atmos::Atmos_t`      the atmosphere struct instance to be used.
    """
    function estimate_timescale_rad!(atmos::atmosphere.Atmos_t)

        # Equation 10.1 from Seager textbook
        @inbounds for i in 1:atmos.nlev_c
            atmos.timescale_rad[i] = atmos.layer_cp[i] * (atmos.pl[i+1] - atmos.pl[i]) /
                                     (atmos.g[i] * 4 * phys.σSB * atmos.tmp[i]^3)
        end

        return nothing
    end

    """
    **Estimate a diagnostic convective timescale in each layer.**

    This quantity must be taken lightly.

    Arguments:
    - `atmos::Atmos_t`      the atmosphere struct instance to be used.
    """
    function estimate_timescale_conv!(atmos::atmosphere.Atmos_t)

        @inbounds for i in 1:atmos.nlev_c
            atmos.timescale_conv[i] = atmos.λ_conv[i] / max(atmos.w_conv[i], 1e-300)
        end

        return nothing
    end

end
