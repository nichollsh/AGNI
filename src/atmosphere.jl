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
    using LoopVectorization
    import Statistics
    import Interpolations: interpolate, Gridded, Linear, Flat, Line, extrapolate, Extrapolation
    import DelimitedFiles:readdlm

    # SOCRATES library
    const SOCRATESjl = abspath(ENV["RAD_DIR"], "julia","src","SOCRATES.jl")
    include(SOCRATESjl)

    # Local modules
    import ..phys
    import ..spectrum

    # Constants
    const AGNI_VERSION::String   = "1.7.3"
    const HYDROGRAV_STEPS::Int64 = 40

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

        # Radiation scheme performance
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
        p_boa::Float64          # Pressure at bottom [Pa]
        p_toa::Float64          # Pressure at top [Pa]
        tmp::Array{Float64,1}   # cc temperature [K]
        tmpl::Array{Float64,1}  # ce temperature [K]
        p::Array{Float64,1}     # cc pressure [Pa]
        pl::Array{Float64,1}    # ce pressure [Pa]
        r::Array{Float64,1}     # cc radius [m]
        rl::Array{Float64,1}    # ce radius [m]

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

        # Gas variables (incl gases which are not in spectralfile)
        gas_num::Int                                # Number of gases
        gas_names::Array{String,1}                  # List of gas names
        gas_vmr::Dict{String, Array{Float64,1}}     # Layer volume mixing ratios in dict, (key,value) = (gas_name,array)
        gas_sat::Dict{String, Array{Bool, 1}}       # Layer is saturated or cold-trapped
        gas_dat::Dict{String, phys.Gas_t}           # struct variables containing thermodynamic data
        cond_yield::Dict{String, Array{Float64,1}}  # condensate yield [kg/m^2] at each level (can be negative, representing evaporation)
        cond_surf::Dict{String, Float64}            # condensate accumulation left after evaporation (implicit surface liquid) [kg/m^2]
        gas_ovmr::Dict{String, Array{Float64,1}}    # original VMR values at model initialisation
        condensates::Array{String, 1}               # List of condensing gases (strings)
        condense_any::Bool                          # length(condensates)>0 ?

        # Ocean variables (surface liquid layering)
        ocean_calc::Bool                # INPUT: enable ocean calculations
        ocean_ob_frac::Float64          # INPUT: ocean basin area, as fraction of planet surface
        ocean_cs_height::Float64        # INPUT: continental shelf height [m]
        ocean_layers::Array{Tuple,1}    # OUTPUT: layer structure of surface liquids
        ocean_maxdepth::Float64         # OUTPUT: ocean depth at deepest part [m]
        ocean_areacov::Float64          # OUTPUT: fraction of planet surface covered by oceans
        ocean_topliq::String            # OUTPUT: name of top-most ocean component

        # Gases (only those in spectralfile)
        gas_soc_num::Int
        gas_soc_names::Array{String,1}

        # Layers' average properties
        real_gas::Bool                      # use real-gas equations of state where possible
        thermo_funct::Bool                  # use temperature-dependent evaluation of thermodynamic properties
        layer_ρ::Array{Float64,1}           # mass density [kg m-3]
        layer_μ::Array{Float64,1}           # mean molecular weight [kg mol-1]
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

        # Surface planck emission (incl. emissivity)
        surf_flux::Array{Float64, 1}

        # Contribution function (to outgoing flux) per-band
        contfunc_band::Array{Float64,2}     # LW only, not normalised

        # RFM line-by-line calculation
        rfm_fl::Array{Float64,1}            # upward flux [erg/(s cm2 cm-1)]
        rfm_wn::Array{Float64,1}            # wavenumber array [cm-1]
        rfm_npts::Int                       # number of points

        # Sensible heating
        C_d::Float64                        # Turbulent exchange coefficient [dimensionless]
        U::Float64                          # Wind speed [m s-1]
        flux_sens::Float64                  # Turbulent flux

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
        phs_tau_mix::Float64                # Time scale (mixed composition)
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
        ediv_add::Array{Float64, 1}     # Additional energy dissipation inside each cell [W m-3] (e.g. from advection)

        # Total energy flux
        flux_dif::Array{Float64,1}      # Flux lost at each level [W m-2]
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
        fastchem_maxiter::Int           # Maximum FC iterations
        fastchem_xtol::Float64          # FC solver tolerance
        fastchem_exec::String           # Path to fastchem executable
        fastchem_work::String           # Path to fastchem working directory

        # RFM radiative transfer
        flag_rfm::Bool                  # RFM enabled?
        rfm_exec::String                # Path to rfm executable
        rfm_work::String                # Path to rfm working directory
        rfm_parfile::String             # Path to rfm parfile. If empty, do not run RFM.

        # Observing properties
        transspec_p::Float64            # Pressure level probed in transmission [Pa]
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
    - `surface_material::String`        surface material (default is "greybody", but can point to file instead).
    - `albedo_s::Float64`               grey surface albedo used when `surface_material="greybody"`.
    - `tmp_floor::Float64`              temperature floor [K].
    - `C_d::Float64`                    turbulent heat exchange coefficient [dimensionless].
    - `U::Float64`                      surface wind speed [m s-1].
    - `Kzz_floor::Float64`              eddy diffusion coefficient, min value [cm2 s-1]
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
    - `real_gas::Bool`                  use real gas EOS where possible
    - `thermo_functions::Bool`          use temperature-dependent thermodynamic properties
    - `use_all_gases::Bool`             store information on all supported gases, incl those not provided in cfg
    - `fastchem_work::String`           working directory for fastchem
    - `fastchem_floor::Float64`         temperature floor on profile provided to fastchem
    - `fastchem_maxiter::Float64`       maximum solver iterations allowed by fastchem
    - `fastchem_xtol::Float64`          solution tolerance required of fastchem
    - `rfm_parfile::String`             path to LbL .par file provided to RFM

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
                    surface_material::String =  "greybody",
                    albedo_s::Float64 =         0.0,
                    tmp_floor::Float64 =        2.0,
                    C_d::Float64 =              0.001,
                    U::Float64 =                2.0,
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

                    real_gas::Bool =            true,
                    thermo_functions::Bool =    true,
                    use_all_gases::Bool =       false,
                    check_integrity::Bool =     true,

                    fastchem_work::String =     "",
                    fastchem_floor::Float64 =   273.0,
                    fastchem_maxiter::Float64 = 2e4,
                    fastchem_xtol::Float64 =    1.0e-4,

                    rfm_parfile::String =       "",

                    ocean_calc::Bool =          true,
                    ocean_ob_frac::Float64 =    0.6,
                    ocean_cs_height::Float64 =  3000.0
                    )::Bool

        if !isdir(OUT_DIR) && !isfile(OUT_DIR)
            mkdir(OUT_DIR)
        end

        @info  "Setting-up a new atmosphere struct"

        # Code versions
        atmos.SOCRATES_VERSION = readchomp(joinpath(ENV["RAD_DIR"],"version"))
        atmos.AGNI_VERSION = AGNI_VERSION
        @debug "AGNI VERSION = "*AGNI_VERSION
        @debug "Using SOCRATES at $(ENV["RAD_DIR"])"
        @debug "SOCRATES VERSION = "*atmos.SOCRATES_VERSION

        # Check SOCRATES version is valid
        SOCVER_minimum = 2407.2
        if parse(Float64, atmos.SOCRATES_VERSION) < SOCVER_minimum
            @error "SOCRATES is out of date and cannot be used!"
            @error "    found at $(ENV["RAD_DIR"])"
            @error "    version is "*atmos.SOCRATES_VERSION
        end

        atmos.benchmark   = false
        atmos.num_rt_eval = 0
        atmos.tim_rt_eval = 0.0

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
        atmos.FRAMES_DIR  =     joinpath(atmos.OUT_DIR, "frames")
        atmos.THERMO_DIR  =     joinpath(atmos.ROOT_DIR, "res", "thermodynamics")
        atmos.spectral_file =   abspath(spfile)
        atmos.all_channels =    all_channels
        atmos.overlap_method =  overlap_method

        atmos.real_gas      =   real_gas
        atmos.thermo_funct  =   thermo_functions

        atmos.tmp_floor =       max(0.1,tmp_floor)
        atmos.tmp_ceiling =     2.0e4

        if nlev_centre < 25
            nlev_centre = 25
            @warn "Adjusted number of levels to $nlev_centre"
        end
        atmos.nlev_c         =  nlev_centre
        atmos.nlev_l         =  atmos.nlev_c + 1
        atmos.tmp_surf =        max(tmp_surf, atmos.tmp_floor)
        atmos.grav_surf =       max(1.0e-7, gravity)
        atmos.zenith_degrees =  max(min(zenith_degrees,89.8), 0.2)
        atmos.surface_material= surface_material
        atmos.albedo_s =        max(min(albedo_s, 1.0 ), 0.0)
        atmos.instellation =    max(instellation, 0.0)
        atmos.albedo_b =        max(min(albedo_b,1.0), 0.0)
        atmos.s0_fact =         max(s0_fact,0.0)
        atmos.toa_heating =     atmos.instellation * (1.0 - atmos.albedo_b) *
                                    s0_fact * cosd(atmos.zenith_degrees)

        atmos.flux_int =        flux_int
        atmos.target_olr =      max(1.0e-10, target_olr)

        atmos.C_d =             max(0.0, C_d)
        atmos.U =               max(0.0, U)

        atmos.Kzz_floor =       max(0.0, Kzz_floor / 1e4)  # convert to SI units
        atmos.Kzz_ceiling =     1.0e20 / 1e4
        atmos.Kzz_pbreak =      1e5 # 1 bar as default location for break point
        atmos.Kzz_kbreak =      max(0.0, Kzz_floor)
        atmos.mlt_asymptotic =  mlt_asymptotic
        atmos.mlt_criterion =   mlt_criterion

        if atmos.real_gas && (atmos.mlt_criterion == 'l')
            @warn "Ledoux criterion not supported for real gases"
            @warn "    Switching criterion to Schwarzschild, neglecting MMW gradients"
            atmos.mlt_criterion = 's'
        end
        if !(atmos.mlt_criterion in ['s','l'])
            @error "Invalid choice for mlt_criterion: $(atmos.mlt_criterion)"
            @error "    Must be: 's' or 'l' only"
            return false
        end

        atmos.tmp_magma =       max(atmos.tmp_floor, tmp_magma)
        atmos.skin_d =          max(1.0e-9, skin_d)
        atmos.skin_k =          max(1.0e-9, skin_k)

        if p_top > p_surf
            @error "p_top must be less than p_surf"
            return false
        end

        atmos.p_toa =           p_top * 1.0e5 # Convert bar -> Pa
        atmos.p_boa =           p_surf * 1.0e5
        atmos.rp =              max(1.0, radius)

        # derived statistics
        atmos.interior_mass =   atmos.grav_surf * atmos.rp^2 / phys.G_grav
        atmos.interior_rho  =   3.0 * atmos.interior_mass / ( 4.0 * pi * atmos.rp^3)
        atmos.transspec_p   =   2e3     # 20 mbar = 2000 Pa

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

        # Initialise mesh geometry
        atmos.r             = zeros(Float64, atmos.nlev_c) # radii at cell centres [m]
        atmos.rl            = zeros(Float64, atmos.nlev_l) # radii at cell edges [m]
        atmos.layer_thick   = zeros(Float64, atmos.nlev_c) # geometric thickness [m]
        atmos.layer_mass    = zeros(Float64, atmos.nlev_c) # mass per unit area [kg m-2]
        atmos.layer_grav    = ones(Float64, atmos.nlev_c) * atmos.grav_surf

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
        atmos.phs_tau_mix = 1.0e6   # mixed composition case

        # Evaporation efficiency
        atmos.evap_efficiency = 0.05

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

        # The values will be stored in a dict of arrays
        atmos.gas_names =   Array{String}(undef, 0)           # list of names
        atmos.gas_dat =     Dict{String, phys.Gas_t}()        # dict of data structures
        atmos.gas_vmr  =    Dict{String, Array{Float64,1}}()  # dict of VMR arrays
        atmos.gas_ovmr  =   Dict{String, Array{Float64,1}}()  # ^ backup of initial values
        atmos.gas_sat  =    Dict{String, Array{Bool, 1}}()    # dict for saturation
        atmos.cond_yield =  Dict{String, Array{Float64,1}}()  # dict of condensate yield
        atmos.cond_surf =  Dict{String, Float64}()           # dict of ocean masses
        atmos.gas_num   =   0                                 # number of gases
        atmos.condensates   =   Array{String}(undef, 0)       # list of condensates

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
        end # end read VMR from dict

        # File input case
        if mf_source == 1
            # check file
            if !isfile(mf_path)
                @error "Could not read VMR file '$mf_path'"
                return false
            end
            @debug "Composition set by file"

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
                    push!(arr_p, atmos.p_boa*1.1)
                    push!(arr_x, arr_x[end])
                end

                # Set up interpolator using file data
                itp = interpolate((log10.(arr_p),),arr_x, Gridded(Linear()))

                # Set values in atmos struct
                for i in 1:atmos.nlev_c
                    atmos.gas_vmr[atmos.gas_names[gidx]][i] = itp(log10(atmos.p[i]))
                end
            end

        end # end read VMR from file

        # add extra gases if required
        if use_all_gases
            for gas in phys.gases_standard
                if !(gas in atmos.gas_names)
                    atmos.gas_vmr[gas] = zeros(Float64, atmos.nlev_c)
                    push!(atmos.gas_names, gas)
                    atmos.gas_num += 1
                end
            end
        end

        # backup mixing ratios from current state
        for k in keys(atmos.gas_vmr)
            atmos.gas_ovmr[k] = zeros(Float64, atmos.nlev_c)
            @turbo @. atmos.gas_ovmr[k] = atmos.gas_vmr[k]
        end

        # set condensation mask and yield values [kg]
        for g in atmos.gas_names
            atmos.gas_sat[g]    = falses(atmos.nlev_c)
            atmos.cond_yield[g] = zeros(Float64, atmos.nlev_c)
            atmos.cond_surf[g] = 0.0
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
        gas_fail = false
        @info "Loading thermodyamic data"
        for g in atmos.gas_names
            atmos.gas_dat[g] = phys.load_gas(atmos.THERMO_DIR, g,
                                                atmos.thermo_funct, atmos.real_gas;
                                                check_integrity=check_integrity)

            gas_fail = gas_fail || atmos.gas_dat[g].fail
        end
        if gas_fail
            @error "Problem when loading thermodynamic data"
            @error "Try downloading them again and/or updating AGNI."
            return false
        end

        # store condensates
        for c in condensates
            if !atmos.gas_dat[c].stub && !atmos.gas_dat[c].no_sat && !(c == "H2")
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

        # Must have at least one non-condensable gas
        if (length(condensates) == atmos.gas_num)
            @error "There must be at least one non-condensable gas"
            return false
        end

        # Ocean params
        atmos.ocean_calc  =     ocean_calc && atmos.condense_any
        atmos.ocean_ob_frac  =  max(0.0, min(1.0, ocean_ob_frac))
        atmos.ocean_cs_height = max(0.0, ocean_cs_height)
        atmos.ocean_maxdepth  = 0.0
        atmos.ocean_areacov   = 0.0
        atmos.ocean_topliq    = "_unset"
        atmos.ocean_layers    = Tuple[(1,"_unset",0.0,0.0),]

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

        # Fastchem
        atmos.flag_fastchem = false
        atmos.fastchem_work = joinpath(atmos.OUT_DIR, "fastchem/")  # default path
        if ("FC_DIR" in keys(ENV))

            @debug "FastChem env has been set"

            # check fastchem installation folder
            atmos.FC_DIR = abspath(ENV["FC_DIR"])
            if !isdir(atmos.FC_DIR)
                @error "Could not find fastchem folder at '$(atmos.FC_DIR)'"
                @error "Install FastChem with `\$ ./src/get_fastchem.sh`"
                return false
            end

            # check executable
            atmos.fastchem_exec = abspath(atmos.FC_DIR,"fastchem")
            atmos.flag_fastchem = isfile(atmos.fastchem_exec)
            if !atmos.flag_fastchem
                @error "Could not find fastchem executable inside '$(atmos.FC_DIR)'"
                @error "Install FastChem with `\$ ./src/get_fastchem.sh`"
                return false
            else
                @debug "Found FastChem executable"
            end

            # working directory for FC runtime files
            if !isempty(fastchem_work)
                atmos.fastchem_work = abspath(fastchem_work)
                @debug "Fastchem working dir set to $(atmos.fastchem_work)"
            else
                @debug "Fastchem working dir defaulting to $(atmos.fastchem_work)"
            end

            # make working directory
            rm(atmos.fastchem_work,force=true,recursive=true)
            mkdir(atmos.fastchem_work)
        else
            @debug "FastChem env variable not set"
        end
        # other parameters for FC
        atmos.fastchem_maxiter = fastchem_maxiter
        atmos.fastchem_floor   = fastchem_floor
        atmos.fastchem_xtol    = fastchem_xtol

        # RFM
        atmos.flag_rfm = !isempty(rfm_parfile)
        if atmos.flag_rfm
            atmos.rfm_parfile = abspath(rfm_parfile)
            @debug "RFM parfile set: $(atmos.rfm_parfile)"

            atmos.rfm_work = joinpath(atmos.OUT_DIR, "rfm/")
            rm(atmos.rfm_work,force=true,recursive=true)
            mkdir(atmos.rfm_work)
        end

        # Record that the parameters are set
        atmos.is_param = true
        atmos.is_solved = false
        atmos.is_converged = false

        @debug "Setup complete"
        return true
    end # end function setup


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
        atmos.transspec_r = atmos.r[idx]

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
    **Calculate properties within each layer of the atmosphere (e.g. density, mmw).**

    Arguments:
        - `atmos::Atmos_t`          the atmosphere struct instance to be used.
        - `ignore_errors::Bool`     do not generate errors from hydrostatic integrator.

    Returns:
        - `ok::Bool`                function result is ok
    """
    function calc_layer_props!(atmos::atmosphere.Atmos_t; ignore_errors::Bool=false)::Bool
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
        ok = ok && calc_profile_radius!(atmos, ignore_errors=ignore_errors)

        # Pass arrays to SOCRATES
        atmos.atm.p[1, :]           .= atmos.p[:]
        atmos.atm.p_level[1, 0:end] .= atmos.pl[:]
        atmos.atm.r_layer[1,:]      .= atmos.r[:]
        atmos.atm.r_level[1,0:end]  .= atmos.rl[:]
        atmos.atm.mass[1, :]        .= atmos.layer_mass[:]
        atmos.atm.density[1,:]      .= atmos.layer_ρ[:]

        return ok
    end

    """
    **Calculate radius and gravity for all layers.**

    Performs hydrostatic integration from the ground upwards.
    Requires density, temperature, pressure to have already been set.

    Arguments:
        - `atmos::Atmos_t`          the atmosphere struct instance to be used.
        - `ignore_errors::Bool`     do not generate errors from hydrostatic integrator.

    Returns:
        - `ok::Bool`                function result is ok
    """
    function calc_profile_radius!(atmos::atmosphere.Atmos_t;
                                    ignore_errors::Bool=false)::Bool

        # Reset arrays
        fill!(atmos.r         ,   atmos.rp)
        fill!(atmos.rl        ,   atmos.rp)
        fill!(atmos.layer_grav,   atmos.grav_surf)
        fill!(atmos.layer_thick,  1.0)
        fill!(atmos.layer_mass ,  1.0)

        # Temporary values
        grav::Float64       = atmos.grav_surf   # gravity at current level
        mass_encl::Float64  = atmos.interior_mass # mass enclosed within current level

        # Integrate from surface upwards
        for i in range(start=atmos.nlev_c, stop=1, step=-1)

            # Calculate gravity at lower edge
            grav = phys.G_grav * mass_encl / atmos.rl[i+1]^2

            # Integrate from lower edge to centre
            atmos.r[i] = integrate_hydrograv(atmos.rl[i+1], grav, atmos.pl[i+1], atmos.p[i], atmos.layer_ρ[i])

            # Calculate gravity at cell centre
            grav = phys.G_grav * mass_encl / atmos.r[i]^2

            # Integrate from centre to upper edge
            atmos.rl[i] = integrate_hydrograv(atmos.r[i], grav, atmos.p[i], atmos.pl[i], atmos.layer_ρ[i])

            # Store: Layer-centre gravity [m s-2]
            atmos.layer_grav[i] = grav

            # Store: Layer geometrical thickness [m]
            atmos.layer_thick[i] = atmos.rl[i] - atmos.rl[i+1]

            # Store: Layer-centre mass per unit area [kg m-2]
            atmos.layer_mass[i] = (atmos.pl[i+1] - atmos.pl[i])/atmos.layer_grav[i]

            # Add cumulative mass [kg]
            mass_encl += atmos.layer_mass[i] * 4 * pi * atmos.r[i]^2
        end

        return true
    end

    """
    **Integrate hydrostatic and gravity equations across a pressure interval.**

    Uses the classic fourth-order Runge-Kutta method.

    Arguments:
        - `r0::Float64`     radius   at start of interval [m]
        - `g0::Float64`     gravity  at start of interval [m s-2]
        - `p0::Float64`     pressure at start of interval [Pa]
        - `p1::Float64`     pressure at end   of interval [Pa]
        - `rho::Float64`    density throughout interval, constant [kg m-3]

    Returns:
        - `r1::Float64`     radius at end of interval [m]
    """
    function integrate_hydrograv(r0::Float64, g0::Float64,
                                    p0::Float64, p1::Float64, rho::Float64)::Float64

        # Gravity at given radius (neglecting mass within the interval)
        function _grav(r)
            return g0 * (r0/r)^2
        end

        # Derivative to integrate
        #   dr/dp = -1 / (rho * g(r))
        function _drdp(p,r)
            return -1 / (rho * _grav(r))
        end

        # Parameters
        dp::Float64  = (p1-p0)/HYDROGRAV_STEPS # this will be negative
        dp2::Float64 = dp/2
        k1::Float64  = 0.0; k2::Float64 = 0.0
        k3::Float64  = 0.0; k4::Float64 = 0.0

        # Integrate over pressure space
        pj::Float64 = p0    # rolling pressure (decreasing)
        rj::Float64 = r0    # rolling radius   (increasing)
        for _ in range(p0, stop=p1, step=dp)

            # gradient terms
            k1 = _drdp(pj,       rj)
            k2 = _drdp(pj + dp2, rj + k1*dp2)
            k3 = _drdp(pj + dp2, rj + k2*dp2)
            k4 = _drdp(pj + dp,  rj + k3*dp)

            # step height (increase)
            rj += dp/6 * (k1 + 2*k2 + 2*k3 + k4)

            # step pressure (decrease)
            pj += dp
        end

        return rj
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
            @turbo @. atmos.layer_μ += atmos.gas_vmr[gas] * atmos.gas_dat[gas].mmw
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

        # Shrink top-most layer to avoid doing too much extrapolation
        p_fact = 0.8
        atmos.p[1] = atmos.pl[1]*p_fact + atmos.p[1]*(1-p_fact)

        # Finally, convert arrays to actual pressure units [Pa]
        @turbo @. atmos.p  = 10.0 ^ atmos.p
        @turbo @. atmos.pl = 10.0 ^ atmos.pl

        return nothing
    end


    """
    **Allocate atmosphere arrays, prepare spectral files, and final steps.**

    Will not modify spectral file if `stellar_spectrum` is an empty string.
    Will treat star as blackbody with photospheric effective temperature of `stellar_Teff`
    if the parameter `stellar_spectrum` has value of `"blackbody"`.

    Arguments:
    - `atmos::Atmos_t`                 the atmosphere struct instance to be used.
    - `stellar_spectrum::String`       path to stellar spectrum csv file
    - `stellar_Teff::Float64`          star effective temperature if blackbody
    """
    function allocate!(atmos::atmosphere.Atmos_t, stellar_spectrum::String;
                        stellar_Teff::Float64=-1.0)::Bool

        @debug "Allocate atmosphere"
        if !atmos.is_param
            @error "Atmosphere struct has not been setup"
            return false
        end

        atmos.atm.n_profile = 1

        #########################################
        # spectral data
        #########################################

        # Validate files
        if !isfile(atmos.spectral_file)
            @error "Spectral file '$(atmos.spectral_file)' does not exist"
            @error "Try running `\$ ./src/get_data.sh`"
            @error "    e.g. to get CodenameXX you would run `\$ ./src/get_data.sh anyspec Codename XX`"
            return false
        end

        spectral_file_run::String  = joinpath([atmos.OUT_DIR, "runtime.sf"])
        spectral_file_runk::String = joinpath([atmos.OUT_DIR, "runtime.sf_k"])

        # Setup spectral file
        socstar::String = joinpath([atmos.OUT_DIR, "socstar.dat"])
        if !isempty(stellar_spectrum)
            @debug "Inserting stellar spectrum into spectral file"

            # Spectral file to be loaded, created in output folder
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
                    @error "Try running `\$ ./src/get_data.sh stellar`"
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

        # Allocate arrays, etc.
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

        # VMRs are provided to SOCRATES when radtrans is called
        # For now, they are just stored inside the atmos struct

        # Print info on the gases
        @info "Allocating atmosphere with composition:"
        gas_flags::String = ""
        g::String = ""
        for i in 1:atmos.gas_num
            g = atmos.gas_names[i]
            gas_flags = ""
            if !(g in atmos.gas_soc_names) # flag as not included in radtrans
                gas_flags *= "NO_OPACITY "
            end
            if g in atmos.condensates       # flag as condensable
                gas_flags *= "COND "
            end
            if atmos.gas_dat[g].stub        # flag as containing stub thermo data
                gas_flags *= "NO_THERMO "
            end
            gas_flags *= String(Symbol(atmos.gas_dat[g].eos))*" "
            if !isempty(gas_flags)
                gas_flags = "($(gas_flags[1:end-1]))"
            end
            @info @sprintf("    %3d %-7s %6.2e %s", i, g, atmos.gas_vmr[g][end], gas_flags)
        end


        # Calc layer properties using initial temperature profile.
        #    Can generate weird issues since the TOA temperature may be large
        #    large but pressure small, which gives it a low density. With the
        #    hydrostatic integrator, this can cause dr to blow up, especially
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

        atmos.control.i_angular_integration = SOCRATES.rad_pcf.ip_two_stream

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
                @error "Get these data with `\$ ./src/get_data.sh surfaces`"
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
            @turbo @. _srf_w = _srf_w / 1e9

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
                @turbo @. atmos.surf_e_arr = 1.0 - atmos.surf_r_arr

            elseif _srf_head[1] == "e"
                @debug "Reading surface data, hemispherical emissivity (e_h)"
                for i in 1:atmos.nbands
                    atmos.surf_e_arr[i] = _srf_i(atmos.bands_cen[i])
                end
                @turbo @. atmos.surf_r_arr = 1.0 - atmos.surf_e_arr


            elseif _srf_head[1] == "w"
                @debug "Reading surface data, single scattering albedo (w)"
                γ::Float64 = 0.0
                for i in 1:atmos.nbands
                    γ = sqrt( 1 - _srf_i(atmos.bands_cen[i]) )
                    atmos.surf_r_arr[i] = (1−γ)/(1+γ) *(1 - γ/(3+3γ) )
                end
                @turbo @. atmos.surf_e_arr = 1.0 - atmos.surf_r_arr

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

        atmos.flux_tot =          zeros(Float64, atmos.nlev_l)
        atmos.flux_dif =          zeros(Float64, atmos.nlev_c)
        atmos.ediv_add =          zeros(Float64, atmos.nlev_c)
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
        atmos.p_boa = atmos.p_toa * 1.1
        atmos.transspec_p = atmos.p_boa
        generate_pgrid!(atmos)

        # Set temperatures to be small, except the surface
        fill!(atmos.tmp[1:end-1],  atmos.tmp_floor)
        fill!(atmos.tmpl[1:end-1], atmos.tmp_floor)

        # Turn off gas opacity and rayleigh scattering
        atmos.control.l_continuum = false
        atmos.control.l_cont_gen  = false
        atmos.control.l_gas       = false
        atmos.control.l_rayleigh  = false

        # Turn off oceans
        atmos.ocean_maxdepth  = 0.0
        atmos.ocean_areacov   = 0.0
        atmos.ocean_topliq    = "_unset"
        atmos.ocean_layers    = Tuple[(1,"_unset",0.0,0.0),]

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
        @turbo @. κ = phys.calc_therm_diffus(atmos.layer_kc, atmos.layer_ρ, atmos.layer_cp)

        # One over beta
        ooβ::Float64 = 1.0 / phys.βRa

        # Estimate Rayleigh number
        @turbo @. atmos.diagnostic_Ra = ( atmos.w_conv * atmos.λ_conv / $κ) ^ $ooβ

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
                                     (atmos.layer_grav[i] * 4 * phys.σSB * atmos.tmp[i])
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

        @turbo @. atmos.timescale_conv = atmos.λ_conv / max(atmos.w_conv, 1e-300)

        return nothing
    end

end
