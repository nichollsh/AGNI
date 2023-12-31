# Contains the atmosphere module, which contains all of the core code 
# for setting-up the atmosphere and running SOCRATES. A lot of this code was 
# adapted from the SOCRATES-Julia examples.

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 


module atmosphere

    # Import libraries
    include("../socrates/julia/src/SOCRATES.jl")

    using Revise
    using Printf
    using NCDatasets
    using DataStructures
    using DelimitedFiles
    using PCHIPInterpolation
    using LinearAlgebra
    using Dates
    
    import moving_average
    import phys

    AGNI_VERSION     = "0.1"
    SOCRATES_VERSION = "2306"

    # Contains data pertaining to the atmosphere (fluxes, temperature, etc.)
    mutable struct Atmos_t

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
        dimen::SOCRATES.StrDim
        control::SOCRATES.StrCtrl
        spectrum::SOCRATES.StrSpecData
        atm::SOCRATES.StrAtm
        cld::SOCRATES.StrCld
        aer::SOCRATES.StrAer
        bound::SOCRATES.StrBound
        radout::SOCRATES.StrOut

        # Radiation parameters
        all_channels::Bool              # Use all bands?
        spectral_file::String           # Path to spectral file
        star_file::String               # Path to star spectrum
        albedo_b::Float64               # Enforced bond albedo (DO NOT USE ALONGSIDE CLOUD RADIATION EFFECTS)
        albedo_s::Float64               # Surface albedo
        zenith_degrees::Float64         # Solar zenith angle [deg]
        toa_heating::Float64            # Derived downward shortwave radiation flux at topmost level [W m-2]
        instellation::Float64           # Solar flux at top of atmopshere [W m-2]
        s0_fact::Float64                # Scale factor to instellation (cronin+14)
        tstar::Float64                  # Surface brightness temperature [K]
        grav_surf::Float64              # Surface gravity [m s-2]
        overlap_method::Int             # Absorber overlap method to be used

        # Band edges 
        nbands::Int
        bands_min::Array    # Lower wavelength [m]
        bands_max::Array    # Upper wavelength [m]

        # Pressure-temperature grid (with i=1 at the top of the model)
        nlev_c::Int         # Cell centre (count)
        nlev_l::Int         # Cell edge (count)

        p_boa::Float64      # Pressure at bottom [Pa]
        p_toa::Float64      # Pressure at top [Pa]
        rp::Float64         # planet radius [m]
        res_switching::Bool # resolution switching at high pressures

        tmp::Array          # cc temperature [K]
        tmpl::Array         # ce temperature
        p::Array            # cc pressure 
        pl::Array           # ce pressure
        z::Array            # cc height 
        zl::Array           # ce height 
        time::Array         # time [days]

        tmp_floor::Float64      # Temperature floor to prevent numerics [K]
        tmp_ceiling::Float64    # Temperature ceiling to prevent numerics [K]

        # Conductive skin
        skin_d::Float64         # skin thickness [m]
        skin_k::Float64         # skin thermal conductivity [W m-1 K-1] (You can find reasonable values here: https://doi.org/10.1016/S1474-7065(03)00069-X)
        tmp_magma::Float64      # Mantle temperature [K]

        # Mole fractions (= VMR)
        gases::Array            # List of gas names 
        input_x::Dict           # Layer mole fractions in dict format, incl gases not in spfile (key,value) = (gas_name,array)
        layer_x::Array          # Layer mole fractions in matrix format, excl gases not in spfile [lvl, gas_idx]

        # Layers' average properties
        thermo_funct::Bool      # use temperature-dependent evaluation
        layer_density::Array    # density [kg m-3]
        layer_mmw::Array        # mean molecular weight [kg mol-1]
        layer_cp::Array         # heat capacity at const-p [J K-1 kg-1]
        layer_grav::Array       # gravity [m s-2]
        layer_mass::Array       # mass per unit area [kg m-2]

        # Calculated bolometric radiative fluxes (W m-2)
        flux_d_lw::Array  # down component, lw 
        flux_u_lw::Array  # up component, lw
        flux_n_lw::Array  # net upward, lw

        flux_d_sw::Array  # down component, sw 
        flux_u_sw::Array  # up component, sw
        flux_n_sw::Array  # net upward, sw

        flux_d::Array    # down component, lw+sw 
        flux_u::Array    # up component, lw+sw 
        flux_n::Array    # net upward, lw+sw 

        # Calculated per-band radiative fluxes (W m-2)
        band_d_lw::Array  # down component, lw 
        band_u_lw::Array  # up component, lw
        band_n_lw::Array  # net upward, lw

        band_d_sw::Array  # down component, sw 
        band_u_sw::Array  # up component, sw
        band_n_sw::Array  # net upward, sw

        # Contribution function (to outgoing flux) per-band
        contfunc_norm::Array    # LW+SW, and normalised by maximum value at each wavelength

        # Sensible heating
        C_d::Float64        # Turbulent exchange coefficient [dimensionless]
        U::Float64          # Wind speed [m s-1]
        flux_sens::Float64  # Turbulent flux

        # Convection 
        mask_c::Array       # Layers which are (or were recently) convective (value is set to >0)
        mask_decay::Int     # How long is 'recent'
        flux_c::Array       # Dry convective fluxes from MLT
        Kzz::Array          # Eddy diffusion coefficient from MLT

        # Cloud and condensation
        mask_p::Array       # Layers which are (or were recently) condensing liquid
        re::Array           # Effective radius of the droplets [m] (drizzle forms above 20 microns)
        lwm::Array          # Liquid water mass fraction [kg/kg] - how much liquid vs. gas is there upon cloud formation? 0 : saturated water vapor does not turn liquid ; 1 : the entire mass of the cell contributes to the cloud
        clfr::Array         # Water cloud fraction - how much of the current cell turns into cloud? 0 : clear sky cell ; 1 : the cloud takes over the entire area of the cell (just leave at 1 for 1D runs)

        # Total energy flux
        flux_tot::Array     # Total upward-directed flux at cell edges

        # Heating rate 
        heating_rate::Array # radiative heating rate [K/day]

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
    - `tstar::Float64`                  effective surface temperature to provide upward longwave flux at the bottom of the atmosphere [K].
    - `gravity::Float64`                gravitational acceleration at the surface [m s-2].
    - `radius::Float64`                 planet radius at the surface [m].
    - `nlev_centre::Int`                number of model levels.
    - `p_surf::Float64`                 total surface pressure [bar].
    - `p_top::Float64`                  total top of atmosphere pressure [bar].
    - `mf_dict=nothing`                 dictionary of mole fractions in the format (key,value)=(gas,mf).
    - `mf_path=nothing`                 path to file containing mole fractions at each level.
    - `albedo_s::Float64=0.0`           surface albedo.
    - `tmp_floor::Float64=50.0`         temperature floor [K].
    - `C_d::Float64=0.001`              turbulent heat exchange coefficient [dimensionless].
    - `U::Float64=10.0`                 surface wind speed [m s-1].
    - `tmp_magma::Float64=3000.0`       mantle temperature [K].
    - `skin_d::Float64=0.05`            skin thickness [m].
    - `skin_k::Float64=2.0`             skin thermal conductivity [W m-1 K-1].
    - `all_channels::Bool=true`         use all channels available for RT?
    - `overlap_method::Int=4`           SOCRATES gaseous overlap scheme (2: rand overlap, 4: equiv extinct, 8: ro+resort+rebin).
    - `flag_rayleigh::Bool=false`       include rayleigh scattering?
    - `flag_gcontinuum::Bool=false`     include generalised continuum absorption?
    - `flag_continuum::Bool=false`      include continuum absorption?
    - `flag_aerosol::Bool=false`        include aersols?
    - `flag_cloud::Bool=false`          include clouds?
    - `res_switching::Bool=false`       use resolution switching at high pressures?
    - `thermo_functions::Bool=true`     use temperature-dependent thermodynamic properties
    """
    function setup!(atmos::atmosphere.Atmos_t, 
                    ROOT_DIR::String, OUT_DIR::String, 
                    spfile::String, 
                    instellation::Float64, s0_fact::Float64, albedo_b::Float64, zenith_degrees::Float64,
                    tstar::Float64,
                    gravity::Float64, radius::Float64,
                    nlev_centre::Int, p_surf::Float64, p_top::Float64;
                    mf_dict=                    nothing,
                    mf_path =                   nothing,
                    albedo_s::Float64 =         0.0,
                    tmp_floor::Float64 =        80.0,
                    C_d::Float64 =              0.001,
                    U::Float64 =                2.0,
                    tmp_magma::Float64 =        3000.0,
                    skin_d::Float64 =           0.05,
                    skin_k::Float64 =           2.0,
                    overlap_method::Int =       4,
                    all_channels::Bool  =       true,
                    flag_rayleigh::Bool =       false,
                    flag_gcontinuum::Bool =     false,
                    flag_continuum::Bool =      false,
                    flag_aerosol::Bool =        false,
                    flag_cloud::Bool =          false,
                    res_switching::Bool =       false,
                    thermo_functions::Bool =    true
                    )

        if !isdir(OUT_DIR) && !isfile(OUT_DIR)
            mkdir(OUT_DIR)
        end

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
        atmos.overlap_method =  Int(overlap_method)
       
        if res_switching && (p_surf < 5.0)
            println("WARNING: Resolution switching is enabled, but surface pressure is quite low")
        end 
        atmos.res_switching = res_switching
        atmos.thermo_funct  = thermo_functions

        atmos.tmp_floor =       max(0.1,tmp_floor)
        atmos.tmp_ceiling =     5000.0

        atmos.nlev_c         =  max(nlev_centre,45)
        atmos.nlev_l         =  atmos.nlev_c + 1
        atmos.tstar =           max(tstar, atmos.tmp_floor)
        atmos.grav_surf =       max(1.0e-4, gravity)
        atmos.zenith_degrees =  max(min(zenith_degrees,89.5), 0.5)
        atmos.albedo_s =        max(min(albedo_s, 1.0 ), 0.0)
        atmos.instellation =    max(instellation, 0.0)
        atmos.albedo_b =        max(min(albedo_b,1.0), 0.0)
        atmos.s0_fact =         max(s0_fact,0.0)
        atmos.toa_heating =     atmos.instellation * (1.0 - atmos.albedo_b) * s0_fact * cosd(atmos.zenith_degrees)
        
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

        # store contribution function?
        # atmos.control.l_contrib_func_band = false

        # Initialise temperature grid to be isothermal
        atmos.tmpl = ones(Float64, atmos.nlev_l) .* atmos.tstar
        atmos.tmp =  ones(Float64, atmos.nlev_c) .* atmos.tstar

        # Initialise pressure grid with current p_toa and p_boa
        generate_pgrid!(atmos, switch=atmos.res_switching)
        atmos.z          = zeros(Float64, atmos.nlev_c)
        atmos.zl         = zeros(Float64, atmos.nlev_l)
        atmos.layer_grav = ones(Float64, atmos.nlev_c) * atmos.grav_surf

        # Initialise cloud properties
        atmos.re         = zeros(Float64, atmos.nlev_c) 
        atmos.lwm        = zeros(Float64, atmos.nlev_c)
        atmos.clfr       = zeros(Float64, atmos.nlev_c) 

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
        
        # The values will be stored in a dict of arrays regardless, because 
        # we do not yet know which order the gases should be indexed in.

        atmos.input_x = Dict()  # dict of arrays 
        
        # Dict input case
        if mf_source == 0
            for (key, value) in mf_dict  # store as arrays
                gas_valid = strip(key, ' ')
                if key in SOCRATES.input_head_pcf.header_gas
                    atmos.input_x[gas_valid] = ones(Float64, atmos.nlev_c) * value
                end
            end
        end

        # File input case 
        if mf_source == 1
            # check file 
            if !isfile(mf_path)
                error("Could not read file '$mf_path'")
            end

            # get header
            mf_head = readlines(abspath(mf_path))[1]
            mf_head = mf_head[2:end]  # remove comment symbol at start
            mf_head = replace(mf_head, " " => "")  # remove whitespace
            heads = split(mf_head, ",") # split by columm

            # create arrays 
            for h in heads
                gas_valid = strip(h, ' ')
                if gas_valid in SOCRATES.input_head_pcf.header_gas
                    atmos.input_x[gas_valid] = zeros(Float64, atmos.nlev_c)
                end 
            end 

            # get body
            mf_body = readdlm(abspath(mf_path), ',', Float64; header=false, skipstart=2)
            mf_body = transpose(mf_body)

            # set composition by interpolating with pressure array 
            for li in 4:lastindex(heads)

                gas = heads[li]

                if !(gas in keys(atmos.input_x)) # skip tem, pre, hgt
                    continue
                end

                arr_p = mf_body[1,:]
                arr_x = mf_body[li,:]

                # Extend loaded profile to lower pressures (prevent domain error)
                if arr_p[1] > atmos.p_toa
                    pushfirst!(arr_p,   atmos.p_toa/1.1)
                    pushfirst!(arr_x,   arr_x[1] )
                end

                # Extend loaded profile to higher pressures 
                if arr_p[end] < atmos.p_boa
                    push!(arr_p,   atmos.p_boa*1.1)
                    push!(arr_x,   arr_x[end])
                end
                
                # Set up interpolator
                itp = Interpolator(arr_p, arr_x)
                
                # Set values 
                for i in 1:atmos.nlev_c
                    atmos.input_x[gas][i] = itp(atmos.p[i])
                end 
            end 

        end

        # Check that we actually stored some values
        if length(keys(atmos.input_x)) == 0
            error("No mole fractions were stored")
        end

        # Record that the parameters are set
        atmos.is_param = true
        atmos.is_solved = false 
        atmos.is_converged = false
        return nothing
    end

    """
    **Get the mole fraction of a gas within the atmosphere.**

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    - `gas::String`             name of the gas (e.g. "H2O").
    - `lvl::Int`                model level to measure mole fraction (-1 => use input_x value)

    Returns:
    - `x::Float64`              mole fraction of `gas`.
    """
    function get_x(atmos::atmosphere.Atmos_t, gas::String, lvl::Int)::Float64

        gas_valid = strip(gas, ' ')

        x = 0.0
        if gas_valid in keys(atmos.input_x)
            if lvl > 0
                i_gas = findfirst(==(gas), atmos.gases)
                x = atmos.layer_x[lvl,i_gas]
            else 
                x = atmos.input_x[gas_valid][1]
            end
        end 

        return x
    end

    """
    **Calculate properties within each layer of the atmosphere (e.g. density, mmw).**

    Assumes that the atmosphere may be treated as an ideal gas.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    """
    function calc_layer_props!(atmos::atmosphere.Atmos_t)
        if !atmos.is_param
            error(" atmosphere parameters have not been set")
        end

        # Set pressure arrays in SOCRATES 
        atmos.atm.p[1, :] .= atmos.p[:]
        atmos.atm.p_level[1, 0:end] .= atmos.pl[:]

        # mmw, cp, rho
        atmos.layer_mmw     = zeros(Float64, atmos.nlev_c)
        atmos.layer_density = zeros(Float64, atmos.nlev_c)
        atmos.layer_cp      = zeros(Float64, atmos.nlev_c)
        atmos.layer_mass    = zeros(Float64, atmos.nlev_c)
        for i in 1:atmos.atm.n_layer

            # for each gas
            for (i_gas, gas) in enumerate(atmos.gases)

                # set mmw
                atmos.layer_mmw[i] += atmos.layer_x[i,i_gas] * phys.lookup_safe("mmw",gas)

                # set cp
                t = -1.0
                if atmos.thermo_funct 
                    t = atmos.tmp[i]
                end
                atmos.layer_cp[i] += atmos.layer_x[i,i_gas] * phys.lookup_safe("cp",gas,tmp=t)
            end

            # density (assumes ideal gas)
            atmos.layer_density[i] = ( atmos.p[i] * atmos.layer_mmw[i] )  / (phys.R_gas * atmos.tmp[i]) 
            atmos.atm.density[1,i] = atmos.layer_density[i]
        end

        # geometrical height and gravity
        # dz = -dp / (rho * g)
        atmos.z[:] .= 0.0
        atmos.zl[:] .= 0.0
        atmos.layer_grav[:] .= 0.0
        for i in range(atmos.nlev_c, 1, step=-1)

            # Technically, g and z should be integrated as coupled equations,
            # but here they are not. This loose integration has been found to 
            # be reasonable in all of my tests.

            # Integrate from lower edge to centre
            g1 = atmos.grav_surf * (atmos.rp^2.0) / ((atmos.rp + atmos.zl[i+1])^2.0)
            dzc = phys.R_gas * atmos.tmp[i] / (atmos.layer_mmw[i] * g1 * 0.5 * (atmos.pl[i+1] + atmos.p[i])) * (atmos.pl[i+1] - atmos.p[i]) 
            atmos.z[i] = atmos.zl[i+1] + dzc

            # Integrate from centre to upper edge
            g2 = atmos.grav_surf * (atmos.rp^2.0) / ((atmos.rp + atmos.z[i])^2.0)
            dzl = phys.R_gas * atmos.tmp[i] / (atmos.layer_mmw[i] * g2 * 0.5 * (atmos.p[i] + atmos.pl[i] )) * (atmos.p[i]- atmos.pl[i]) 
            atmos.zl[i] = atmos.z[i] + dzl

            atmos.layer_grav[i] = 0.5 * (g1 + g2)

            if (dzl < 1e-20) || (dzc < 1e-20)
                error("Height integration resulted in dz <= 0")
            end
        end 

        # Mass
        for i in 1:atmos.atm.n_layer
            atmos.layer_mass[i] = (atmos.atm.p_level[1, i] - atmos.atm.p_level[1, i-1])/atmos.layer_grav[i]
            atmos.atm.mass[1, i] = atmos.layer_mass[i]
        end

        return nothing
    end 

    """
    **Generate pressure grid.**

    Equally log-spaced between p_boa and p_boa in the nominal configuration.
    
    If switching is enabled, the grid resolution will change at 10% of the 
    surface pressure, with at least switch_f percent of the levels being placed 
    at pressures greater than this value. This enhances the resolution near the 
    surface, which may help with model stability.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    - `switch::Float64`         switch resolution at 0.1*p_surf level
    - `switch_f::Int`           lower region level proportion minimum [%].
    """
    function generate_pgrid!(atmos::atmosphere.Atmos_t; switch::Bool=false, switch_f::Int=25)

        # Try canonical distribution
        atmos.pl = 10 .^ range( log10(atmos.p_toa), stop=log10(atmos.p_boa), length=atmos.nlev_l)

        # Check if resolution switch is required
        if switch 

            # Upper and lower part counts
            switch_p  = 0.10 * atmos.p_boa
            switch_lo = ceil(Int, float(switch_f) / 100 * atmos.nlev_l)
            switch_up = atmos.nlev_l - switch_lo

            # Check requirement
            count_lo = 0
            for i in 1:atmos.nlev_l 
                if atmos.pl[i] > switch_p 
                    count_lo += 1
                end
            end

            # Upper and lower parts
            if count_lo < switch_lo 
                pl_up = 10 .^ range( log10(atmos.p_toa), stop=log10(switch_p),    length=switch_up + 1)
                pl_lo = 10 .^ range( log10(switch_p),    stop=log10(atmos.p_boa), length=switch_lo)
                pl_lo[1] = sqrt(pl_up[end-1]*pl_lo[2])
                atmos.pl = vcat(pl_up[1:end-1], pl_lo[1:end])
            end
        end
    
        # Set bottom and top layers to be quite small  (to avoid giving the extrapolation too much control)
        # atmos.pl[end-1] = max(atmos.pl[end-1], atmos.pl[end] * 0.999)
        # atmos.pl[2]     = min(atmos.pl[2]    , atmos.pl[1]   / 0.999)

        # Set pressure cell-centre array using geometric mean
        atmos.p =    zeros(Float64, atmos.nlev_c)
        for i in 1:atmos.nlev_c
            atmos.p[i] = sqrt(atmos.pl[i]*atmos.pl[i+1])
        end

        return nothing
    end

    
    """
    **Allocate atmosphere arrays, prepare spectral files, and final steps before RT calculations.**

    Arguments:
    - `atmos::Atmos_t`                 the atmosphere struct instance to be used.
    - `spfile_has_star::Bool=false`    does the spectral file already contain the stellar spectrum?
    - `stellar_spectrum::String=""`    path to stellar spectrum csv file (used when `spfile_has_star==false`)
    - `spfile_noremove::Bool=false`    should the runtime spectral file not be removed after use?
    """
    function allocate!(atmos::atmosphere.Atmos_t;
                        spfile_has_star::Bool=false,
                        stellar_spectrum::String="",
                        spfile_noremove::Bool=false
                        )

        if !atmos.is_param
            error("atmosphere parameters have not been set")
        end

        atmos.time = zeros(Float64, atmos.nlev_c)
        
        atmos.atm.n_profile = 0

        #########################################
        # spectral data
        #########################################

        # Validate files
        if !isfile(atmos.spectral_file)
            error("Spectral file '$(atmos.spectral_file)' does not exist")
        end
        if !spfile_has_star && (stellar_spectrum == "")
            error("No stellar spectrum provided")
        end
        if spfile_has_star && (stellar_spectrum != "")
            error("Stellar spectrum provided, but spectral file already contains this block")
        end 
        if !isfile(stellar_spectrum) && !spfile_has_star
            error("Stellar spectrum file '$(stellar_spectrum)' does not exist")
        end

        # Record location
        if !spfile_has_star
            atmos.star_file = abspath(stellar_spectrum)
        else 
            atmos.star_file = "IN_SPECTRAL_FILE"
        end

        # Spectral file to be loaded
        spectral_file_run  = joinpath([atmos.OUT_DIR, ".spfile_runtime"])  

        # Insert stellar spectrum 
        if !spfile_has_star
            println("Python: inserting stellar spectrum")
            runfile = joinpath([atmos.ROOT_DIR, "src", "insert_stellar.py"])
            run(`python $runfile $(stellar_spectrum) $(atmos.spectral_file) $(spectral_file_run)`)
        elseif atmos.spectral_file != spectral_file_run
            cp(atmos.spectral_file,      spectral_file_run,      force=true)
            cp(atmos.spectral_file*"_k", spectral_file_run*"_k", force=true)
        end

        # Insert rayleigh scattering (optionally)
        if atmos.control.l_rayleigh
            println("Python: inserting rayleigh scattering")
            co2_x = get_x(atmos, "CO2", -1) 
            n2_x  = get_x(atmos, "N2" , -1)
            h2o_x = get_x(atmos, "H2O", -1)
            runfile = joinpath([atmos.ROOT_DIR, "src", "insert_rayleigh.py"])
            run(`python $runfile $spectral_file_run $co2_x $n2_x $h2o_x`)
        end

        # Validate files
        if occursin("NaN", readchomp(spectral_file_run*"_k"))
            error("Spectral _k file contains NaN values")
        end 
        if occursin("NaN", readchomp(spectral_file_run))
            error("Spectral file contains NaN values")
        end 

        # Read-in spectral file to be used at runtime
        atmos.control.spectral_file = spectral_file_run
        SOCRATES.set_spectrum(spectrum=atmos.spectrum, 
                              spectral_file=atmos.control.spectral_file, 
                              l_all_gasses=true)

        if !spfile_noremove                              
            rm(spectral_file_run; force=true)
            rm(spectral_file_run*"_k"; force=true)
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
        npd_region                 =  3 # Number of regions in a layer

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

        # Set composition
        atmos.gases = Array{String}(undef, atmos.spectrum.Gas.n_absorb)
        atmos.layer_x = zeros(Float64, (atmos.nlev_c,length(atmos.gases)))

        for i_gas in 1:atmos.spectrum.Gas.n_absorb
            ti = atmos.spectrum.Gas.type_absorb[i_gas]
            gas = SOCRATES.input_head_pcf.header_gas[ti]

            atmos.gases[i_gas] = gas

            if gas in keys(atmos.input_x)
                atmos.layer_x[:,i_gas] .= atmos.input_x[gas]
            else
                atmos.layer_x[:,i_gas] .= 0.0
            end

            # Values are provided to SOCRATES when radtrans is called
        end

        # Renormalise mole fractions
        for i in 1:atmos.nlev_c
            tot = 0.0
            for i_gas in 1:atmos.spectrum.Gas.n_absorb 
                tot += atmos.layer_x[i,i_gas]
            end 

            if tot == 0
                error("Layer $i has all-zero mole fractions - check your spectral file")
            end

            for i_gas in 1:atmos.spectrum.Gas.n_absorb
                ti = atmos.spectrum.Gas.type_absorb[i_gas]
                gas = SOCRATES.input_head_pcf.header_gas[ti]

                atmos.layer_x[i,i_gas] /= tot

                if gas in keys(atmos.input_x)
                    atmos.input_x[gas][i] = atmos.layer_x[i,i_gas]
                end

            end 
        end

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
            atmos.cld.n_condensed = 0
        end

        SOCRATES.allocate_cld(  atmos.cld,   atmos.dimen, atmos.spectrum)
        SOCRATES.allocate_cld_prsc(atmos.cld, atmos.dimen, atmos.spectrum)

        if atmos.control.l_cloud
            atmos.cld.n_condensed       = 1
            atmos.cld.type_condensed[1] = SOCRATES.rad_pcf.ip_clcmp_st_water
            atmos.cld.n_cloud_type      = 1
            atmos.cld.i_cloud_type[1]   = SOCRATES.rad_pcf.ip_cloud_type_homogen
            atmos.cld.i_condensed_param[1] = 5
        else
            atmos.cld.n_condensed = 0
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

        atmos.mask_p =            zeros(Float64, atmos.nlev_c)
        atmos.mask_c =            zeros(Float64, atmos.nlev_c)
        atmos.flux_c =            zeros(Float64, atmos.nlev_l)
        atmos.Kzz =               zeros(Float64, atmos.nlev_l)

        atmos.flux_tot =          zeros(Float64, atmos.nlev_l)

        atmos.heating_rate =      zeros(Float64, atmos.nlev_c)

        # Mark as allocated
        atmos.is_alloc = true

        return nothing
    end  # end of allocate

    function radtrans!(atmos::atmosphere.Atmos_t, lw::Bool; calc_cf=false)
        if !atmos.is_alloc
            error("atmosphere arrays have not been allocated")
        end
        if !atmos.is_param
            error("atmosphere parameters have not been set")
        end

        # Longwave or shortwave calculation?
        if lw
            atmos.control.isolir = SOCRATES.rad_pcf.ip_infra_red
        else
            atmos.control.isolir = SOCRATES.rad_pcf.ip_solar
        end

        # Check files are acceptable and set instellation if doing SW flux
        if lw
            if !Bool(atmos.spectrum.Basic.l_present[6])
                error("The spectral file contains no data for the Planckian function.\nCheck that the file contains a stellar spectrum or run AGNI with the `stf` argument.")
            end

            if Bool(atmos.spectrum.Basic.l_present[2])
                atmos.control.l_solar_tail_flux = true
            end

        else
            if !Bool(atmos.spectrum.Basic.l_present[2])
                error("The spectral file contains no solar spectral data.")
            end
             
            atmos.bound.zen_0[1] = 1.0/cosd(atmos.zenith_degrees)   #   Convert the zenith angles to secants.
            atmos.bound.solar_irrad[1] = atmos.toa_heating / cosd(atmos.zenith_degrees)
        end

        atmos.bound.rho_alb[:, SOCRATES.rad_pcf.ip_surf_alb_diff, :] .= atmos.albedo_s

        # Set the two-stream approximation to be used (-t f)
        if lw
            atmos.control.i_2stream = 12 # Practical improved flux method (1985) with Elsasser's diffusivity (D=1.66)
        else
            atmos.control.i_2stream = 16 # Practical improved fl ux method (original form of 1980)
        end

        #####################################
        # Angular integration
        # see src/aux/angular_control_cdf.f
        #####################################

        # Cl_run_cdf +R flag
        atmos.control.l_rescale = true
        if atmos.control.l_rescale
            atmos.control.l_henyey_greenstein_pf = true
        end

        # The internal SOCRATES solver used for the two-stream calculations (-v flag)
        if atmos.control.l_cloud
            # 16 is recommended for cloudy-sky
            # 17 is recommended for cloud with separate stratiform and convective regions
            atmos.control.i_solver = 16 
        else 
            # 13 is recommended for clear-sky
            atmos.control.i_solver = 13  
        end 
        

        #      Arrays of fluxes must be of the full size.
        atmos.dimen.nd_2sg_profile =        atmos.dimen.nd_profile
        atmos.dimen.nd_flux_profile =       atmos.dimen.nd_profile
        atmos.dimen.nd_radiance_profile =   1
        atmos.dimen.nd_j_profile =          1
        atmos.dimen.nd_viewing_level =      1
        atmos.dimen.nd_sph_coeff =          1

        # Reset dimen.nd_max_order to reduce memory requirements
        atmos.dimen.nd_max_order = 1

        #####################################
        # surface properties
        # see src/aux/assign_surface_char_cdf.f
        # IP_surface_char  = 51, file suffix 'surf'
        #####################################

        if atmos.control.i_angular_integration == SOCRATES.rad_pcf.ip_two_stream
            if !lw
                atmos.bound.rho_alb[:, SOCRATES.rad_pcf.ip_surf_alb_dir, :] .= atmos.bound.rho_alb[:, SOCRATES.rad_pcf.ip_surf_alb_diff, :]
            end
        end

        ###################################################
        # Cloud information
        ###################################################

        if atmos.control.l_cloud
            atmos.cld.w_cloud[1,:]               .= atmos.clfr[:]
            atmos.cld.condensed_mix_ratio[1,:,1] .= atmos.lwm[:]
            atmos.cld.condensed_dim_char[1,:,1]  .= atmos.re[:]
        end

        ###################################################
        # Treatment of scattering
        ###################################################

        atmos.control.i_scatter_method = SOCRATES.rad_pcf.ip_scatter_full
        for i in atmos.control.first_band:atmos.control.last_band
            atmos.control.i_scatter_method_band[i] = atmos.control.i_scatter_method
        end

        ####################################################
        # Temperature
        ###################################################

        clamp!(atmos.tmp,  atmos.tmp_floor, atmos.tmp_ceiling)
        clamp!(atmos.tmpl, atmos.tmp_floor, atmos.tmp_ceiling)

        atmos.atm.p[1, :] .= atmos.p[:]
        atmos.atm.t[1, :] .= atmos.tmp[:]
        atmos.atm.p_level[1, 0:end] .= atmos.pl[:]
        atmos.atm.t_level[1, 0:end] .= atmos.tmpl[:]

        # calc_layer_props!(atmos)

        if lw
            atmos.bound.t_ground[1] = atmos.tstar
        end

        if lw
            atmos.control.l_ir_source_quad = true
        end

        ######################################################
        # Run radiative transfer model
        ######################################################   
        
        atmos.control.l_contrib_func_band = calc_cf

        for i_gas in 1:length(atmos.gases)
            for i in 1:atmos.nlev_c 
                # convert mole fraction to mass mixing ratio
                atmos.atm.gas_mix_ratio[1, i, i_gas] = atmos.layer_x[i,i_gas] * phys.lookup_safe("mmw", atmos.gases[i_gas]) / atmos.layer_mmw[i]
            end
        end 

        SOCRATES.radiance_calc(atmos.control, atmos.dimen, atmos.spectrum, 
                                atmos.atm, atmos.cld, atmos.aer, 
                                atmos.bound, atmos.radout)

        # Store new fluxes in atmos struct
        if lw 
            # LW case
            fill!(atmos.flux_u_lw, 0.0)
            fill!(atmos.flux_d_lw, 0.0)
            fill!(atmos.flux_n_lw, 0.0)
            fill!(atmos.band_u_lw, 0.0)
            fill!(atmos.band_d_lw, 0.0)
            fill!(atmos.band_n_lw, 0.0)
            for lv in 1:atmos.nlev_l      # sum over levels
                for ba in 1:atmos.nbands  # sum over bands
                    idx = lv+(ba-1)*atmos.nlev_l

                    atmos.band_d_lw[lv,ba] = atmos.radout.flux_down[idx]
                    atmos.band_u_lw[lv,ba] = atmos.radout.flux_up[idx]
                    atmos.band_n_lw[lv,ba] = atmos.band_u_lw[lv,ba] - atmos.band_d_lw[lv,ba]

                    atmos.flux_d_lw[lv] += max(0.0, atmos.radout.flux_down[idx])
                    atmos.flux_u_lw[lv] += max(0.0, atmos.radout.flux_up[idx])
                end 
                atmos.flux_n_lw[lv] = atmos.flux_u_lw[lv] - atmos.flux_d_lw[lv] 
            end

            # Normalised contribution function (only LW stream contributes)
            fill!(atmos.contfunc_norm,0.0)
            if calc_cf
                for ba in 1:atmos.dimen.nd_channel
                    for lv in 1:atmos.nlev_c               
                        atmos.contfunc_norm[lv,ba] = atmos.radout.contrib_funcf_band[1,lv,ba]/maximum(atmos.radout.contrib_funcf_band[1,:,ba])
                    end 
                end
            end
            atmos.is_out_lw = true 
        else
            # SW case
            fill!(atmos.flux_u_sw, 0.0)
            fill!(atmos.flux_d_sw, 0.0)
            fill!(atmos.flux_n_sw, 0.0)
            fill!(atmos.band_u_sw, 0.0)
            fill!(atmos.band_d_sw, 0.0)
            fill!(atmos.band_n_sw, 0.0)
            for lv in 1:atmos.nlev_l                # sum over levels
                for ba in 1:atmos.dimen.nd_channel  # sum over bands
                    idx = lv+(ba-1)*atmos.nlev_l

                    atmos.band_d_sw[lv,ba] = atmos.radout.flux_down[idx]
                    atmos.band_u_sw[lv,ba] = atmos.radout.flux_up[idx]
                    atmos.band_n_sw[lv,ba] = atmos.band_u_sw[lv,ba] - atmos.band_d_sw[lv,ba]

                    atmos.flux_d_sw[lv] += max(0.0, atmos.radout.flux_down[idx])
                    atmos.flux_u_sw[lv] += max(0.0, atmos.radout.flux_up[idx])
                end 
                atmos.flux_n_sw[lv] = atmos.flux_u_sw[lv] - atmos.flux_d_sw[lv]
            end
            atmos.is_out_sw = true
        end

        # Store net fluxes when we have both SW and LW components
        if (atmos.is_out_lw && atmos.is_out_sw)
            atmos.flux_d = atmos.flux_d_lw .+ atmos.flux_d_sw
            atmos.flux_u = atmos.flux_u_lw .+ atmos.flux_u_sw
            atmos.flux_n = atmos.flux_n_lw .+ atmos.flux_n_sw
        end

        return nothing
    end # end of radtrans

    # Calculate sensible heat flux (turbulence at surface boundary)
    function sensible!(atmos::atmosphere.Atmos_t)
        # TKE scheme for this 1D case
        # transports energy from the surface to the bottom node
        atmos.flux_sens = atmos.layer_cp[end]*atmos.layer_mmw[end]*atmos.p[end]/(phys.R_gas*atmos.tmp[end]) * atmos.C_d * atmos.U * (atmos.tstar-atmos.tmp[end])
        return nothing
    end


    """
    **Calculate dry convective fluxes using mixing length theory.**

    Uses the mixing length formulation outlined by Joyce & Tayar (2023), which 
    was also implemented in Lee et al. (2023), and partially outlined in an 
    earlier paper by Robinson & Marley (2014).

    Convective energy transport fluxes are calculated at every level edge, just 
    like the radiative fluxes. This is not compatible with moist convection. By 
    using MLT to parameterise convection, we can also calculate Kzz directly.

    The mixing length is set to asymptotically approach H (for z>>H) or z (for 
    z<H) as per Blackadar (1962). Or alternatively it can be set equal to H.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `pmin::Float64=0.0`               pressure below which convection is disabled.
    - `mltype::Int=1`                   mixing length (0: fixed, 1: asymptotic)
    """
    function mlt!(atmos; pmin::Float64=0.0, mltype::Int=1)

        # Reset arrays
        atmos.flux_c[:] .= 0.0
        atmos.Kzz[:] .= 0.0

        # Work variables 
        H = 0.0; l = 0.0; w = 0.0
        m1 = 0.0; m2 = 0.0; mt = 0.0
        grav = 0.0; mu = 0.0; c_p = 0.0; rho = 0.0
        grad_ad = 0.0; grad_pr = 0.0

        # Loop from bottom upwards
        for i in range(start=atmos.nlev_l-1 , step=-1, stop=3) 

            # Optionally skip low pressures 
            if atmos.p[i] < pmin
                continue
            end

            m1 = atmos.layer_mass[i-1]
            m2 = atmos.layer_mass[i]
            mt = m1+m2

            grav = (atmos.layer_grav[i] * m2 + atmos.layer_grav[i-1] * m1)/mt
            mu   = (atmos.layer_mmw[i]  * m2 + atmos.layer_mmw[i-1]  * m1)/mt
            c_p  = (atmos.layer_cp[i]   * m2 + atmos.layer_cp[i-1]   * m1)/mt

            grad_ad = (phys.R_gas / mu) / c_p
            grad_pr = ( log(atmos.tmp[i-1]/atmos.tmp[i]) ) / ( log(atmos.p[i-1]/atmos.p[i]) )

            # Check instability
            if (grad_pr > grad_ad)

                # Check if this layer is condensing (this shouldn't ever be
                # true, because the condensation curve dT/dp is too shallow)
                if atmos.mask_p[i] > 0
                    println("    WARNING: Somehow unstable to dry convection in a condensing region!")
                end 

                rho = (atmos.layer_density[i] * m2 + atmos.layer_density[i-1] * m1)/mt

                atmos.mask_c[i]   = atmos.mask_decay
                atmos.mask_c[i-1] = atmos.mask_decay
                
                # Pressure scale height
                H = phys.R_gas * atmos.tmpl[i] / (mu * grav)

                # Mixing length
                if mltype == 0
                    # Fixed
                    l = H
                elseif mltype == 1
                    # Asymptotic 
                    l = phys.k_vk * atmos.zl[i] / (1 + phys.k_vk * atmos.zl[i]/H)
                else 
                    # Otherwise
                    error("Invalid mixing length type selected")
                end

                # Characteristic velocity (from Brunt-Vasalla frequency of parcel oscillations)
                w = l * sqrt(grav/H * (grad_pr-grad_ad))

                # Convective flux
                atmos.flux_c[i] = 0.5 * rho * c_p * w * atmos.tmpl[i] * l/H * (grad_pr-grad_ad)

                # Thermal eddy diffusion coefficient
                atmos.Kzz[i] = w * l
            
            end
        end
        return nothing
    end # end of mlt

    """
    **Calculate heating rates at cell-centres from the total flux.**

    Requires the total flux to have already been set (atmos.flux_tot). Heating 
    rates are calculated in units of kelvin per day.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    """
    function calc_hrates!(atmos::atmosphere.Atmos_t)

        dF::Float64 = 0.0
        dp::Float64 = 0.0

        atmos.heating_rate[:] .= 0.0

        for i in 1:atmos.nlev_c
            dF = atmos.flux_tot[i+1] - atmos.flux_tot[i]
            dp = atmos.pl[i+1] - atmos.pl[i]
            atmos.heating_rate[i] = (atmos.layer_grav[i] / atmos.layer_cp[i]) * dF/dp # K/s
        end

        atmos.heating_rate *= 86400.0 # K/day

        return nothing
    end

     # Dry convective adjustment (returning the temperature tendency without modifying the atmos struct)
     function adjust_dry(atmos::atmosphere.Atmos_t, nsteps::Int)

        tmp_old = zeros(Float64, atmos.nlev_c)  # old temperatures
        tmp_tnd = zeros(Float64, atmos.nlev_c)  # temperature tendency

        tmp_old[:] .+= atmos.tmp[:]

        for i in 1:nsteps

            # Downward pass
            for i in 2:atmos.nlev_c

                T1 = atmos.tmp[i-1]    # upper layer
                p1 = atmos.p[i-1]

                T2 = atmos.tmp[i]  # lower layer
                p2 = atmos.p[i]
                
                cp = 0.5 * ( atmos.layer_cp[i-1] * atmos.layer_mmw[i-1] +  atmos.layer_cp[i] * atmos.layer_mmw[i])
                pfact = (p1/p2)^(phys.R_gas / cp)
                
                # If slope dT/dp is steeper than adiabat (unstable), adjust to adiabat
                if T1 < T2*pfact
                    Tbar = 0.5 * ( T1 + T2 )
                    T2 = 2.0 * Tbar / (1.0 + pfact)
                    T1 = T2 * pfact
                    atmos.tmp[i-1] = T1
                    atmos.tmp[i]   = T2
                end
            end

            # Upward pass
            for i in atmos.nlev_c:2

                T1 = atmos.tmp[i-1]
                p1 = atmos.p[i-1]

                T2 = atmos.tmp[i]
                p2 = atmos.p[i]
                
                cp = 0.5 * ( atmos.layer_cp[i-1] * atmos.layer_mmw[i-1] + atmos.layer_cp[i] * atmos.layer_mmw[i])
                pfact = (p1/p2)^(phys.R_gas / cp)

                if T1 < T2*pfact
                    Tbar = 0.5 * ( T1 + T2 )
                    T2 = 2.0 * Tbar / ( 1.0 + pfact)
                    T1 = T2 * pfact
                    atmos.tmp[i-1] = T1
                    atmos.tmp[i]   = T2 
                end 
            end
        end
        
        # Calculate tendency
        tmp_tnd[:] .= atmos.tmp[:] .- tmp_old[:]

        # Restore temperature array
        atmos.tmp[:] .= tmp_old[:]

        return tmp_tnd
    end

    # Apply condensation according to vapour-liquid coexistance curve (return mask of condensing levels)
    function apply_vlcc!(atmos::atmosphere.Atmos_t, gas::String)

        changed = falses(atmos.nlev_c)

        i_gas = findfirst(==(gas), atmos.gases)

        # Check if each level is condensing. If it is, place on phase curve
        for i in 1:atmos.nlev_c

            x = atmos.layer_x[i,i_gas]
            if x < 1.0e-10 
                continue
            end

            # Cell centre only
            Tsat = phys.calc_Tdew(gas,atmos.p[i] * x )
            if atmos.tmp[i] < Tsat
                atmos.tmp[i] = Tsat
                changed[i] = true 
                
                atmos.re[i]   = 1.0e-5  # 10 micron droplets
                atmos.lwm[i]  = 0.8     # 80% of the saturated vapor turns into cloud
                atmos.clfr[i] = 1.0     # The cloud takes over the entire cell
            else 
                atmos.re[i]   = 0.0
                atmos.lwm[i]  = 0.0
                atmos.clfr[i] = 0.0
            end
        end
        
        return changed
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
    topmost edge. Does not set the bottommost edge.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `limit_change::Bool=false`        limit change relative to old value
    - `back_interp::Bool=false`         interpolate resultant tmpl back to tmp (should be avoided)
    """
    function set_tmpl_from_tmp!(atmos::atmosphere.Atmos_t; limit_change::Bool=false, back_interp::Bool=false)

        top_old_e = atmos.tmpl[1]

        # Interpolate temperature to bulk cell-edge values (log-linear)
        itp = Interpolator(log.(atmos.p), atmos.tmp)
        atmos.tmpl[2:end-1] .= itp.(log.(atmos.pl[2:end-1]))

        # Extrapolate top edge temperature (log-linear)
        grad_dt = atmos.tmp[1] - atmos.tmp[2]
        grad_dp = log(atmos.p[1]/atmos.p[2])
        atmos.tmpl[1] = atmos.tmp[1] + grad_dt/grad_dp * log(atmos.pl[1]/atmos.p[1])
        if limit_change
            atmos.tmpl[1] = 0.9 * atmos.tmpl[1] + 0.1 * top_old_e
        end

        # Second interpolation back to cell-centres.
        # This can help prevent grid-imprinting issues, but in some cases it can 
        # lead to unphysical behaviour, because it gives the interpolator too 
        # much control over the temperature profile at small scales.
        if back_interp
            clamp!(atmos.tmpl, atmos.tmp_floor, atmos.tmp_ceiling)
            itp = Interpolator(log.(atmos.pl), atmos.tmpl)  
            atmos.tmp[:] .= itp.(log.(atmos.p[:]))
        end

        # Limit domain
        clamp!(atmos.tmpl, atmos.tmp_floor, atmos.tmp_ceiling)

        return nothing
    end 

    # Get interleaved cell-centre and cell-edge PT grid
    function get_interleaved_pt(atmos::atmosphere.Atmos_t)
        arr_T = zeros(Float64, atmos.nlev_c + atmos.nlev_l)
        arr_P = zeros(Float64, atmos.nlev_c + atmos.nlev_l)

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

    # Write current interleaved pressure/temperature grid to a CSV file
    function write_pt(atmos::atmosphere.Atmos_t, fname::String)

        arr_P, arr_T = atmosphere.get_interleaved_pt(atmos)

        # Remove old file if exists
        rm(fname, force=true)

        open(fname, "w") do f
            write(f, "# pressure  , temperature \n")
            write(f, "# [Pa]      , [K] \n")
            for i in 1:atmos.nlev_l+atmos.nlev_c
                @printf(f, "%1.5e , %1.5e \n", arr_P[i], arr_T[i])
            end
        end
        return nothing
    end

    # Write current cell-edge fluxes to a CSV file
    function write_fluxes(atmos::atmosphere.Atmos_t, fname::String)

        # Remove old file if exists
        rm(fname, force=true)

        open(fname, "w") do f
            write(f, "# pressure  , U_LW        , D_LW        , N_LW        , U_SW        , D_SW        , N_SW        , U           , D           , N           , C           , tot      \n")
            write(f, "# [Pa]      , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]  \n")
            for i in 1:atmos.nlev_l
                @printf(f, "%1.5e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e, %+1.4e \n", 
                          atmos.pl[i], 
                          atmos.flux_u_lw[i], atmos.flux_d_lw[i], atmos.flux_n_lw[i],
                          atmos.flux_u_sw[i], atmos.flux_d_sw[i], atmos.flux_n_sw[i],
                          atmos.flux_u[i],    atmos.flux_d[i],    atmos.flux_n[i],
                          atmos.flux_c[i],    atmos.flux_tot[i]
                          )
            end
        end
        return nothing
    end

    # Write atmosphere data to a NetCDF file
    function write_ncdf(atmos::atmosphere.Atmos_t, fname::String)

        # See the tutorial at:
        # https://github.com/Alexander-Barth/NCDatasets.jl#create-a-netcdf-file

        # Note that the content of the NetCDF file is designed to be compatible
        # with what AEOLUS writes. As a result, they can both be integrated 
        # into PROTEUS without compatibility issues.

        # Create dataset handle
        fname = abspath(fname)
        rm(fname, force=true)
        ds = Dataset(fname,"c")

        # Global attributes
        ds.attrib["description"]        = "AGNI atmosphere data"
        ds.attrib["date"]               = Dates.format(now(), "yyyy-u-dd HH:MM:SS")
        ds.attrib["hostname"]           = gethostname()
        ds.attrib["username"]           = ENV["USER"]
        ds.attrib["AGNI_version"]       = AGNI_VERSION
        ds.attrib["SOCRATES_version"]   = SOCRATES_VERSION 

        plat = "Generic"
        if Sys.isapple()
            plat = "Darwin"
        elseif Sys.iswindows()
            plat = "Windows"
        elseif Sys.islinux()
            plat = "Linux"
        end
        ds.attrib["platform"] = plat
        
        # ----------------------
        # Create dimensions
        nlev_c = Int(atmos.nlev_c)
        nlev_l = nlev_c + 1
        ngases = length(atmos.gases)
        nchars = 16

        defDim(ds, "nlev_c", nlev_c)        # Cell centres
        defDim(ds, "nlev_l", nlev_l)        # Cell edges
        defDim(ds, "ngases", ngases)        # Gases
        defDim(ds, "nchars", nchars)        # Length of string containing gas names
        defDim(ds, "nbands", atmos.nbands)  # Number of spectral bands

        # ----------------------
        # Scalar quantities  
        #    Create variables
        var_tstar =     defVar(ds, "tstar",         Float64, (), attrib = OrderedDict("units" => "K"))      # BOA LW BC
        var_inst =      defVar(ds, "instellation",  Float64, (), attrib = OrderedDict("units" => "W m-2"))  # Solar flux at TOA
        var_s0fact =    defVar(ds, "inst_factor",   Float64, ())                                            # Scale factor applied to instellation
        var_albbond =   defVar(ds, "bond_albedo",   Float64, ())                                            # Bond albedo used to scale-down instellation
        var_toah =      defVar(ds, "toa_heating",   Float64, (), attrib = OrderedDict("units" => "W m-2"))  # TOA SW BC
        var_tmagma =    defVar(ds, "tmagma",        Float64, (), attrib = OrderedDict("units" => "K"))      # Magma temperature
        var_tmin =      defVar(ds, "tfloor",        Float64, (), attrib = OrderedDict("units" => "K"))      # Minimum temperature
        var_tmax =      defVar(ds, "tceiling",      Float64, (), attrib = OrderedDict("units" => "K"))      # Maximum temperature
        var_plrad =     defVar(ds, "planet_radius", Float64, (), attrib = OrderedDict("units" => "m"))      # Value taken for planet radius
        var_gsurf =     defVar(ds, "surf_gravity",  Float64, (), attrib = OrderedDict("units" => "m s-2"))  # Surface gravity
        var_albsurf =   defVar(ds, "surf_albedo",   Float64, ())                                            # Surface albedo
        var_fray =      defVar(ds, "flag_rayleigh", Char, ())                                               # Includes rayleigh scattering?
        var_fcon =      defVar(ds, "flag_continuum",Char, ())                                               # Includes continuum absorption?
        var_fcld =      defVar(ds, "flag_cloud"    ,Char, ())                                               # Includes clouds?
        var_tfun =      defVar(ds, "thermo_funct"  ,Char, ())                                               # Using temperature-dependent thermodynamic functions
        var_ssol =      defVar(ds, "solved"        ,Char, ())                                               # Has a solver been used?
        var_scon =      defVar(ds, "converged"     ,Char, ())                                               # Did the solver converge?
        var_znth =      defVar(ds, "zenith_angle"  ,Float64, (), attrib = OrderedDict("units" => "deg"))    # Zenith angle of direct stellar radiation
        var_sknd =      defVar(ds, "cond_skin_d"   ,Float64, (), attrib = OrderedDict("units" => "m"))      # Conductive skin thickness
        var_sknk =      defVar(ds, "cond_skin_k"   ,Float64, (), attrib = OrderedDict("units" => "W m-1 K-1"))    # Conductive skin thermal conductivity
        var_specfile =  defVar(ds, "specfile"      ,String, ())     # Path to spectral file when read
        var_starfile =  defVar(ds, "starfile"      ,String, ())     # Path to star file when read

        #     Store data
        var_tstar[1] =      atmos.tstar 
        var_inst[1] =       atmos.instellation
        var_s0fact[1] =     atmos.s0_fact
        var_albbond[1] =    atmos.albedo_b
        var_znth[1] =       atmos.zenith_degrees
        var_toah[1] =       atmos.toa_heating
        var_tmagma[1] =     atmos.tmp_magma
        var_tmin[1] =       atmos.tmp_floor
        var_tmax[1] =       atmos.tmp_ceiling
        var_plrad[1]  =     atmos.rp
        var_gsurf[1] =      atmos.grav_surf
        var_albsurf[1] =    atmos.albedo_s

        if atmos.control.l_rayleigh
            var_fray[1] = 'y'
        else
            var_fray[1] = 'n'
        end 

        if atmos.control.l_cont_gen
            var_fcon[1] = 'y'
        else
            var_fcon[1] = 'n'
        end 

        if atmos.control.l_cloud
            var_fcld[1] = 'y'
        else
            var_fcld[1] = 'n'
        end 

        if atmos.thermo_funct
            var_tfun[1] = 'y'
        else
            var_tfun[1] = 'n'
        end 

        if atmos.is_solved
            var_ssol[1] = 'y'
        else 
            var_ssol[1] = 'n'
        end

        if atmos.is_converged
            var_scon[1] = 'y'
        else 
            var_scon[1] = 'n'
        end

        var_sknd[1] = atmos.skin_d
        var_sknk[1] = atmos.skin_k

        var_specfile[1] = atmos.spectral_file
        var_starfile[1] = atmos.star_file

        # ----------------------
        # Vector quantities
        #    Create variables
        var_p =         defVar(ds, "p",         Float64, ("nlev_c",), attrib = OrderedDict("units" => "Pa"))
        var_pl =        defVar(ds, "pl",        Float64, ("nlev_l",), attrib = OrderedDict("units" => "Pa"))
        var_tmp =       defVar(ds, "tmp",       Float64, ("nlev_c",), attrib = OrderedDict("units" => "K"))
        var_tmpl =      defVar(ds, "tmpl",      Float64, ("nlev_l",), attrib = OrderedDict("units" => "K"))
        var_z =         defVar(ds, "z",         Float64, ("nlev_c",), attrib = OrderedDict("units" => "m"))
        var_zl =        defVar(ds, "zl",        Float64, ("nlev_l",), attrib = OrderedDict("units" => "m"))
        var_grav =      defVar(ds, "gravity",   Float64, ("nlev_c",), attrib = OrderedDict("units" => "m s-2"))
        var_mmw =       defVar(ds, "mmw",       Float64, ("nlev_c",), attrib = OrderedDict("units" => "kg mol-1"))
        var_gases =     defVar(ds, "gases",     Char,    ("nchars", "ngases")) # Transposed cf AEOLUS because of how Julia stores arrays
        var_x =         defVar(ds, "x_gas",     Float64, ("ngases", "nlev_c")) # ^^
        var_fdl =       defVar(ds, "fl_D_LW",   Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
        var_ful =       defVar(ds, "fl_U_LW",   Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
        var_fnl =       defVar(ds, "fl_N_LW",   Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
        var_fds =       defVar(ds, "fl_D_SW",   Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
        var_fus =       defVar(ds, "fl_U_SW",   Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
        var_fns =       defVar(ds, "fl_N_SW",   Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
        var_fd =        defVar(ds, "fl_D",      Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
        var_fu =        defVar(ds, "fl_U",      Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
        var_fn =        defVar(ds, "fl_N",      Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
        var_fc =        defVar(ds, "fl_C",      Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
        var_ft =        defVar(ds, "fl_tot",    Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
        var_hr =        defVar(ds, "hrate",     Float64, ("nlev_c",), attrib = OrderedDict("units" => "K day-1"))
        var_kzz =       defVar(ds, "Kzz",       Float64, ("nlev_l",), attrib = OrderedDict("units" => "m2 s-1"))
        var_bmin =      defVar(ds, "bandmin",   Float64, ("nbands",), attrib = OrderedDict("units" => "m"))
        var_bmax =      defVar(ds, "bandmax",   Float64, ("nbands",), attrib = OrderedDict("units" => "m"))
        var_bdl =       defVar(ds, "ba_D_LW",   Float64, ("nbands","nlev_l"), attrib = OrderedDict("units" => "W m-2"))
        var_bul =       defVar(ds, "ba_U_LW",   Float64, ("nbands","nlev_l"), attrib = OrderedDict("units" => "W m-2"))
        var_bnl =       defVar(ds, "ba_N_LW",   Float64, ("nbands","nlev_l"), attrib = OrderedDict("units" => "W m-2"))
        var_bds =       defVar(ds, "ba_D_SW",   Float64, ("nbands","nlev_l"), attrib = OrderedDict("units" => "W m-2"))
        var_bus =       defVar(ds, "ba_U_SW",   Float64, ("nbands","nlev_l"), attrib = OrderedDict("units" => "W m-2"))
        var_bns =       defVar(ds, "ba_N_SW",   Float64, ("nbands","nlev_l"), attrib = OrderedDict("units" => "W m-2"))
        var_cfn =       defVar(ds, "contfunc",  Float64, ("nbands","nlev_c"))

        #     Store data
        var_p[:] =      atmos.p
        var_pl[:] =     atmos.pl
        var_tmp[:] =    atmos.tmp
        var_tmpl[:] =   atmos.tmpl
        var_z[:]    =   atmos.z
        var_zl[:]   =   atmos.zl
        var_mmw[:]  =   atmos.layer_mmw
        var_grav[:]  =  atmos.layer_grav

        for i_gas in 1:ngases 
            for i_char in 1:nchars 
                var_gases[i_char, i_gas] = ' '
            end 
            for i_char in 1:length(atmos.gases[i_gas])
                var_gases[i_char,i_gas] = atmos.gases[i_gas][i_char]
            end 
        end 
        
        for i_gas in 1:ngases 
            for i_lvl in 1:nlev_c
                var_x[i_gas, i_lvl] = atmos.layer_x[i_lvl, i_gas]
            end 
        end 

        var_fdl[:] =    atmos.flux_d_lw
        var_ful[:] =    atmos.flux_u_lw
        var_fnl[:] =    atmos.flux_n_lw

        var_fds[:] =    atmos.flux_d_sw
        var_fus[:] =    atmos.flux_u_sw
        var_fns[:] =    atmos.flux_n_sw
        
        var_fd[:] =     atmos.flux_d
        var_fu[:] =     atmos.flux_u
        var_fn[:] =     atmos.flux_n

        var_fc[:] =     atmos.flux_c

        var_ft[:] =     atmos.flux_tot
        var_hr[:] =     atmos.heating_rate

        var_kzz[:] =    atmos.Kzz

        var_bmin[:] =   atmos.bands_min
        var_bmax[:] =   atmos.bands_max

        for lv in 1:atmos.nlev_l 
            for ba in 1:atmos.nbands
                var_bul[ba, lv] = atmos.band_u_lw[lv, ba]
                var_bdl[ba, lv] = atmos.band_d_lw[lv, ba]
                var_bnl[ba, lv] = atmos.band_n_lw[lv, ba]
                var_bus[ba, lv] = atmos.band_u_sw[lv, ba]
                var_bds[ba, lv] = atmos.band_d_sw[lv, ba]
                var_bns[ba, lv] = atmos.band_n_sw[lv, ba]
            end 
        end 

        for lc in 1:atmos.nlev_c
            for ba in 1:atmos.nbands
                var_cfn[ba, lc] = atmos.contfunc_norm[lc, ba]
            end 
        end 

        # ----------------------
        # Close
        close(ds)
        return nothing
    end 

end
