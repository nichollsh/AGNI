# Contains the atmosphere module, which contains all of the core code 
# for setting-up the atmosphere and running SOCRATES

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module atmosphere

    # Import libraries
    include("../socrates/julia/src/SOCRATES.jl")

    using Printf
    using Revise

    import phys

    # Contains data pertaining to the atmosphere (fluxes, temperature, etc.)
    mutable struct Atmos_t

        # Track state of atmos struct
        is_param::Bool  # Params have been set
        is_alloc::Bool  # Arrays have been allocated
        is_out_lw::Bool # Contains data for LW calculation
        is_out_sw::Bool # Contains data for SW calculation

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

        # SOCRATES parameters
        all_channels::Bool 
        spectral_file::String 
        albedo_s::Float64 
        zenith_degrees::Float64 
        toa_heating::Float64 
        tstar::Float64 
        grav_surf::Float64

        # Pressure-temperature grid 
        nlev_c::Int32  # Cell centre (count)
        nlev_l::Int32  # Cell edge (count)

        p_boa::Float64 # bottom 
        p_toa::Float64 # top 

        tmp::Array      
        tmpl::Array     
        p::Array        
        pl::Array      

        T_floor::Float64        # Temperature floor [K]

        # Mixing ratios 
        mr_hom::Dict      # Dictonary with scalar quantities (well-mixed)
        mr_het::Dict      # Dictionary with per-level values (not yet implemented)

        # Layers' average properties
        layer_density::Array  # density [kg m-3]
        layer_mmw::Array      # mean molecular weight [kg mol-1]
        layer_cp::Array       # heat capacity at const-p [J K-1 mol-1]
        layer_grav::Array     # gravity [m s-2]

        # Calculated fluxes (W m-2)
        flux_d_lw::Array  # down component, lw 
        flux_u_lw::Array  # up component, lw
        flux_n_lw::Array  # net upward, lw

        flux_d_sw::Array  # down component, sw 
        flux_u_sw::Array  # up component, sw
        flux_n_sw::Array  # net upward, sw

        flux_d::Array    # down component, lw+sw 
        flux_u::Array    # up component, lw+sw 
        flux_n::Array    # net upward, lw+sw 

        # Heating rate 
        heating_rate::Array # radiative heating rate [K/day]

        Atmos_t() = new()
    end

    # Deallocate atmosphere arrays
    function deallocate!(atmos::atmosphere.Atmos_t)
        println("Atmosphere: dellocating arrays")

        SOCRATES.deallocate_atm(     atmos.atm)

        SOCRATES.deallocate_cld(     atmos.cld)
        SOCRATES.deallocate_cld_prsc(atmos.cld)

        SOCRATES.deallocate_aer(     atmos.aer)
        SOCRATES.deallocate_aer_prsc(atmos.aer)

        SOCRATES.deallocate_bound(   atmos.bound)
        SOCRATES.deallocate_out(     atmos.radout)

        atmos.is_alloc = false
    end

    """
    **Set parameters of the atmosphere.**

    This is used to set up the struct immediately after creation. It must be 
        done before `allocate!()` is called, and before any RT calculcations
        are performed.

    Arguments:
    - `atmos::Atmos_t`                  the atmosphere struct instance to be used.
    - `ROOT_DIR::String`                AGNI root directory 
    - `OUT_DIR::String`                 Output directory
    - `spfile_name::String`             name of spectral file
    - `toa_heating::Float64`            downward shortwave flux at the top of the atmosphere [W m-2].
    - `tstar::Float64`                  effective surface temperature to provide upward longwave flux at the bottom of the atmosphere [K].
    - `gravity::Float64`                gravitational acceleration at the surface [m s-2].
    - `nlev_centre::Int64`              number of cells.
    - `p_surf::Float64`                 total surface pressure [bar].
    - `p_top::Float64`                  total top of atmosphere pressure [bar].
    - `mixing_ratios::Dict`             dictionary of mixing ratios in the format (key,value)=(gas,mixing_ratio).
    - `zenith_degrees::Float64=54.74`   angle of radiation from the star, relative to the zenith [degrees].
    - `albedo_s::Float64=0.0`           surface albedo.
    - `all_channels::Bool=true`         use all channels available for RT?
    - `flag_rayleigh::Bool=false`       include rayleigh scattering?
    - `flag_gcontinuum::Bool=false`     include generalised continuum absorption?
    - `flag_continuum::Bool=false`      include continuum absorption?
    - `flag_aerosol::Bool=false`        include aersols?
    - `flag_cloud::Bool=false`          include clouds?
    """
    function setup!(atmos::atmosphere.Atmos_t, 
                    ROOT_DIR::String, OUT_DIR::String, 
                    spfile_name::String, 
                    toa_heating::Float64, tstar::Float64,
                    gravity::Float64, nlev_centre::Int64, p_surf::Float64, p_top::Float64,
                    mixing_ratios::Dict;
                    zenith_degrees::Float64 =   54.74,
                    albedo_s::Float64 =         0.0,
                    all_channels::Bool  =       true,
                    flag_rayleigh::Bool =       false,
                    flag_gcontinuum::Bool =     false,
                    flag_continuum::Bool =      false,
                    flag_aerosol::Bool =        false,
                    flag_cloud::Bool =          false
                    )

        atmos.dimen =       SOCRATES.StrDim()
        atmos.control =     SOCRATES.StrCtrl()
        atmos.spectrum =    SOCRATES.StrSpecData()
        atmos.atm =         SOCRATES.StrAtm()
        atmos.cld =         SOCRATES.StrCld()
        atmos.aer =         SOCRATES.StrAer()
        atmos.bound =       SOCRATES.StrBound()
        atmos.radout =      SOCRATES.StrOut()

        # Set parameters
        atmos.ROOT_DIR = abspath(ROOT_DIR)
        atmos.OUT_DIR = abspath(OUT_DIR)

        if spfile_name == "_runtime"
            atmos.spectral_file =  joinpath([atmos.OUT_DIR, ".spfile_runtime"])
        else
            atmos.spectral_file =  joinpath([atmos.ROOT_DIR, "res", "spectral_files", spfile_name, spfile_name])
        end
        atmos.all_channels =    all_channels

        atmos.T_floor =          5.0 

        atmos.nlev_c         =  max(nlev_centre,10)
        atmos.nlev_l         =  nlev_centre + 1
        atmos.zenith_degrees =  max(min(zenith_degrees,90.0), 0.0)
        atmos.albedo_s =        max(min(albedo_s, 1.0 ), 0.0)
        atmos.toa_heating =     max(toa_heating, 0.0)
        atmos.tstar =           max(tstar, atmos.T_floor)
        atmos.grav_surf =       gravity

        if p_top > p_surf 
            error("p_top must be less than p_surf")
        end

        atmos.p_toa =           p_top * 1.0e+5 # Convert bar -> Pa
        atmos.p_boa =           p_surf * 1.0e+5

        atmos.control.l_gas =       true
        atmos.control.l_rayleigh =  flag_rayleigh
        atmos.control.l_continuum = flag_continuum
        atmos.control.l_cont_gen =  flag_gcontinuum
        atmos.control.l_aerosol =   flag_aerosol
        atmos.control.l_cloud =     flag_cloud

        # Normalise and store gas mixing ratios in dictionary
        if isempty(mixing_ratios)
            error("No mixing ratios provided")
        end

        atmos.mr_hom = Dict()
        norm_factor = sum(values(mixing_ratios))
        
        for (key, value) in mixing_ratios
            if key in SOCRATES.input_head_pcf.header_gas
                atmos.mr_hom[key] = value / norm_factor
            else
                error("Invalid gas '$key'")
            end
        end

        # Record that the parameters are set
        atmos.is_param = true
    end

    """
    **Get the mass mixing ratio of a gas within the atmosphere.**

    The atmosphere is assumed to be well-mixed.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    - `gas::String`             name of the gas (e.g. "H2O").

    Returns:
    - `mr::Float64`             mixing ratio of `gas`.
    """
    function get_mr(atmos::atmosphere.Atmos_t, gas::String)::Float64

        gas_valid = strip(gas, ' ')
        gas_valid = uppercase(gas_valid)

        mr = 0.0
        if gas_valid in keys(atmos.mr_hom)
            mr = atmos.mr_hom[gas_valid]
        end 

        return mr
    end

    """
    **Calculate properties within each layer of the atmosphere (e.g. density, mmw).**

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    """
    function calc_layer_props!(atmos::atmosphere.Atmos_t)
        if !atmos.is_param
            error(" atmosphere parameters have not been set")
        end

        # Gravity (TODO: Make height-dependent)
        atmos.layer_grav = ones(Float64, atmos.nlev_c) * atmos.grav_surf

        # Mass
        for i in 1:atmos.atm.n_layer
            atmos.atm.mass[1, i] = (atmos.atm.p_level[1, i] - atmos.atm.p_level[1, i-1])/atmos.layer_grav[i]
        end

        # mmw, cp, rho
        atmos.layer_mmw     = zeros(Float64, atmos.nlev_c)
        atmos.layer_density = zeros(Float64, atmos.nlev_c)
        atmos.layer_cp      = zeros(Float64, atmos.nlev_c)
        for i in 1:atmos.atm.n_layer

            # for each gas
            for gas in keys(atmos.mr_hom) 

                # set mmw
                if gas in keys(atmos.mr_hom) 
                    atmos.layer_mmw[i] += atmos.mr_hom[gas] * phys.lookup_mmw[gas]
                end

                # set cp
                if gas in keys(atmos.mr_hom) 
                    atmos.layer_cp[i] += atmos.mr_hom[gas] * phys.lookup_cp[gas] * phys.lookup_mmw[gas]
                end

            end

            # density
            atmos.layer_density[i] = ( atmos.p[i] * atmos.layer_mmw[i] )  / (phys.R_gas * atmos.tmp[i]) 
            atmos.atm.density[1,i] = atmos.layer_density[i]
        end

    end 

    """
    **Generate a log-spaced pressure grid.**

    Uses the values of `atmos.nlev_l`, `atmos.p_toa`, and `atmos.p_boa`.

    Arguments:
    - `atmos::Atmos_t`          the atmosphere struct instance to be used.
    """
    function generate_pgrid!(atmos::atmosphere.Atmos_t)
        # Set pressure cell-edge array by logspace
        atmos.pl = 10 .^ range( log10(atmos.p_toa), stop=log10(atmos.p_boa), length=atmos.nlev_l)

        # Set pressure cell-centre array by interpolation
        atmos.p =    zeros(Float64, atmos.nlev_c)
        atmos.p[:] .= 0.5 .* (atmos.pl[2:end] + atmos.pl[1:end-1])
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

        atmos.control.i_cloud_representation = SOCRATES.rad_pcf.ip_cloud_type_homogen
        atmos.cld.n_condensed = 0
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
            co2_mr = get_mr(atmos, "co2")
            n2_mr  = get_mr(atmos, "n2")
            h2o_mr = get_mr(atmos, "h2o")
            runfile = joinpath([atmos.ROOT_DIR, "src", "insert_rayleigh.py"])
            run(`python $runfile $spectral_file_run $co2_mr $n2_mr $h2o_mr`)
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
        SOCRATES.allocate_cld(  atmos.cld,   atmos.dimen, atmos.spectrum)
        SOCRATES.allocate_aer(  atmos.aer,   atmos.dimen, atmos.spectrum)
        SOCRATES.allocate_bound(atmos.bound, atmos.dimen, atmos.spectrum)

        atmos.bound.rho_alb[:, SOCRATES.rad_pcf.ip_surf_alb_diff, :] .= atmos.albedo_s   #   Set surface albedo.
        atmos.bound.zen_0[1] = 1.0/cosd(atmos.bound.zen_0[1])   #   Convert the zenith angles to secants.

        # atm sizes and coordinates 
        atmos.atm.n_layer = npd_layer
        atmos.atm.n_profile = 1
        atmos.atm.lat[1] = 0.0
        atmos.atm.lon[1] = 0.0

        #########################################
        # Temperature and pressure grids
        #########################################

        # Initialise temperature grid to be isothermal
        atmos.tmpl = ones(Float64, atmos.nlev_l) .* atmos.tstar
        atmos.tmp =  ones(Float64, atmos.nlev_c) .* atmos.tstar

        # Initialise pressure grid with current p_toa and p_boa
        generate_pgrid!(atmos)

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

        atmos.control.i_gas_overlap = SOCRATES.rad_pcf.ip_overlap_random # = 2
        for j in atmos.control.first_band:atmos.control.last_band
            atmos.control.i_gas_overlap_band[j] = atmos.control.i_gas_overlap
        end

        # Set mixing ratios (assumed to be well-mixed)
        for i_gas in 1:atmos.spectrum.Gas.n_absorb
            ti = atmos.spectrum.Gas.type_absorb[i_gas]
            gas = SOCRATES.input_head_pcf.header_gas[ti]

            if gas in keys(atmos.mr_hom)
                atmos.atm.gas_mix_ratio[1, :, i_gas] .= atmos.mr_hom[gas]
            else
                atmos.atm.gas_mix_ratio[:, :, i_gas] .= 0.0
            end
        end

        # Warn user if we are missing gases 
        missing_gases = false
        for i in 1:atmos.nlev_c
            lvl_tot = 0.0
            for i_gas in 1:atmos.spectrum.Gas.n_absorb
                lvl_tot += atmos.atm.gas_mix_ratio[1, i, i_gas]
            end 
            missing_gases = missing_gases || (lvl_tot < 1.0)
        end 
        if missing_gases
            println("WARNING: mixing ratios do not sum to unity (possibly due to missing gases in spectral file)")
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
            error("Clouds not implemented")
        else
            atmos.control.i_cloud = SOCRATES.rad_pcf.ip_cloud_off # 5 (clear sky)
        end

        SOCRATES.allocate_cld_prsc(atmos.cld, atmos.dimen, atmos.spectrum)

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

        atmos.heating_rate =      zeros(Float64, atmos.nlev_c)

        # Mark as allocated
        atmos.is_alloc = true
    end  # end of allocate

    function radtrans!(atmos::atmosphere.Atmos_t, lw::Bool)
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
            Bool(atmos.spectrum.Basic.l_present[6]) ||
                error("The spectral file contains no data for the Planckian function." )

            if Bool(atmos.spectrum.Basic.l_present[2])
                atmos.control.l_solar_tail_flux = true
            end

        else
            Bool(atmos.spectrum.Basic.l_present[2]) ||
                error("The spectral file contains no solar spectral data.")
             
            atmos.bound.zen_0[1] = atmos.zenith_degrees # Assign the solar zenith angle
            atmos.bound.solar_irrad[1] = atmos.toa_heating   # The file of solar irradiances.
        end

        # Set the two-stream approximation to be used
        if lw
            atmos.control.i_2stream = 12 # -t 12, as per UKMO recommendation
        else
            atmos.control.i_2stream = 16 # as per Cl_run_cdf 
        end

        #####################################
        # Angular integration
        # see src/aux/angular_control_cdf.f
        #####################################

        atmos.control.l_rescale = false # Cl_run_cdf default  (+R for true)
        if atmos.control.l_rescale
            atmos.control.l_henyey_greenstein_pf = true
        end

        atmos.control.i_solver = 13 # -v 13: 
        # the solver used for the two-stream calculations. 
        # 13 is recommended for clear-sky, 
        # 16 is recommended for cloudy-sky,
        # 17 is recommended for cloud with separate stratiform and convective regions.

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
        # Treatment of scattering
        ###################################################

        atmos.control.i_scatter_method = SOCRATES.rad_pcf.ip_scatter_full
        for i in atmos.control.first_band:atmos.control.last_band
            atmos.control.i_scatter_method_band[i] = atmos.control.i_scatter_method
        end

        ####################################################
        # Surface temperatures
        ####################################################

        if lw
            atmos.bound.t_ground[1] = atmos.tstar
        end

        #####################################################
        # Variation of the temperature within layers
        ####################################################

        if lw
            atmos.control.l_ir_source_quad = true
        end

        ######################################################
        # Run radiative transfer model
        ######################################################

        clamp!(atmos.tmp, atmos.T_floor, Inf)
        clamp!(atmos.tmpl, atmos.T_floor, Inf)

        atmos.atm.p[1, :] .= atmos.p[:]
        atmos.atm.t[1, :] .= atmos.tmp[:]
        atmos.atm.p_level[1, 0:end] .= atmos.pl[:]
        atmos.atm.t_level[1, 0:end] .= atmos.tmpl[:]
        
        calc_layer_props!(atmos)

        SOCRATES.radiance_calc(atmos.control, atmos.dimen, atmos.spectrum, 
                                atmos.atm, atmos.cld, atmos.aer, 
                                atmos.bound, atmos.radout)
    
        # Store new fluxes in atmos struct
        if lw 
            # LW case
            atmos.flux_u_lw[:] .= 0.0
            atmos.flux_d_lw[:] .= 0.0
            atmos.flux_n_lw[:] .= 0.0
            for lv in 1:atmos.nlev_l                # sum over levels
                for ch in 1:atmos.dimen.nd_channel  # sum over channels
                    idx = lv+(ch-1)*atmos.nlev_l
                    atmos.flux_d_lw[lv] += max(0.0, atmos.radout.flux_down[idx])
                    atmos.flux_u_lw[lv] += max(0.0, atmos.radout.flux_up[idx])
                end 
                atmos.flux_n_lw[lv] = atmos.flux_u_lw[lv] - atmos.flux_d_lw[lv] 
            end
            atmos.is_out_lw = true 
        else
            # SW case
            atmos.flux_u_sw[:] .= 0.0
            atmos.flux_d_sw[:] .= 0.0
            atmos.flux_n_sw[:] .= 0.0
            for lv in 1:atmos.nlev_l                # sum over levels
                for ch in 1:atmos.dimen.nd_channel  # sum over channels
                    idx = lv+(ch-1)*atmos.nlev_l
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

    end # end of radtrans

    # Get interleaved cell-centre and cell-edge PT grid
    function get_interleaved_pt(atmos)
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

    # Write current interleaved pressure/temperature grid to a file
    function write_pt(atmos, fname)

        arr_P, arr_T = atmosphere.get_interleaved_pt(atmos)
        arr_P *= 1e-5
        
        rm(fname, force=true)

        open(fname, "w") do f
            write(f, "# pressure  , temperature \n")
            write(f, "# [bar]     , [K] \n")
            for i in 1:atmos.nlev_l+atmos.nlev_c
                @printf(f, "%1.5e , %1.5e \n", arr_P[i], arr_T[i])
            end
        end

    end

    # Write current cell-edge fluxes to a file
    function write_fluxes(atmos, fname)

        arr_P = atmos.pl * 1e-5
        
        rm(fname, force=true)

        open(fname, "w") do f
            write(f, "# pressure  , U_LW        , D_LW        , N_LW        , U_SW        , D_SW        , N_SW        , U           , D           , N       \n")
            write(f, "# [bar]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2] \n")
            for i in 1:atmos.nlev_l
                @printf(f, "%1.5e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e \n", 
                          arr_P[i], 
                          atmos.flux_u_lw[i], atmos.flux_d_lw[i], atmos.flux_n_lw[i],
                          atmos.flux_u_sw[i], atmos.flux_d_sw[i], atmos.flux_n_sw[i],
                          atmos.flux_u[i],    atmos.flux_d[i],    atmos.flux_n[i],
                          )
            end
        end

    end

end
