# Contains the radtrans module, which contains all of the core code 
# for setting-up and running SOCRATES

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module radtrans

    # Import libraries
    import NCDatasets
    include("../socrates/julia/src/SOCRATES.jl")

    example_dir     = "socrates/examples/netcdf/CIRC_case6"
    base_name = "case6"

    # Struct for holding data pertaining to the atmosphere
    mutable struct Atmos_t

        # Track state of atmos struct
        is_param::Bool  # Params have been set
        is_alloc::Bool  # Arrays have been allocated
        is_out_lw::Bool # Contains data for LW calculation
        is_out_sw::Bool # Contains data for SW calculation

        # SOCRATES objects
        dimen::SOCRATES.StrDim
        control::SOCRATES.StrCtrl
        spectrum::SOCRATES.StrSpecData
        atm::SOCRATES.StrAtm
        cld::SOCRATES.StrCld
        aer::SOCRATES.StrAer
        bound::SOCRATES.StrBound
        radout::SOCRATES.StrOut

        # Other parameters
        all_channels::Bool 

        spectral_file::String 
        surface_albedo::Float64 
        zenith_degrees::Float64 
        toa_heating::Float64 
        tstar::Float64 
        grav_surf::Float64

        flag_rayleigh::Bool 
        flag_continuum::Bool 
        flag_aerosol::Bool 
        flag_cloud::Bool

        nlev_c::Int32  # Cell centre
        nlev_l::Int32  # Cell edge

        # Temperature/pressure profile 
        tmp::Array
        tmpl::Array
        p::Array 
        pl::Array

        # Calculated fluxes
        flux_d_lw::Array  # downward component, lw 
        flux_u_lw::Array  # upward component, lw
        flux_n_lw::Array  # net upward, lw

        flux_d_sw::Array
        flux_u_sw::Array
        flux_n_sw::Array

        flux_d::Array
        flux_u::Array
        flux_n::Array

        Atmos_t() = new()
    end

    # Deallocate atmosphere arrays
    function dealloc_atmos!(atmos)
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

    # Set parameters of atmosphere
    function setup_atmos!(atmos, spectral_file, all_channels,
                            flag_rayleigh,flag_continuum,flag_aerosol,flag_cloud,
                            zenith_degrees,toa_heating,T_surf,
                            gravity, nlev_centre)

        println("Atmosphere: instantiating SOCRATES objects")
        atmos.dimen =       SOCRATES.StrDim()
        atmos.control =     SOCRATES.StrCtrl()
        atmos.spectrum =    SOCRATES.StrSpecData()
        atmos.atm =         SOCRATES.StrAtm()
        atmos.cld =         SOCRATES.StrCld()
        atmos.aer =         SOCRATES.StrAer()
        atmos.bound =       SOCRATES.StrBound()
        atmos.radout =      SOCRATES.StrOut()
        
        println("Atmosphere: setting parameters")

        # Set parameters
        atmos.spectral_file =   spectral_file
        atmos.flag_rayleigh =   flag_rayleigh
        atmos.flag_continuum =  flag_continuum
        atmos.flag_aerosol =    flag_aerosol
        atmos.flag_cloud =      flag_cloud
        atmos.all_channels =    all_channels

        atmos.nlev_c         =  nlev_centre
        atmos.nlev_l         =  nlev_centre + 1
        atmos.zenith_degrees =  zenith_degrees
        atmos.toa_heating =     toa_heating
        atmos.tstar =           T_surf
        atmos.grav_surf =       gravity

        # Consequential things
        atmos.control.l_gas = true
        atmos.control.l_rayleigh = atmos.flag_rayleigh
        atmos.control.l_continuum = atmos.flag_continuum
        atmos.control.l_aerosol = atmos.flag_aerosol

        # Record that the parameters are set
        atmos.is_param = true
    end
    
    # Allocate atmosphere arrays
    function alloc_atmos!(atmos)

        println("Atmosphere: allocate arrays")

        if !atmos.is_param
            error(" atmosphere parameters have not been set")
        end

        atmos.control.i_cloud_representation = SOCRATES.rad_pcf.ip_cloud_type_homogen
        atmos.cld.n_condensed = 0
        atmos.atm.n_profile = 0

        #########################################
        # spectral data
        #########################################

        atmos.control.spectral_file = atmos.spectral_file
        SOCRATES.set_spectrum(spectrum=atmos.spectrum, 
                              spectral_file=atmos.control.spectral_file, 
                              l_all_gasses=true)

        gas_index = zeros(Int, SOCRATES.gas_list_pcf.npd_gases) # pointers to gases in spectral file
        for i in 1:atmos.spectrum.Gas.n_absorb
            ti = atmos.spectrum.Gas.type_absorb[i]
            gas_index[ti] = i
            println("spectrum.Gas $i type_absorb $ti $(SOCRATES.gas_list_pcf.name_absorb[ti])")
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

        # atm sizes and coordinates 
        atmos.atm.n_layer = npd_layer
        atmos.atm.n_profile = 1
        atmos.atm.lat[1] = 0.0
        atmos.atm.lon[1] = 0.0

        #########################################
        # Temperature and pressure grids
        #########################################

        atmos.tmpl = zeros(Float64, atmos.nlev_l)
        atmos.pl =   zeros(Float64, atmos.nlev_l)
        atmos.tmp =  zeros(Float64, atmos.nlev_c)
        atmos.p =    zeros(Float64, atmos.nlev_c)

        # TEMP
        nc_t = NCDatasets.NCDataset(joinpath(example_dir, base_name*".t"))
        nc_tl = NCDatasets.NCDataset(joinpath(example_dir, base_name*".tl"))

        nc_open = [nc_t, nc_tl]

        atmos.atm.p[1, :] .= nc_t["plev"][:]
        atmos.atm.t[1, :] .= nc_t["t"][1, 1, :]

        atmos.atm.p_level[1, 0:end] .= nc_tl["plev"][:]
        atmos.atm.t_level[1, 0:end] .= nc_tl["tl"][1, 1, :]

        atmos.p[:]    .= nc_t["plev"][:]
        atmos.tmp[:]  .= nc_t["t"][1, 1, :]

        atmos.pl[:]   .= nc_tl["plev"][:]
        atmos.tmpl[:] .= nc_tl["tl"][1, 1, :]

        close.(nc_open)  # close netcdf files

        # /TEMP

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

        atmos.control.l_cloud = atmos.flag_cloud 


        ################################
        # Gaseous absorption
        #################################

        atmos.control.i_gas_overlap = SOCRATES.rad_pcf.ip_overlap_random # = 2
        for j in atmos.control.first_band:atmos.control.last_band
            atmos.control.i_gas_overlap_band[j] = atmos.control.i_gas_overlap
        end

        mr_gases = Dict()
        for i_gas in 1:atmos.spectrum.Gas.n_absorb
            # Read gas mixing ratios
            ti = atmos.spectrum.Gas.type_absorb[i_gas]
            sfx = SOCRATES.input_head_pcf.gas_suffix[ti]
            fn = joinpath(example_dir, base_name*"."*sfx)
            println("Reading mixing ratio for gas $ti $(SOCRATES.gas_list_pcf.name_absorb[ti])    $fn")

            if isfile(fn)
                NCDatasets.NCDataset(fn, "r") do nc
                    mr_gases[sfx] = nc[sfx][1, 1, :]
                    atmos.atm.gas_mix_ratio[1, :, i_gas] .= mr_gases[sfx]
                end
            else
                println("  no file found - setting mixing ratio to 0.0")
                atmos.atm.gas_mix_ratio[:, :, i_gas] .= 0.0
            end
        end

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

        # Mark as allocated
        atmos.is_alloc = true
    end


    function calc_fluxes!(atmos, lw::Bool)

        if lw
            println("Atmosphere: Calculating LW fluxes")
        else
            println("Atmosphere: Calculating SW fluxes")
        end

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

        #   Convert the zenith angles to secants.
        atmos.bound.zen_0[1] = 1.0/cosd(atmos.bound.zen_0[1])

        # Reset dimen.nd_max_order to reduce memory requirements
        atmos.dimen.nd_max_order = 1


        #####################################
        # surface properties
        # see src/aux/assign_surface_char_cdf.f
        # IP_surface_char  = 51, file suffix 'surf'
        #####################################

        atmos.bound.rho_alb[:, SOCRATES.rad_pcf.ip_surf_alb_diff, :] .= atmos.surface_albedo

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
            # 'Is the ir-source function to be ' //      &
            # 'taken as linear or quadratic in tau? (l/q)'
            atmos.control.l_ir_source_quad = false
            # control.l_ir_source_quad = true  # Cl_run_cdf -q
        end


        ######################################################
        # Run radiative transfer model
        ######################################################

        # Gas constant for dry air (FIX ME)
        r_gas_dry          = 287.026 # J K-1 kg-1
        for i in range(atmos.atm.n_layer, 1, step=-1)
            atmos.atm.mass[1, i] = (atmos.atm.p_level[1, i] - atmos.atm.p_level[1, i-1])/atmos.grav_surf
        end

        # dry air density (FIX ME)
        r = r_gas_dry
        for i in 1:atmos.atm.n_layer
            atmos.atm.density[1, i] = atmos.atm.p[1, i]/(r*atmos.atm.t[1, i])
        end

        # Calculation of fluxes (result is stored in atmos.radout )
        SOCRATES.radiance_calc(atmos.control, atmos.dimen, atmos.spectrum, 
                                atmos.atm, atmos.cld, atmos.aer, 
                                atmos.bound, atmos.radout)
    
        # Store new fluxes in atmos struct
        if lw 
            # LW case
            for lv in 1:atmos.nlev_l                # sum over levels
                for ch in 1:atmos.dimen.nd_channel  # sum over channels
                    idx = lv+(ch-1)*atmos.nlev_l
                    atmos.flux_d_lw[lv] += atmos.radout.flux_down[idx]
                    atmos.flux_u_lw[lv] += atmos.radout.flux_up[idx]
                end 
                atmos.flux_n_lw[lv] = atmos.flux_u_lw[lv] - atmos.flux_d_lw[lv] 
            end
            atmos.is_out_lw = true 
        else
            # SW case include direct beam
            for lv in 1:atmos.nlev_l                # sum over levels
                for ch in 1:atmos.dimen.nd_channel  # sum over channels
                    idx = lv+(ch-1)*atmos.nlev_l
                    atmos.flux_d_sw[lv] += atmos.radout.flux_down[idx] + atmos.radout.flux_direct[idx]
                    atmos.flux_u_sw[lv] += atmos.radout.flux_up[idx]
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

    end # end of calc_fluxes

end
