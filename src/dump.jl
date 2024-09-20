# Write atmosphere to NetCDF or CSV file

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end


module dump

    import ..atmosphere

    using DataStructures
    using NCDatasets
    using LoggingExtras
    using Printf
    using Dates

    """
    Write Pressure vs Temperature & Height profile to a CSV file
    """
    function write_ptz(atmos::atmosphere.Atmos_t, fname::String)

        arr_P, arr_T, arr_Z = atmosphere.get_interleaved_ptz(atmos)

        # Remove old file if exists
        rm(fname, force=true)

        @debug "Writing T(p) csv to $fname"

        open(fname, "w") do f
            write(f, "# pressure  , temperature, height \n")
            write(f, "# [Pa]      , [K]        , [m]  \n")
            for i in 1:atmos.nlev_l+atmos.nlev_c
                @printf(f, "%1.5e , %1.5e , %1.5e \n", arr_P[i], arr_T[i], arr_Z[i])
            end
        end
        return nothing
    end


    """
    Write cell-edge energy fluxes to a CSV file
    """
    function write_fluxes(atmos::atmosphere.Atmos_t, fname::String)

        # Remove old file if exists
        rm(fname, force=true)

        @debug "Writing fluxes CSV to $fname"

        open(fname, "w") do f
            write(f, "# pressure  , U_LW        , D_LW        , N_LW        , U_SW        , D_SW        , N_SW        , U           , D           , N           , convect     , latent      , tot      \n")
            write(f, "# [Pa]      , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]     , [W m-2]  \n")
            for i in 1:atmos.nlev_l
                @printf(f, "%1.5e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e , %+1.4e \n",
                          atmos.pl[i],
                          atmos.flux_u_lw[i], atmos.flux_d_lw[i], atmos.flux_n_lw[i],
                          atmos.flux_u_sw[i], atmos.flux_d_sw[i], atmos.flux_n_sw[i],
                          atmos.flux_u[i],    atmos.flux_d[i],    atmos.flux_n[i],
                          atmos.flux_cdry[i], atmos.flux_l[i],    atmos.flux_tot[i]
                          )
            end
        end
        return nothing
    end

    """
    Write verbose atmosphere data to a NetCDF file
    """
    function write_ncdf(atmos::atmosphere.Atmos_t, fname::String)

        # Create dataset handle
        fname = abspath(fname)
        rm(fname, force=true)

        @debug "Writing NetCDF to $fname"

        # See the tutorial at:
        # https://github.com/Alexander-Barth/NCDatasets.jl#create-a-netcdf-file

        # Note that the content of the NetCDF file is designed to be compatible
        # with what JANUS writes. As a result, they can both be integrated
        # into PROTEUS without compatibility issues.

        # Absorb output from these calls, because they spam the Debug logger
        @debug "ALL DEBUG SUPPRESSED"
        with_logger(MinLevelLogger(current_logger(), Logging.Info-200)) do

            ds = Dataset(fname,"c")

            # Global attributes
            ds.attrib["description"]        = "AGNI atmosphere data"
            ds.attrib["date"]               = Dates.format(now(), "yyyy-u-dd HH:MM:SS")
            ds.attrib["hostname"]           = gethostname()
            ds.attrib["username"]           = ENV["USER"]
            ds.attrib["AGNI_version"]      = atmos.AGNI_VERSION
            ds.attrib["SOCRATES_version"]  = atmos.SOCRATES_VERSION

            plat::String = "Generic"
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
            ngases = atmos.gas_num
            nchars = 16

            defDim(ds, "nlev_c", nlev_c)        # Cell centres
            defDim(ds, "nlev_l", nlev_l)        # Cell edges
            defDim(ds, "ngases", ngases)        # Gases
            defDim(ds, "nchars", nchars)        # Length of string containing gas names
            defDim(ds, "nbands", atmos.nbands)  # Number of spectral bands
            defDim(ds, "nchannels", atmos.dimen.nd_channel)  # Number of spectral channels used for calculations

            # ----------------------
            # Scalar quantities
            #    Create variables
            var_tmp_surf =  defVar(ds, "tmp_surf",      Float64, (), attrib = OrderedDict("units" => "K"))      # Surface brightness temperature [K]
            var_flux_int =  defVar(ds, "flux_int",      Float64, (), attrib = OrderedDict("units" => "W m-2"))  # Internal flux [W m-2]
            var_inst =      defVar(ds, "instellation",  Float64, (), attrib = OrderedDict("units" => "W m-2"))  # Solar flux at TOA
            var_s0fact =    defVar(ds, "inst_factor",   Float64, ())                                            # Scale factor applied to instellation
            var_albbond =   defVar(ds, "bond_albedo",   Float64, ())                                            # Bond albedo used to scale-down instellation
            var_toah =      defVar(ds, "toa_heating",   Float64, (), attrib = OrderedDict("units" => "W m-2"))  # TOA SW BC
            var_tmagma =    defVar(ds, "tmagma",        Float64, (), attrib = OrderedDict("units" => "K"))      # Magma temperature
            var_tmin =      defVar(ds, "tfloor",        Float64, (), attrib = OrderedDict("units" => "K"))      # Minimum temperature
            var_tmax =      defVar(ds, "tceiling",      Float64, (), attrib = OrderedDict("units" => "K"))      # Maximum temperature
            var_plrad =     defVar(ds, "planet_radius", Float64, (), attrib = OrderedDict("units" => "m"))      # Value taken for planet radius
            var_gsurf =     defVar(ds, "surf_gravity",  Float64, (), attrib = OrderedDict("units" => "m s-2"))  # Surface gravity
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
            var_tmp_surf[1] =   atmos.tmp_surf
            var_flux_int[1] =   atmos.flux_int
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
            var_thick =     defVar(ds, "dz",        Float64, ("nlev_c",), attrib = OrderedDict("units" => "m"))
            var_grav =      defVar(ds, "gravity",   Float64, ("nlev_c",), attrib = OrderedDict("units" => "m s-2"))
            var_cp =        defVar(ds, "cp",        Float64, ("nlev_c",), attrib = OrderedDict("units" => "J K-1 kg-1"))
            var_mmw =       defVar(ds, "mmw",       Float64, ("nlev_c",), attrib = OrderedDict("units" => "kg mol-1"))
            var_gases =     defVar(ds, "gases",     Char,    ("nchars", "ngases")) # Transposed cf JANUS because of how Julia stores arrays
            var_x =         defVar(ds, "x_gas",     Float64, ("ngases", "nlev_c"), attrib = OrderedDict("units" => "mol mol-1")) # ^^
            var_cldl  =     defVar(ds, "cloud_mmr", Float64, ("nlev_c",), attrib = OrderedDict("units" => "kg kg-1"))
            var_fdl =       defVar(ds, "fl_D_LW",   Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
            var_ful =       defVar(ds, "fl_U_LW",   Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
            var_fnl =       defVar(ds, "fl_N_LW",   Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
            var_fds =       defVar(ds, "fl_D_SW",   Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
            var_fus =       defVar(ds, "fl_U_SW",   Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
            var_fns =       defVar(ds, "fl_N_SW",   Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
            var_fd =        defVar(ds ,"fl_D",      Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
            var_fu =        defVar(ds, "fl_U",      Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
            var_fn =        defVar(ds, "fl_N",      Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
            var_fcd =       defVar(ds, "fl_cnvct",  Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
            var_fcc =       defVar(ds, "fl_cndct",  Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
            var_fla =       defVar(ds, "fl_latent", Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
            var_ft =        defVar(ds, "fl_tot",    Float64, ("nlev_l",), attrib = OrderedDict("units" => "W m-2"))
            var_fdiff =     defVar(ds, "fl_dif",    Float64, ("nlev_c",), attrib = OrderedDict("units" => "W m-2"))
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
            var_albr =      defVar(ds, "albedo_r",  Float64, ("nbands",))
            var_albe =      defVar(ds, "albedo_e",  Float64, ("nbands",))

            #     Store data
            var_p[:] =      atmos.p
            var_pl[:] =     atmos.pl
            var_tmp[:] =    atmos.tmp
            var_tmpl[:] =   atmos.tmpl
            var_z[:]    =   atmos.z
            var_zl[:]   =   atmos.zl
            var_mmw[:]  =   atmos.layer_mmw
            var_cp[:]  =    atmos.layer_cp
            var_grav[:]  =  atmos.layer_grav
            var_thick[:]  = atmos.layer_thick

            # Composition
            for (i_gas,gas) in enumerate(atmos.gas_names)
                # Fill gas names
                for i_char in 1:nchars
                    var_gases[i_char, i_gas] = ' '
                end
                for i_char in 1:length(atmos.gas_names[i_gas])
                    var_gases[i_char,i_gas] = atmos.gas_names[i_gas][i_char]
                end

                # Fill VMR
                for i_lvl in 1:nlev_c
                    var_x[i_gas, i_lvl] = atmos.gas_vmr[gas][i_lvl]
                end
            end

            # Clouds
            var_cldl[:] =   atmos.cloud_arr_l

            # Kzz mixing
            var_kzz[:] =    atmos.Kzz

            # Bolometric fluxes
            var_fdl[:] =    atmos.flux_d_lw
            var_ful[:] =    atmos.flux_u_lw
            var_fnl[:] =    atmos.flux_n_lw

            var_fds[:] =    atmos.flux_d_sw
            var_fus[:] =    atmos.flux_u_sw
            var_fns[:] =    atmos.flux_n_sw

            var_fd[:] =     atmos.flux_d
            var_fu[:] =     atmos.flux_u
            var_fn[:] =     atmos.flux_n

            var_fcd[:] =    atmos.flux_cdry
            var_fcc[:] =    atmos.flux_cdct
            var_fla[:] =    atmos.flux_l

            var_ft[:] =     atmos.flux_tot
            var_hr[:] =     atmos.heating_rate
            var_fdiff[:] =  atmos.flux_dif

            # Spectral fluxes
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

            # Surface spectral albedo
            var_albr[:] = atmos.surf_r_arr
            var_albe[:] = atmos.surf_e_arr

            close(ds)

        end # suppress output
        @debug "ALL DEBUG RESTORED"

        return nothing
    end # end write_ncdf


end
