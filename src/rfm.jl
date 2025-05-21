# Contains module and functions for performing radiative transfer with RFM

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module rfm

    # Import modules
    using LoggingExtras

    # Local files
    import ..atmosphere

    # Species supported by RFM (https://eodg.atm.ox.ac.uk/RFM/rfm_spec.html#gaslist)
    const RFM_GASES::Array{String} = [
        # HITRAN
        "H2O", "CO2", "O3", "N2O", "CO", "CH4", "O2", "NO", "SO2", "NO2", "NH3", "HNO3",
        "OH", "HF", "HCl", "HBr", "HI", "ClO", "OCS", "H2CO", "HOCl", "N2", "HCN", "CH3Cl",
        "H2O2", "C2H2", "C2H6", "PH3", "COF2", "SF6", "H2S", "HCOOH", "HO2", "O", "ClONO2",
        "HOBr", "C2H4", "CH3OH", "CH3Br", "CH3CN", "CF4", "C4H2", "HC3N", "H2", "CS",
        "SO3", "C2N2", "COCl2", "SO", "CH3F", "GeH4", "CS2", "CH3I", "NF3", "C3H4", "CH3",
        # Extra
        "ClONO2", "N2O5", "SF6", "CCl4", "HNO4", "SF5CF3", "BrONO2", "ClOOCl", "CH3OH",
        ]

    """
    **Write atmospheric profile for RFM from current state.**

    Parameters:
    - `atmos::atmosphere.Atmos_t`       atmosphere data struct
    """
    function write_profile(atmos::atmosphere.Atmos_t)
        @debug "Write profile for RFM"

        # Convert arrays
        t_arr::Array{Float64} = reverse(atmos.tmpl) # K
        p_arr::Array{Float64} = reverse(atmos.pl) / 1e5  * 1e3 # mbar
        z_arr::Array{Float64} = reverse(atmos.zl) / 1e3 # km

        # file counters
        ir::Int = 0
        nr::Int = 5

        # header
        outstr::String =  ""
        outstr *= "         %d ! Profile Levels\n"%nlev

        # height profile
        outstr *= "*HGT [km]\n"
        for i in 1:atmos.nlev_l
            # line break
            if ir >= nr
                ir=1
                outstr *= "\n"
            else
                ir += 1
            end
            # write next value
            outstr *= @sprintf(" %11.7f",z_arr[i])
        end

        # pressure profile
        ir=0
        outstr *="\n"
        outstr *= "*PRE [mb]\n"
        for i in 1:atmos.nlev_l
            if ir >= nr
                ir=1
                outstr *= "\n"
            else
                ir += 1
            end
            outstr *= @sprintf(" %.6e",p_arr[i])
        end

        # temperature profile
        ir=0
        outstr *= "\n"
        outstr *= "*TEM [K]\n"
        for i in 1:atmos.nlev_l
            if ir >= nr
                ir=1
                outstr*="\n"
            else
                ir += 1
            end
            outstr *= @sprintf(" %8.3f", t_arr[i])
        end


        # gases
        for g in atmos.gas_names

            # get mixing ratio and extend
            x_arr = reverse(atmos.gas_vmr[g])
            push!(x_arr, x_arr[end])

            # clamp range
            clamp!(x_arr, 0.0, 1.0)

            # convert units to ppmv
            x_arr *= 1e6

            # reduce VMR by small fraction so that total VMR is not >1
            x_arr *= 1-1e-10

            # write profile for this  gas
            ir=0
            outstr *="\n"
            outstr *= "*$g [ppmv]\n"
            for i in 1:atmos.nlev_l
                if ir >= nr
                    ir=1
                    outstr*="\n"
                else
                    ir += 1
                end
                outstr *= @sprintf(" %9.7e",x_arr[i])
            end
        end
        outstr *= "\n "

        # write profile to file
        open(joinpath(atmos.rfm_work,"profile.atm"),"w") do hdl
            write(hdl,outstr)
        end
    end

    """
    **Write driver for RFM from current state.**

    See documentation for explanation of each section:
        https://eodg.atm.ox.ac.uk/RFM/sum/rfm_sections.html

    Parameters:
    - `atmos::atmosphere.Atmos_t`       atmosphere data struct
    """
    function write_driver(atmos::atmosphere.Atmos_t)
        @debug "Write driver for RFM"

        # Header
        outstr::String = "*HDR\n  RFM flux observing downwards from TOA \n"

        # Flags
        outstr *= "*FLG\n  NAD RAD SFC \n"

        # Spectral regions
        outstr *= "*SPC\n  1000 7000 1 \n"

        # Gases
        outstr *= "*GAS\n  "
        for g in atmos.gas_names
            if g in RFM_GASES
                outstr *= "$g(1) "
            end
        end
        outstr *= "\n"

        # Path to profile
        outstr *= "*ATM\n  $(joinpath(atmos.rfm_work,"profile.atm")) \n"

        # Secant of viewing angle
        outstr *= "*SEC\n  1"

        # Surface properties
        outstr *= "*SFC\n"
        outstr *= "  temsfc=$(atmos.tmp_surf) \n"
        outstr *= "  emssfc=1.0 \n"
        outstr *= "  hgtsfc=0.0 \n"

        # LbL par file (HITRAN format)
        hitran_path = "/home/n/nichollsh/RFM/par/h2o+co2/00000-12000.par"
        outstr *= "*HIT\n  $hitran_path \n"

        # Output file
        outstr *= "*OUT\n  radfil=$(joinpath(atmos.rfm_work, "fluxes.asc")) \n"

        # End driver table
        outstr *= "*END\n "

        # write driver content to file
        driver_path = joinpath(atmos.rfm_work,"config.drv")
        open(driver_path,"w") do hdl
            write(hdl,outstr)
        end
    end

    """
    **Run RFM radiative transfer with nadir-viewing geometry**
    """
    function run_rfm(atmos::atmosphere.Atmos_t)
        # Write required files
        write_profile(atmos)
        write_driver(atmos)

        # Run subprocess
        @debug "Run RFM"
        execpath = joinpath(atmos.RFM_DIR, "rfm", "rfm")
        # cmd = pipeline(`$execpath`, stdout=devnull)
        cmd = `$execpath`
        run(setenv(cmd, dir=atmos.rfm_work))

        # Get output from RFM
        read_fluxes(atmos::atmosphere.Atmos_t)
    end

    """
    **Read radiative fluxes calculated by RFM**
    """
    function read_fluxes(atmos::atmosphere.Atmos_t)
        @debug "Read fluxes from RFM"

        outpath = joinpath(atmos.rfm_work, "fluxes.asc")
        open(outpath,"r") do f
            while !eof(f)
               # read next line
               s = readline(f)
               println(s)
            end
        end

    end


end # end module
