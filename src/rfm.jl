# Contains module and functions for performing radiative transfer with RFM

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module rfm

    # Import modules
    using LoggingExtras
    using Printf

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

    # Executable paths
    const RFM_LINUX::String = abspath(dirname(@__FILE__),
                                            "..", "res", "blobs", "rfm-amd64-linux")
    const RFM_MACOS::String = abspath(dirname(@__FILE__),
                                            "..", "res", "blobs", "rfm-arm64-macos")

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
        z_arr::Array{Float64} = reverse(atmos.rl .- atmos.rl[end]) / 1e3 # km

        # file counters
        ir::Int = 0
        nr::Int = 5

        # header
        outstr::String =  ""
        outstr *= "         $(atmos.nlev_l) ! Profile Levels\n"

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
            x_arr *= 1-1e-7

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
    - `numin::Float64`                  minimum wavenumber [cm-1]
    - `numax::Float64`                  maximum wavenumber [cm-1]
    - `nures::Float64`                  resolution on wavenumber [cm-1]
    """
    function write_driver(atmos::atmosphere.Atmos_t,
                            numin::Float64, numax::Float64, nures::Float64)

        @debug "Write driver file for RFM"

        # Validate wavenumber parameters
        if numin > numax
            numin, numax = numax, numin
        end
        if nures > 1.0
            @warn "Clipping wavenumber resolution to 1 cm^-1"
            nures = 1.0
        end
        numax = max(numax, numin+nures*2) # must do at least two points

        # Header
        outstr::String = "*HDR\n  RFM flux observing downwards from TOA \n"

        # Flags
        outstr *= "*FLG\n  NAD RAD SFC \n"

        # Spectral regions
        outstr *= @sprintf("*SPC\n  %.4f %.4f %.4f \n", numin, numax, nures)

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
        outstr *= "*SEC\n  $(secd(atmos.zenith_degrees)) \n"

        # Surface properties
        outstr *= "*SFC\n"
        outstr *= "  temsfc=$(atmos.tmp_surf) \n"
        outstr *= "  emssfc=$(1-atmos.albedo_s) \n"
        outstr *= "  hgtsfc=0.0 \n"

        # LbL par file (HITRAN format)
        outstr *= "*HIT\n  $(atmos.rfm_parfile) \n"

        # Output file
        outstr *= "*OUT\n  radfil=$(joinpath(atmos.rfm_work, "fluxes.asc")) \n"

        # End driver table
        outstr *= "*END\n "

        # write driver content to file
        driver_path = joinpath(atmos.rfm_work,"rfm.drv")
        open(driver_path,"w") do hdl
            write(hdl,outstr)
        end
    end

    """
    **Run RFM radiative transfer with nadir-viewing geometry**

    Parameters:
    - `atmos::atmosphere.Atmos_t`       atmosphere data struct
    - `numin::Float64`                  minimum wavenumber [cm-1]
    - `numax::Float64`                  maximum wavenumber [cm-1]
    - `nures::Float64`                  resolution on wavenumber [cm-1]
    """
    function run_rfm(atmos::atmosphere.Atmos_t,
                        numin::Float64, numax::Float64; nures::Float64=1.0)
        # Write required files
        write_profile(atmos)
        write_driver(atmos, numin, numax, nures)

        # Locate executable
        if Sys.isapple()
            @debug "Run RFM (MacOS binary)"
            execpath = RFM_MACOS
        else
            @debug "Run RFM (Linux binary)"
            execpath = RFM_LINUX
        end

        # Construct command
        cmd = `$execpath` # path to RFM executable
        cmd = setenv(cmd, dir=atmos.rfm_work) # run inside working directory
        cmd = pipeline(cmd, stdout=devnull) # hide rfm terminal output

        # Run subprocess
        run(cmd)

        # Get output from RFM
        read_fluxes(atmos::atmosphere.Atmos_t)
    end

    """
    **Read radiative fluxes calculated by RFM**
    """
    function read_fluxes(atmos::atmosphere.Atmos_t)
        @debug "Read fluxes from RFM"

        # output file path
        outpath = joinpath(atmos.rfm_work, "fluxes.asc")
        if !isfile(outpath)
            @error "Could not find RFM output file: $outpath"
            return
        end

        # file header info
        numin::Float64  = 0.0       # min wavenumber from file [cm-1]
        nures::Float64  = 0.0       # wavenumber resolution from file [cm-1]

        # read the file
        headstr::String = ""   # ASCII, header info
        datastr::String = ""   # ASCII, full flux data
        open(outpath,"r") do f
            # loop through lines
            idx::Int = -1
            linestr::String = ""   # string containing content from line
            while !eof(f)
                # read next line
                linestr = readline(f)

                # skip comments
                if startswith(linestr,'!')
                    continue
                end
                idx += 1

                # first line is header
                if idx == 0
                    headstr = linestr
                    continue
                end

                # all other lines are data ...
                datastr *= chomp(strip(linestr))
            end
        end # close fluxes file

        # parse header
        atmos.rfm_npts = parse(Int, split(headstr)[1])
        numin = parse(Float64, split(headstr)[2])
        nures = parse(Float64, split(headstr)[3])

        # check grid regularity
        if nures <= 0
            @error "Cannot parse irregular wavenumber grid"
            return
        end

        # set wavenumber array
        atmos.rfm_wn = collect(Float64, range(start=numin,length=atmos.rfm_npts,step=nures))

        # set flux array
        atmos.rfm_fl = zero(atmos.rfm_wn)
        for (idx,val) in enumerate(split(datastr))
            atmos.rfm_fl[idx] = parse(Float64, val)
        end

        # convert flux units to [erg/(s cm2 cm-1)]
        atmos.rfm_fl *= 1e-9 * 1e7 * pi
    end

end # end module
