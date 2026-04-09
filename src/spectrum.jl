# Contains the spectrum module, for updating the SOCRATES spectral file at
# runtime with the solar flux / thermal source function, and scattering.

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end


module spectrum

    using LoggingExtras
    using Printf
    using LinearAlgebra
    import Interpolations: interpolate, Gridded, Linear, Flat, extrapolate, Extrapolation
    import DelimitedFiles:readdlm

    import ..phys

    include(joinpath(ENV["RAD_DIR"], "julia", "gen", "input_head_pcf.jl"))

    # Constants
    const FLOAT_SML = 1.0e-45
    const FLOAT_BIG = 1.0e45

    """
    **Get the version of SOCRATES being used.**

    Returns:
    - `version::String`  Version string from SOCRATES `version` file
    """
    function get_socrates_version()::String
        return readchomp(joinpath(ENV["RAD_DIR"],"version"))
    end

    """
    **Count the number of gaseous absorbers in a SOCRATES spectral file.**

    Arguments:
    - `spec_file::String`   Path to spectral file

    Returns:
    - `num_gases::Int64`    Number of gaseous absorbers, or -1 if not found
    """
    function count_gases(spec_file::String)::Int

        if !isfile(spec_file)
            @warn "Spectral file not found: '$spec_file'"
            return -1
        end

        # Read file and search for the gas count line
        @debug "Counting gaseous absorbers in spectral file: '$spec_file'"
        lines = readlines(spec_file)
        for line in lines
            if contains(line, "Total number of gaseous absorbers")
                # Extract the number from the line
                # Format: "Total number of gaseous absorbers =     6"
                parts = split(line, "=")
                if length(parts) == 2
                    num_str = strip(parts[2])
                    try
                        return parse(Int, num_str)
                    catch
                        @warn "Could not parse gas count from line: '$line'"
                        return -1
                    end
                else
                    @warn "Unexpected format for gas count line: '$line'"
                    return -1
                end
            end
        end

        @warn "Could not find 'Total number of gaseous absorbers' in spectral file"
        return -1
    end

    """
    **Insert aerosol header information into a SOCRATES spectral file.**

    Arguments:
    - `work_file::String`          Path to spectral file to modify in-place
    - `species::Array{String,1}`   List of aerosol species names

    Returns:
    - `success::Bool`              function executed successfully
    """
    function insert_aerosol_header(work_file::String, species::Array{String,1})::Bool

        # Read file into memory
        lines = readlines(work_file)

        # Check that file does not already contain aerosol header information
        for l in lines
            if contains(l, "aerosols")
                @warn "Spectral file already contains aerosol header information"
                return false
            end
        end

        # Loop through aerosols
        block0_lines = [
            "List of indexing numbers of aerosols.",
            "Index       Aerosol(type number and name)",
        ]
        num_aer = 0
        for (i,s) in enumerate(input_head_pcf.aerosol_suffix)
            if s in species
                num_aer += 1
                push!(block0_lines,
                        @sprintf("   %2d          %2d       %s ",
                        num_aer, i, strip(input_head_pcf.aerosol_title[i]))
                    )
            end
        end

        # List number of aerosols near the top of file
        insert!(lines, 4, "Total number of aerosols =    $num_aer")

        # Match previous insertion before first `*BLOCK: TYPE =    1`.
        block_idx = -1
        for l in lines
            if startswith(l, "*END")
                block_idx = findfirst(==(l), lines)
                break
            end
        end
        if block_idx < 1
            @warn "Could not find first '*END' when injecting aerosol header"
            return false
        end

        # Insert aerosol header lines into file
        foreach(l -> insert!(lines, block_idx, l), reverse(block0_lines))

        # Overwrite with modified file
        open(work_file, "w") do f
            write(f, join(lines, "\n"))
            write(f, "\n")
        end
        return true
    end

    """
    **Generate blackbody stellar spectrum.**

    Arguments:
    - `Teff::Float64`       Star's photospheric effective temperature
    - `S0::Float64`         Bolometric instellation received by planet [W m-2]

    Returns:
    - `wl::Array`           Wavelength array [nm]
    - `fl::Array`           Flux array [erg s-1 cm-2 nm-1]
    """
    function blackbody_star(Teff::Float64, S0::Float64)::Tuple{Array{Float64,1}, Array{Float64,1}}

        @debug "Generate blackbody stellar spectrum"

        # Generate wavelength array (1 nm linear spacing from 0.5 nm to 100 μm)
        wl::Array{Float64, 1} = collect(Float64, range(start=1.0, stop=100e3, step=1.0))
        fl::Array{Float64, 1} = zero(wl)

        # Calculate fluxes at star surface [W m-2 nm-1]
        fl[:] .= phys.evaluate_planck.(wl[:], Teff)

        # Convert to [erg s-1 cm-2 nm-1]
        fl *= 1e7

        # Calculate scaling factor to account for Rstar and Separation
        #    S0 = σT^4 * (R/a)^2
        # Apply scaling factor to fl array, to get spectrum at planet's TOA
        fl *= S0 / ( phys.σSB * Teff^4)

        # Limit range on fl to prevent numerical problems
        clamp!(fl, FLOAT_SML, FLOAT_BIG)

        return wl, fl
    end

    """
    **Load stellar spectrum from a text file.**

    Doesn't matter where the flux is scaled to, because SOCRATES will normalise it.

    Arguments:
    - `path::String`        Path to the file

    Returns:
    - `wl::Array`           Wavelength array [nm]
    - `fl::Array`           Flux array [erg s-1 cm-2 nm-1]
    """
    function load_from_file(path::String)::Tuple{Array{Float64,1}, Array{Float64,1}}

        @debug "Read stellar spectrum from file"

        if isfile(path)
            spec_data = readdlm(abspath(path), Float64; header=false, skipstart=2)
        else
            @warn "Cannot find stellar spectrum at '$path'"
            return zeros(Float64, 1), zeros(Float64, 1)
        end

        return spec_data[:,1], spec_data[:,2]
    end


    """
    **Write a stellar spectrum in the SOCRATES format.**

    Will down-bin spectrum if the resolution is too high.

    Arguments:
    - `wl::Array`           Wavelength array [nm]
    - `fl::Array`           Flux array [erg s-1 cm-2 nm-1]
    - `star_file::String`   Path to output file
    - `nbins_max::Int64`    Maximum number of points in the spectrum

    Returns:
    - `success::Bool`       function executed successfully
    """
    function write_to_socrates_format(wl::Array{Float64,1}, fl::Array{Float64,1},
                                        star_file::String, nbins_max::Int64=99900)::Bool

        len_wl::Int64 = length(wl)
        len_fl::Int64 = length(fl)
        socrates_nbins_max::Int64 = Int(1e5 - 3)  # do not change this
        tgt_bins::Int64 = min(nbins_max, len_wl, socrates_nbins_max)

        # Validate
        if len_wl != len_fl
            @warn "Stellar wavelength and flux arrays have different lengths"
            return false
        end
        if len_wl < 2
            @warn "Loaded stellar spectrum is too short!"
            return false
        end
        if minimum(wl) < FLOAT_SML
            @warn "Minimum wavelength is too small"
            return false
        end
        if wl[2] < wl[1]
            @warn "Stellar wavelength array must be strictly ascending"
            return false
        end
        clamp!(fl, FLOAT_SML, FLOAT_BIG)  # Clamp values

        # warn about non-critical problems
        if len_wl < 500
            @warn "Loaded stellar spectrum is very short!"
        end

        # Ensure that wl array has no duplicates
        # https://discourse.julialang.org/t/unique-indices-method-similar-to-matlab/34446/3
        unique_mask::Array{Int64} = unique(z -> wl[z], 1:length(wl))
        wl = wl[unique_mask]
        fl = fl[unique_mask]

        # Bin data to required number of bins...

        # Log data first
        lfl = log10.(fl)
        lwl = log10.(wl)

        # Downsample spectrum onto N=tgt_bins, logwavelength space
        itp = extrapolate(interpolate((lwl,),lfl, Gridded(Linear())), Flat())
        bin_c = collect(range(start=lwl[1],  stop=lwl[end], length=tgt_bins))
        bin_v = zeros(Float64, tgt_bins)
        @. bin_v = itp(bin_c)

        # new arrrays
        owl = zeros(Float64, tgt_bins)
        ofl = zeros(Float64, tgt_bins)
        @. owl = 10.0 .^ bin_c
        @. ofl = 10.0 .^ bin_v

        # Convert units
        owl *= 1.0e-9  # [nm] -> [m]
        ofl *= 1.0e6   # [erg s-1 cm-2 nm-1] -> [W m-3]

        # Write file
        @debug "Writing stellar spectrum to SOCRATES format"
        len_new::Int64 = length(owl)
        open(star_file, "w") do f

            # Header
            write(f, "Star spectrum at TOA. Created by AGNI. \n")
            write(f, "      WAVELENGTH        IRRADIANCE\n")
            write(f, "          (m)               (W/m3)\n")
            write(f, "*BEGIN_DATA\n")

            # Body
            for i in 1:len_new
                @printf(f, "      %1.7e      %1.7e\n", owl[i],ofl[i])
            end

            # Footer
            write(f, "*END\n")
            write(f, " ")
        end

        return true
    end


    """
    **Insert a stellar spectrum, Rayleigh coeffs, and aerosol properties into a SOCRATES spectral file.**

    Will not overwrite the original file.

    Inserts stellar spectrum first, then Rayleigh coefficients, then aerosol properties.
    The aerosol .avg files must already exist. The stellar spectrum must already exist.

    Arguments:
    - `orig_file::String`        Path to original spectral file.
    - `star_file::String`        Path to file containing stellar spectrum in SOC format.
    - `outp_file::String`        Path to output spectral file.
    - `insert_rscatter::Bool`    Calculate Rayleigh scattering coefficients?
    - `insert_aerosol::Bool`     Insert aerosol parametrizations?
    - `aerosol_avg_files::Dict`  Paths to aerosol .avg files (if `insert_aerosol=true`)

    Returns:
    - `success::Bool`            function executed successfully
    """
    function insert_blocks(orig_file::String, star_file::String, outp_file::String,
                            insert_rscatter::Bool, insert_aerosol::Bool;
                            aerosol_avg_files::Dict{String,String}=Dict())::Bool

        # Inputs to prep_spec
        prep_spec = abspath(ENV["RAD_DIR"],"bin","prep_spec")
        @debug "Using prep_spec at: "*prep_spec
        star_inputs = [
            "6","n","T",            # ask prep_spec to tabulate the thermal source function
            "100 4000","1200",      # tmp_min, tmp_max, num_points
            "2","n",star_file,      # ask prep_spec to insert this stellar spectrum
            "y"                     # exit prep_spec
            ]

        # Check files exist
        if !isfile(orig_file)
            @warn "Original spectral file not found: '$orig_file'"
            return false
        end
        if !isfile(star_file)
            @warn "Stellar spectrum file not found: '$star_file'"
            return false
        end

        # Copy original file to output file
        cp(orig_file,      outp_file;      force=true)
        cp(orig_file*"_k", outp_file*"_k"; force=true)

        # Count number of gaseous absorbers in file, for later use
        num_gases = count_gases(outp_file)
        if num_gases < 0
            @warn "Failed to parse gaseous absorbers in spectral file"
            return false
        end

        # Add all aerosol header information into spectral file
        if insert_aerosol
            insert_aerosol_header(outp_file, [s for s in keys(aerosol_avg_files)]) || return false
        end

        # Check that aerosol .avg files exist
        for avgfile in values(aerosol_avg_files)
            if !isfile(avgfile)
                @warn "Aerosol .avg file not found: '$avgfile'"
                return false
            end
        end

        # Write executable
        execpath::String = tempname() * "_agni_insert_stellar.sh"
        @debug "Wrapping script: $execpath"
        rm(execpath, force=true)
        open(execpath, "w") do f

            # exec prep_spec
            write(f, prep_spec*" <<-EOF\n")

            # paths
            write(f, outp_file*" \n")
            write(f, "a \n") # modify output file in-place

            # write(f, "n \n") # write to new file
            # write(f, outp_file*" \n")

            # write thermal source function + stellar spectrum
            for inp in star_inputs
                write(f, inp*" \n")
            end

            # write rayleigh coefficients
            todo_str = "Inserting stellar spectrum"
            if insert_rscatter
                todo_str *= ", Rayleigh coefficients"
                write(f, "3 \n")       #  block 3, please
                write(f, "c \n")       #  custom composition
                write(f, "a \n")       #  all gases

                # If there's only one gas, add an extra 'y' confirmation
                if (num_gases == 1) && !startswith(get_socrates_version(), "24")
                    write(f, "y \n")
                end
            end

            # write aerosol properties
            if insert_aerosol
                todo_str *= ", aerosol properties"

                for avgfile in values(aerosol_avg_files)
                    write(f, "11 \n")
                    write(f, avgfile*" \n")
                    # write(f, "y \n")
                end
            end

            @info todo_str

            # exit prep_spec
            write(f, "-1 \n")
            write(f, "EOF\n")
            write(f, " ")
        end

        # Run executable
        @debug "Running prep_spec now"
        try
            ps = run(pipeline(`bash $execpath`, stdout=devnull))

            if !success(ps)
                @warn "prep_spec failed with exit code $(ps.exitcode)"
                @warn "Command: bash $execpath"
                return false
            end

        catch e
            @warn "Failed to run prep_spec: $e"
            @warn "Command: bash $execpath"
            return false
        end

        # Tidy up
        rm(execpath)
        return true
    end


    """
    **Generate SOCRATES aerosol `.avg` files from monochromatic aerosol scattering files.**

    Uses `scatter_average_90` as documented in the SOCRATES user guide with
    solar-weighted averaging (`-S <solar> -w`). Does not modify the spectral file.

    Arguments:
    - `orig_file::String`             Original spectral file
    - `species::Array{String,1}`      List of aerosol species
    - `output_dir::String`            Directory where output `.avg` files are written
    - `phase_moments::Int64`          Number of phase-function moments to retain
    - `star_file::String`             Solar spectrum file for SW-weighted averaging
    - `scattering_dir::String`        Directory containing scattering files (.mon)

    Returns:
    - `avg_files::Dict{String,String}` Mapping from aerosol species to generated `.avg` file paths.
    """
    function generate_aerosol_avg_files(orig_file::String,
                                        species::Array{String,1},
                                        output_dir::String,
                                        phase_moments::Int64,
                                        star_file::String,
                                        scattering_dir::String)::Dict{String,String}

        avg_files::Dict = Dict{String, String}()

        if !isfile(orig_file)
            @warn("Spectral file not found: '$orig_file'")
            return avg_files
        end
        if phase_moments < 1
            @warn("phase_moments must be >= 1, got $phase_moments")
            return avg_files
        end
        if isempty(star_file)
            @warn("star_file must be provided for scatter_average")
            return avg_files
        end

        if !isfile(star_file)
            @warn("Solar spectrum file not found: '$star_file'")
            return avg_files
        end

        scat_av_90 = abspath(ENV["RAD_DIR"], "bin", "scatter_average_90")

        # check mon files exist
        for s in species
            mon = abspath(joinpath(scattering_dir, s*".mon"))
            if !isfile(mon)
                @warn "Monochromatic aerosol data not found for '$s'"
                @warn "    Expected at: '$mon'"
                return avg_files
            end
        end

        # Write executable
        execpath::String = tempname() * "_agni_aerosol_scatter.sh"
        @debug "Wrapping script: $execpath"
        rm(execpath, force=true)
        open(execpath, "w") do f

            # loop over aerosol species
            for s in species
                @debug "    adding species $s"

                # exec scatter_average_90
                write(f, scat_av_90*" <<-EOF\n")

                # input spectral file
                write(f, orig_file*" \n")

                # input mon file
                write(f, abspath(joinpath(scattering_dir, s*".mon"))*" \n")

                # averaging parameters
                write(f, "3 \n") # spectrum weighting
                write(f, star_file*" \n") # solar spectrum for weighting
                write(f, "1 \n") # thin averaging (-w)
                write(f, "n \n") # no instrumental response required
                write(f, "y \n") # per block?

                # output avg file
                avg_files[s] = abspath(joinpath(output_dir, s*".avg"))
                rm(avg_files[s], force=true)
                write(f, avg_files[s]*" \n") # path to output file

                # number of phase moments
                write(f, "$phase_moments \n")
                write(f, "n \n") # not dependent on radius

                # exit scatter_average_90
                write(f, "EOF\n")
                write(f, " ")
            end
        end

        # Check that any species were processed
        if isempty(avg_files)
            @debug "No aerosol species processed."
            return avg_files
        end

        # Run executable
        @debug "Running scatter_average_90 now"
        try
            ps = run(pipeline(`bash $execpath`, stdout=devnull))

            if !success(ps)
                @warn "scatter_average_90 failed with exit code $(ps.exitcode)"
                @warn "Command: bash $execpath"
                return avg_files
            end

        catch e
            @warn "Failed to run scatter_average_90: $e"
            @warn "Command: bash $execpath"
            return avg_files
        end

        # Tidy up
        rm(execpath)
        return avg_files
    end

end # end module spectrum
