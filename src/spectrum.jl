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
    using PCHIPInterpolation
    using DelimitedFiles
    using LinearAlgebra

    """
    **Load stellar spectrum from a text file.**

    The flux needs to be scaled to the top of the atmosphere.

    Arguments:
    - `path::String`        Path to the file

    Returns:
    - `wl::Array`           Wavelength array [nm]
    - `fl::Array`           Flux array [erg s-1 cm-2 nm-1]
    """
    function load_from_file(path::String)::Tuple{Array{Float64,1}, Array{Float64,1}}

        if isfile(path)
            spec_data = readdlm(abspath(path), '\t', Float64; header=false, skipstart=2)
        else
            error("Cannot find stellar spectrum at '$path'")
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
    - `nbins_max::Int`      Maximum number of points in the spectrum
    """
    function write_to_socrates_format(wl::Array{Float64,1}, fl::Array{Float64,1}, star_file::String, nbins_max::Int=99900)

        len_wl::Int = length(wl)
        len_fl::Int = length(fl)
        socrates_nbins_max::Int = Int(1e5 - 3)  # do not change this
        tgt_bins::Int = min(nbins_max, len_wl, socrates_nbins_max)

        # Validate
        if len_wl != len_fl
            error("Stellar wavelength and flux arrays have different lengths")
        end 
        if len_wl < 500
            @warn "Loaded stellar spectrum is very short!"
        end
        if minimum(wl) < 1.0e-45
            error("Minimum wavelength is too small")
        end
        clamp!(fl, 1.0e-45, 1.0e+45)  # Clamp values

        # Bin data to required number of bins...

        # Log data first 
        lfl = log10.(fl)
        lwl = log10.(wl)

        # Downsample spectrum onto N=tgt_bins, logwavelength space
        itp = Interpolator(lwl, lfl)
        bin_c = collect(range(start=lwl[1],  stop=lwl[end], length=tgt_bins))
        bin_v = zeros(Float64, tgt_bins)
        bin_v[:] .= itp.(bin_c[:])

        owl = 10.0 .^ Array(bin_c[:])
        ofl = 10.0 .^ Array(bin_v[:])

        # Convert units
        owl = owl .* 1.0e-9  # [nm] -> [m]
        ofl = ofl .* 1.0e6   # [erg s-1 cm-2 nm-1] -> [W m-3]

        # Write file
        len_new::Int = length(owl)
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

        return nothing
    end


    """
    **Insert a stellar spectrum and Rayleigh coeffs into a SOCRATES spectral file.**

    Will not overwrite the original file. 

    Arguments:
    - `orig_file::String`        Path to original spectral file.
    - `star_file::String`        Path to file containing stellar spectrum written with this module.
    - `outp_file::String`        Path to output spectral file.
    - `insert_rscatter::Bool`    Calculate Rayleigh scattering coefficients?
    """
    function insert_stellar_and_rscatter(orig_file::String, star_file::String, outp_file::String, insert_rscatter::Bool)

        # Inputs to prep_spec
        prep_spec = joinpath(ENV["RAD_DIR"],"bin","prep_spec")
        @debug "Using prep_spec at: "*prep_spec
        star_inputs = [
            "6","n","T",            # ask prep_spec to tabulate the thermal source function
            "100 4000","1200",      # tmp_min, tmp_max, num_points
            "2","n",star_file,      # ask prep_spec to insert this stellar spectrum
            "y"                     # exit prep_spec
            ]

        # Write executable 
        execpath::String = "/tmp/$(abs(rand(Int,1)[1]))_agni_insert_stellar.sh"
        rm(execpath, force=true)
        open(execpath, "w") do f
            
            # exec prep_spec 
            write(f, prep_spec*" <<-EOF\n")

            # paths
            write(f, orig_file*" \n")
            write(f, "n \n")
            write(f, outp_file*" \n")

            # write thermal source function + stellar spectrum
            for inp in star_inputs 
                write(f, inp*" \n")
            end 

            # write rayleigh coefficients
            if insert_rscatter
                @info "Inserting stellar spectrum and Rayleigh coefficients" 
                write(f, "3 \n")       #  block 3, please
                write(f, "c \n")       #  custom composition
                write(f, "a \n")       #  all gases
            else 
                @info "Inserting stellar spectrum"
            end 

            # exit prep_spec 
            write(f, "-1 \n")
            write(f, "EOF\n")
            write(f, " ")
        end

        # Run executable
        @debug "Running prep_spec now"
        run(pipeline(`bash $execpath`, stdout=devnull))

        # Delete executable 
        rm(execpath)

        return nothing 
    end 

end # end module spectrum
