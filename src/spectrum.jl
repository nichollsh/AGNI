# Contains the spectrum module, for updating the SOCRATES spectral file at 
# runtime with the solar flux / thermal source function, and scattering.

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 


module spectrum 

    using Revise
    using Printf
    using BinnedStatistics
    using DelimitedFiles
    using LinearAlgebra

    """
    **Load stellar spectrum from a text file.**

    The flux needs to be provided at 1 AU from the star's surface.

    Arguments:
    - `path::String`        Path to the file

    Returns:
    - `wl::Array`           Wavelength array [nm]
    - `fl::Array`           Flux array [erg s-1 cm-2 nm-1]
    """
    function load_from_file(path::String)::Tuple{Array, Array}

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
    function write_to_socrates_format(wl::Array, fl::Array, star_file::String, nbins_max::Int=99900)

        len_wl = length(wl)
        len_fl = length(fl)
        socrates_nbins_max = Int(1e5 - 3)  # do not change this
        tgt_bins = min(nbins_max, len_wl, socrates_nbins_max)

        # Validate
        if len_wl != len_fl
            error("Stellar wavelength and flux arrays have different lengths")
        end 
        if len_wl < 500
            println("WARNING: Loaded stellar spectrum is very short!")
        end
        if minimum(wl) < 1.0e-45
            error("Minimum wavelength is too small")
        end
        clamp!(fl, 1.0e-45, 1.0e+45)  # Clamp values

        # Bin data if required 
        if tgt_bins < len_wl

            # Log data first 
            lfl = log10.(fl)
            lwl = log10.(wl)

            # Fit histogram
            bin_e, bin_c, bin_v = binnedStatistic(lwl,lfl,nbins=tgt_bins,statistic=:mean)

            owl = 10 .^ Array(bin_c[:])
            ofl = 10 .^ Array(bin_v[:])

        # No binning required
        else 
            owl = wl 
            ofl = fl
        end 
        
        # Convert units
        owl = owl .* 1.0e-9  # [nm] -> [m]
        ofl = ofl .* 1.0e6   # [erg s-1 cm-2 nm-1] -> [W m-3]

        # Write file
        println("Opening file $star_file")
        open(star_file, "w") do f

            # Header
            write(f, "Star spectrum at 1 AU. Created by AGNI. \n")
            write(f, "      WAVELENGTH        IRRADIANCE\n")
            write(f, "          (m)               (W/m3)\n")
            write(f, "*BEGIN_DATA\n")
            
            # Body
            for i in 1:len_wl
                @printf(f, "      %1.7e      %1.7e\n", owl[i],ofl[i])
            end

            # Footer 
            write(f, "*END\n")
            write(f, " ")
        end

        return nothing
    end


    """
    **Insert a stellar spectrum into a SOCRATES spectral file.**

    Will not overwrite the original file.

    Arguments:
    - `orig_file::String`        Path to original spectral file.
    - `star_file::String`        Path to file containing stellar spectrum written with this module.
    - `outp_file::String`        Path to output spectral file.
    """
    function insert_stellar_spectrum(orig_file::String, star_file::String, outp_file::String)
        # k files
        orig_filek::String = orig_file*"_k"
        outp_filek::String = outp_file*"_k"

        # Delete "new" files if they already exist
        rm(outp_file , force=true)
        rm(outp_filek, force=true)

        # Copy original files to new location (retain old files)
        cp(orig_file,  outp_file )
        cp(orig_filek, outp_filek)

        # prep_spec inputs
        inputs = [outp_file,"a","6","n","T","100 4000","1000","2","n",star_file,"y","-1","EOF"]

        # Write executable 
        execpath::String = "/tmp/$(rand(Int,1)[1])_agni_insert_stellar.sh"
        rm(execpath, force=true)
        open(execpath, "w") do f
            write(f, "socrates/bin/prep_spec <<EOF\n")
            for inp in inputs 
                write(f, inp*" \n")
            end 
        end

        # Run executable
        println("Inserting stellar spectrum")
        run(pipeline(`bash $execpath`, stdout=devnull))

        # Delete executable 
        rm(execpath)

        return nothing 
    end 

end # end module spectrum
