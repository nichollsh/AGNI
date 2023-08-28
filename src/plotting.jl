# Contains functions for plotting/diagnosing the atmosphere 

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module plotting

    # Import stuff
    using Plots
    include("atmosphere.jl")

    function plot_pt(atmos, output_dir)
        """
        Plot the temperature-pressure profile.
        """

        println("Plotting PT profile")

        # Interleave cell-centre and cell-edge arrays
        arr_P, arr_T = atmosphere.get_interleaved_pt(atmos)
        arr_P *= 1e-5

        # Plot PT profile
        plt = plot(legend=false)
        plot!(plt, arr_T, arr_P, lc="black", lw=2)
        xlabel!(plt, "Temperature [K]")
        ylabel!(plt, "Pressure [bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)
                
        savefig(plt, joinpath(output_dir,"plot_pt.pdf"))

    end

    function plot_fluxes(atmos, output_dir)
        """
        Plot the fluxes at each pressure level
        """

        println("Plotting fluxes")

        arr_P = atmos.pl .* 1.0e-5 # Convert Pa to bar

        
        w = 2
        plt = plot(legend=:top)
        
        if atmos.is_out_lw
            c = "brown3"
            plot!(plt, -1.0.*atmos.flux_d_lw, arr_P, label="LW d", lw=w, lc=c, ls=:dot)
            plot!(plt, atmos.flux_u_lw, arr_P,       label="LW u", lw=w, lc=c, ls=:dash)
            plot!(plt, atmos.flux_n_lw, arr_P,       label="LW n", lw=w, lc=c, ls=:solid)
        end

        if atmos.is_out_sw
            c = "royalblue3"
            plot!(plt, -1.0.*atmos.flux_d_sw, arr_P, label="SW d", lw=w, lc=c, ls=:dot)
            plot!(plt, atmos.flux_u_sw, arr_P,       label="SW u", lw=w, lc=c, ls=:dash)
            plot!(plt, atmos.flux_n_sw, arr_P,       label="SW n", lw=w, lc=c, ls=:solid)
        end 

        if atmos.is_out_lw && atmos.is_out_sw
            c = "seagreen"
            plot!(plt, atmos.flux_n, arr_P, label="NET", lw=w, lc=c, ls=:solid)
        end 

        xlabel!(plt, "Upward directed flux [W m-2]")
        ylabel!(plt, "Pressure [bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)

        savefig(plt, joinpath(output_dir,"plot_fluxes.pdf"))

    end
    

end

