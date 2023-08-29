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


    function plot_pt(atmos, fname)
        """
        Plot the temperature-pressure profile.
        """

        # Interleave cell-centre and cell-edge arrays
        arr_P, arr_T = atmosphere.get_interleaved_pt(atmos)
        arr_P *= 1e-5

        # Plot PT profile
        plt = plot(legend=false, framestyle=:box)
        plot!(plt, arr_T, arr_P, lc="black", lw=2)
        xlabel!(plt, "Temperature [K]")
        ylabel!(plt, "Pressure [bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)
                
        savefig(plt, fname)

    end

    function plot_fluxes(atmos, fname)
        """
        Plot the fluxes at each pressure level
        """

        arr_P = atmos.pl .* 1.0e-5 # Convert Pa to bar

        w = 2
        plt = plot(legend=:topleft, framestyle=:box)

        col_u = "brown3"
        col_d = "seagreen"
        col_n = "black"
        
        if atmos.is_out_lw
            plot!(plt, -1.0.*atmos.flux_d_lw, arr_P, label="DN LW", lw=w, lc=col_d, ls=:dash)
            plot!(plt, atmos.flux_u_lw, arr_P,       label="UP LW", lw=w, lc=col_u, ls=:dash)
        end

        if atmos.is_out_sw
            plot!(plt, -1.0.*atmos.flux_d_sw, arr_P, label="DN SW", lw=w, lc=col_d, ls=:dot)
            plot!(plt, atmos.flux_u_sw, arr_P,       label="UP SW", lw=w, lc=col_u, ls=:dot)
        end 

        if atmos.is_out_lw && atmos.is_out_sw
            plot!(plt, atmos.flux_u, arr_P,         label="UP",  lw=w, lc=col_u, ls=:solid)
            plot!(plt, -1.0*atmos.flux_d, arr_P,    label="DN",  lw=w, lc=col_d, ls=:solid)
            plot!(plt, atmos.flux_n, arr_P,         label="NET", lw=w, lc=col_n, ls=:solid)
        end 

        xlabel!(plt, "Upward directed flux [W m-2]")
        ylabel!(plt, "Pressure [bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)

        savefig(plt, fname)

    end

end

