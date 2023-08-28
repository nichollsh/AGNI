# Contains functions for plotting/diagnosing the atmosphere 

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module plotting

    # Import stuff
    using Plots

    function plot_pt(atmos, output_dir)
        """
        Plot the temperature-pressure profile.
        """

        println("Plotting PT profile")

        arr_T = atmos.tmpl
        arr_P = atmos.pl .* 1e-5

        plt = plot(arr_T, arr_P, 
                legend=false, lc="black", lw=2)
        xlabel!(plt, "Temperature [K]")
        ylabel!(plt, "Pressure [Bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)
                
        savefig(joinpath(output_dir,"plot_pt.pdf"))

    end

    function plot_fluxes(atmos, output_dir)
        """
        Plot the fluxes at each pressure level
        """

        println("Plotting fluxes")

        arr_P = atmos.pl .* 1.0e-5 # Convert Pa to bar

        c = "red"
        w = 2
        plt = plot()
        plot!(plt, -1.0.*atmos.flux_d_lw, arr_P, label="LW d", lw=w, lc=c, ls=:dot)
        plot!(plt, atmos.flux_u_lw, arr_P, label="LW u", lw=w, lc=c, ls=:dash)
        plot!(plt, atmos.flux_n_lw, arr_P, label="LW n", lw=w, lc=c, ls=:solid)

        c = "blue"
        plot!(plt, -1.0.*atmos.flux_d_sw, arr_P, label="SW d", lw=w, lc=c, ls=:dot)
        plot!(plt, atmos.flux_u_sw, arr_P, label="SW u", lw=w, lc=c, ls=:dash)
        plot!(plt, atmos.flux_n_sw, arr_P, label="SW n", lw=w, lc=c, ls=:solid)

        xlabel!(plt, "Upward directed flux [W m-2]")
        ylabel!(plt, "Pressure [Bar]")
        xlabel!(plt, "Temperature [K]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)

        savefig(plt, joinpath(output_dir,"plot_fluxes.pdf"))

    end
    

end

