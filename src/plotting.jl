# Contains functions for plotting/diagnosing the atmosphere 

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module plotting

    # Import stuff
    using Plots
    using LaTeXStrings
    using Glob

    import atmosphere

    """
    Plot the temperature-pressure profile.
    """
    function plot_pt(atmos, fname)
        
        # Interleave cell-centre and cell-edge arrays
        arr_P, arr_T = atmosphere.get_interleaved_pt(atmos)
        arr_P *= 1e-5

        ylims  = (arr_P[1], arr_P[end])
        yticks = 10 .^ floor.(range( log10(ylims[1]), stop=log10(ylims[2]), step=1))

        # Create plot
        plt = plot(framestyle=:box, ylims=ylims, yticks=yticks)

        # Plot temperature
        scatter!(plt, [atmos.tstar], [atmos.pl[end]*1e-5], color="brown3", label=L"T_*")
        plot!(plt, arr_T, arr_P, lc="black", lw=2, label=L"T(p)")
        xlabel!(plt, "Temperature [K]")
        ylabel!(plt, "Pressure [bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)
                
        savefig(plt, fname)

    end

    """
    Plot the current temperature-pressure profile, current heating rates, and
    optionally the previous states that the atmosphere has taken.
    """
    function plot_solver(atmos, fname; hist_tmpl::Array=[])

        dpi=250
        lw=1.5

        # Interleave cell-centre and cell-edge arrays for current atmosphere
        arr_P, arr_T = atmosphere.get_interleaved_pt(atmos)
        arr_P *= 1e-5

        ylims  = (arr_P[1], arr_P[end])
        yticks = 10 .^ floor.(range( log10(ylims[1]), stop=log10(ylims[2]), step=1))

        # Optionally get past atmosphere
        plot_hist = !isempty(hist_tmpl)
        if plot_hist 
            len_hist = size(hist_tmpl, 1)-1
        end

        # Create plot 1
        plt1 = plot(framestyle=:box, legend=:topright, ylims=ylims, yticks=yticks)

        # Plot temperature
        scatter!(plt1, [atmos.tstar], [atmos.pl[end]*1e-5], color="brown3", label=L"T_*") 
        plot!(plt1, arr_T, arr_P, lc="black", lw=lw, label=L"T_n")
        if plot_hist 
            for i in range(len_hist-1,1,step=-1)
                xvals = hist_tmpl[i,1:end]
                if count(x->x!=0.0, xvals) > 0
                    alpha = Float64(i)/(len_hist)
                    n = len_hist-i
                    plot!(plt1, xvals, atmos.pl*1e-5, lc="black", linealpha=alpha, lw=lw, label=L"T_{n-%$n}")
                end
            end 
        end
        

        xlabel!(plt1, "Temperature [K]")
        ylabel!(plt1, "Pressure [bar]")
        yflip!(plt1)
        yaxis!(plt1, yscale=:log10)

        # Create plot 2
        plt2 = plot(framestyle=:box, legend=:topleft, ylims=ylims, yticks=yticks)

        # Plot heating rate
        p =  atmos.p*1e-5
        abshr = zeros(Float64, atmos.nlev_c)
        poshr = trues(atmos.nlev_c)
        for i in 1:atmos.nlev_c
            abshr[i] = abs(atmos.heating_rate[i])
            poshr[i] = (atmos.heating_rate[i] >= 0)
        end 
        plot!(plt2, abshr, p, lc="brown3", lw=lw, label=L"|H_n|")
        scatter!(plt2, abshr[poshr],    p[poshr],    markershape=:diamond, markeralpha=0.8, label=L"H_n>0")
        scatter!(plt2, abshr[.!poshr],  p[.!poshr],  markershape=:circle,  markeralpha=0.8, label=L"H_n<0")
        xlabel!(plt2, "Heating rate [K/day]")
        yflip!(plt2)
        yaxis!(plt2, yscale=:log10)
        xaxis!(plt2, xscale=:log10, xlims=(1e-2, maximum(abshr)))
        
        # Combine subplots and save
        plt = plot(plt1, plt2, layout=(1,2), dpi=dpi)
        savefig(plt, fname)

    end

    """
    Plot the fluxes at each pressure level
    """
    function plot_fluxes(atmos, fname)

        arr_P = atmos.pl .* 1.0e-5 # Convert Pa to bar

        ylims  = (arr_P[1], arr_P[end])
        yticks = 10 .^ floor.(range( log10(ylims[1]), stop=log10(ylims[2]), step=1))

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

    """
    Animate output frames from solver, using ffmpeg
    """
    function anim_solver(atmos)

        # Command line format:
        # bash> ffmpeg -framerate 16 -i out/solver_monitor_%04d.png -y out/anim.mp4

        # Find output files
        out = atmos.OUT_DIR
        frames = glob("solver_monitor_*.png",out)
        nframes = length(frames)
        if nframes < 1
            println("WARNING: Cannot animate solver because no output frames were found")
        end

        # Create animation
        runtime = 15.0
        fps = max(nframes/runtime, 5)
        run(`ffmpeg -loglevel quiet -framerate $fps -i $out/solver_monitor_%04d.png -y $out/anim.mp4`)

    end


end

