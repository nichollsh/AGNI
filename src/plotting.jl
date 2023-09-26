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
    using Revise

    import atmosphere

    """
    Plot the temperature-pressure profile.
    """
    function plot_pt(atmos, fname)
        
        # Interleave cell-centre and cell-edge arrays
        arr_P, arr_T = atmosphere.get_interleaved_pt(atmos)
        arr_P *= 1e-5

        ylims  = (arr_P[1]/1.5, arr_P[end]*1.5)
        yticks = 10.0 .^ round.(Int,range( log10(ylims[1]), stop=log10(ylims[2]), step=1))

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

        ylims  = (arr_P[1]/1.5, arr_P[end]*1.5)
        yticks = 10.0 .^ round.(Int,range( log10(ylims[1]), stop=log10(ylims[2]), step=1))

        # Optionally get past atmosphere
        plot_hist = !isempty(hist_tmpl)
        if plot_hist 
            len_hist = size(hist_tmpl, 1)-1
        end

        # Create plot 1
        plt1 = plot(framestyle=:box, legend=:topright, ylims=ylims, yticks=yticks)

        # Plot surface temperature(s)
        scatter!(plt1, [atmos.tmp_magma], [atmos.pl[end]*1e-5], color="cornflowerblue", label=L"T_m") 
        scatter!(plt1, [atmos.tstar],     [atmos.pl[end]*1e-5], color="brown3",         label=L"T_*") 
        
        # Plot temperature profiles 
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
        
        # Highlight convective regions 
        atmosphere.dryconvection_check!(atmos, tmp_eps=10.0)
        region_p = []
        region_T = []
        region_new = true
        for i in 1:atmos.nlev_c
            if atmos.conv_inst[i]
                if region_new
                    region_new = false
                    if length(region_p) > 0
                        plot!(plt1, region_T, region_p*1e-5, color="orchid2", linealpha=0.4, lw=4*lw, label="")
                    end
                    region_p = []
                    region_T = []
                end 
                push!(region_p, atmos.p[i])
                push!(region_p, atmos.pl[i])
                push!(region_T, atmos.tmp[i])
                push!(region_T, atmos.tmpl[i])
            else 
                region_new = true
            end
        end 

        xlabel!(plt1, "Temperature [K]")
        ylabel!(plt1, "Pressure [bar]")
        yflip!(plt1)
        yaxis!(plt1, yscale=:log10)

        # Process heating rates 
        p =  atmos.p*1e-5
        abshr = zeros(Float64, atmos.nlev_c)
        poshr = trues(atmos.nlev_c)
        for i in 1:atmos.nlev_c
            abshr[i] = abs(atmos.heating_rate[i])
            poshr[i] = (atmos.heating_rate[i] >= 0)
        end 

        xlims  = (1e-2, maximum(abshr))
        xticks = 10.0 .^ round.(Int,range( log10(xlims[1]), stop=log10(xlims[2]), step=1))

        # Create plot 2
        plt2 = plot(framestyle=:box, legend=:topleft, ylims=ylims, yticks=yticks, xlims=xlims, xticks=xticks)

        # Plot heating rate
        plot!(plt2, abshr, p, lc="brown3", lw=lw, label=L"|H_{n-1}|")
        scatter!(plt2, abshr[poshr],    p[poshr],    markershape=:diamond, markeralpha=0.8, label=L">0")
        scatter!(plt2, abshr[.!poshr],  p[.!poshr],  markershape=:circle,  markeralpha=0.8, label=L"<0")
        xlabel!(plt2, "Heating rate [K/day]")
        yflip!(plt2)
        yaxis!(plt2, yscale=:log10)
        xaxis!(plt2, xscale=:log10)
        
        # Combine subplots and save
        plt = plot(plt1, plt2, layout=(1,2), dpi=dpi)
        savefig(plt, fname)

    end

    """
    Plot the fluxes at each pressure level
    """
    function plot_fluxes(atmos, fname)

        arr_P = atmos.pl .* 1.0e-5 # Convert Pa to bar

        ylims  = (arr_P[1]/1.5, arr_P[end]*1.5)
        yticks = 10.0 .^ round.(Int,range( log10(ylims[1]), stop=log10(ylims[2]), step=1))

        w = 2
        plt = plot(legend=:topleft, framestyle=:box, ylims=ylims, yticks=yticks)

        col_u = "brown3"
        col_d = "seagreen"
        col_n = "black"

        max_fl = 1e2
        
        # LW component
        if atmos.is_out_lw
            y = abs.(atmos.flux_d_lw)
            max_fl = max(max_fl, maximum(y))
            plot!(plt, y, arr_P, label="RAD DN LW", lw=w, lc=col_d, ls=:dash)

            y = abs.(atmos.flux_u_lw)
            max_fl = max(max_fl, maximum(y))
            plot!(plt, y, arr_P, label="RAD UP LW", lw=w, lc=col_u, ls=:dash)
        end

        # SW component
        if atmos.is_out_sw
            y = abs.(atmos.flux_d_sw)
            max_fl = max(max_fl, maximum(y))
            plot!(plt, y, arr_P, label="RAD DN SW", lw=w, lc=col_d, ls=:dot)

            y = abs.(atmos.flux_u_sw)
            max_fl = max(max_fl, maximum(y))
            plot!(plt, y, arr_P, label="RAD UP SW", lw=w, lc=col_u, ls=:dot)
        end 

        # Net fluxes
        if atmos.is_out_lw && atmos.is_out_sw
            y = abs.(atmos.flux_u)
            max_fl = max(max_fl, maximum(y))
            plot!(plt, y, arr_P,    label="RAD UP",  lw=w, lc=col_u, ls=:solid)

            y = abs.(atmos.flux_d)
            max_fl = max(max_fl, maximum(y))
            plot!(plt, y, arr_P,    label="RAD DN",  lw=w, lc=col_d, ls=:solid)

            absnet = zeros(Float64, atmos.nlev_l)
            posnet = trues(atmos.nlev_l)
            for i in 1:atmos.nlev_l
                absnet[i] = abs(atmos.flux_n[i])
                posnet[i] = (atmos.flux_n[i] >= 0)
            end 
            max_fl = max(max_fl, maximum(absnet))
            plot!(plt, absnet, arr_P, label="RAD NET", lw=w, lc=col_n, ls=:solid)
            scatter!(plt, absnet[  posnet], arr_P[  posnet],  markershape=:diamond, markeralpha=0.8, label=L">0")
            scatter!(plt, absnet[.!posnet], arr_P[.!posnet],  markershape=:circle,  markeralpha=0.8, label=L"<0")
        end 

        # Sensible heat
        if atmos.flux_sens > 0
            scatter!(plt, [atmos.flux_sens],      [arr_P[end]], markershape=:utriangle, markercolor=col_u, label="SENS")
        else
            scatter!(plt, [abs(atmos.flux_sens)], [arr_P[end]], markershape=:dtriangle, markercolor=col_d, label="SENS")
        end 

        xlims  = (1e-1, max_fl * 1.5)
        xticks = 10.0 .^ round.(Int,range( log10(xlims[1]), stop=log10(xlims[2]), step=1))

        xlabel!(plt, "Unsigned flux [W m-2]")
        ylabel!(plt, "Pressure [bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)
        xaxis!(plt, xscale=:log10, xlims=xlims, xticks=xticks)

        savefig(plt, fname)

    end

    """
    Animate output frames from solver, using ffmpeg
    """
    function anim_solver(atmos)

        # Command line format:
        # bash> ffmpeg -framerate 16 -i out/solver_monitor_%04d.png -y out/anim.mp4

        runtime = 15.0 # seconds

        # Find output files
        out = atmos.OUT_DIR
        frames = glob("solver_monitor_*.png",out)
        nframes = length(frames)

        # Create animation
        if nframes < 1
            println("WARNING: Cannot animate solver because no output frames were found")
        else 
            fps = max(nframes/runtime, 5)
            run(`ffmpeg -loglevel quiet -framerate $fps -i $out/solver_monitor_%04d.png -y $out/anim.mp4`)
        end

    end


end

