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
    import phys

    """
    Plot the temperature-pressure profile.
    """
    function plot_pt(atmos, fname; dpi::Int=250, incl_magma::Bool=false)
        
        # Interleave cell-centre and cell-edge arrays
        arr_P, arr_T = atmosphere.get_interleaved_pt(atmos)
        arr_P *= 1e-5

        ylims  = (arr_P[1]/1.5, arr_P[end]*1.5)
        yticks = 10.0 .^ round.(Int,range( log10(ylims[1]), stop=log10(ylims[2]), step=1))

        # Create plot
        plt = plot(framestyle=:box, ylims=ylims, yticks=yticks, legend=:outertopright, dpi=dpi, size=(500,400), guidefontsize=9)

        # Plot temperature
        if incl_magma
            scatter!(plt, [atmos.tmp_magma], [atmos.pl[end]*1e-5], color="cornflowerblue", label=L"T_m") 
        end
        scatter!(plt, [atmos.tstar], [atmos.pl[end]*1e-5], color="brown3", label=L"T_*")
        plot!(plt, arr_T, arr_P, lc="black", lw=2, label=L"T(p)")
        xlabel!(plt, "Temperature [K]")
        ylabel!(plt, "Pressure [bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)
                
        savefig(plt, fname)
        return nothing
    end

    """
    Plot the composition of the atmosphere at each cell-centre location.
    """
    function plot_x(atmos, fname)
        
        arr_P = atmos.p .* 1.0e-5 # Convert Pa to bar
        ylims  = (arr_P[1]/1.5, arr_P[end]*1.5)
        yticks = 10.0 .^ round.(Int,range( log10(ylims[1]), stop=log10(ylims[2]), step=1))

        
        # Create plot
        plt = plot(framestyle=:box, ylims=ylims, yticks=yticks, legend=:outertopright, size=(500,400), guidefontsize=9)

        # Plot mole fractions for each gas
        min_x = 1.0e-3
        for (i_gas,gas) in enumerate(atmos.gases)  
            arr_x = atmos.layer_x[1:end,i_gas]
            if minimum(arr_x) > 0.0
                min_x = min(min_x, minimum(arr_x))
            end
            plot!(arr_x, arr_P, label=gas, lw=2)
        end

        xlims  = (max(min_x, 1.0e-20)*0.5, 1.2)
        xticks = 10.0 .^ round.(Int,range( log10(xlims[1]), stop=0, step=1))

        # Set figure properties
        xlabel!(plt, "Mole fraction")
        xaxis!(plt, xscale=:log10, xlims=xlims, xticks=xticks)

        ylabel!(plt, "Pressure [bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)
                
        savefig(plt, fname)
    end

    """
    Plot the current temperature-pressure profile, current heating rates, and
    optionally the previous states that the atmosphere has taken.
    """
    function plot_solver(atmos, fname; hist_tmpl::Array=[], incl_magma::Bool=false, dpi::Int=250, step::Int=-1)

        lw=1.5
        
        # Header info
        title = ""
        if step > 0
            title = "Step $step"
        end

        # Interleave cell-centre and cell-edge arrays for current atmosphere
        arr_P, arr_T = atmosphere.get_interleaved_pt(atmos)
        arr_P *= 1e-5

        ylims  = (arr_P[1]/1.5, arr_P[end]*1.5)
        yticks = 10.0 .^ round.(Int,range( log10(ylims[1]), stop=log10(ylims[2]), step=1))

        # Optionally get past atmosphere
        plot_hist = !isempty(hist_tmpl)
        if plot_hist 
            len_hist = size(hist_tmpl, 1)-1
            len_hist = min(len_hist, 4)  # don't plot more than 4 previous iterations
        end

        # Create plot 1
        plt1 = plot(framestyle=:box, legend=:topright, ylims=ylims, yticks=yticks, title=title, titlefontsize=9)

        # Plot surface temperature(s)
        if incl_magma
            scatter!(plt1, [atmos.tmp_magma], [atmos.pl[end]*1e-5], color="cornflowerblue", label=L"T_m") 
        end
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

        # Plot convection mask 
        convective_p = []
        convective_t = []
        for i in 1:atmos.nlev_c
            if atmos.mask_c[i] > 0
                append!(convective_p, atmos.p[i]*1e-5)
                append!(convective_t, atmos.tmp[i])
            end 
        end
        if length(convective_p) > 0
            scatter!(plt1, convective_t, convective_p, color="goldenrod2", label="Cnvct", markersize=2, markeralpha=0.8) 
        end 

        # Plot phase change mask 
        pchange_p = []
        pchange_t = []
        for i in 1:atmos.nlev_c
            if atmos.mask_p[i] > 0
                append!(pchange_p, atmos.p[i]*1e-5)
                append!(pchange_t, atmos.tmp[i])
            end 
        end
        if length(pchange_p) > 0
            scatter!(plt1, pchange_t, pchange_p, color="dodgerblue", label="Phase", markersize=2, markeralpha=0.8) 
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

        xlims  = (1e-5, maximum(abshr))
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

        return nothing
    end

    """
    Plot the fluxes at each pressure level
    """
    function plot_fluxes(atmos, fname; dpi::Int=250)

        arr_P = atmos.pl .* 1.0e-5 # Convert Pa to bar

        ylims  = (arr_P[1]/1.5, arr_P[end]*1.5)
        yticks = 10.0 .^ round.(Int,range( log10(ylims[1]), stop=log10(ylims[2]), step=1))

        w = 2.8
        plt = plot(legend=:outertopright, framestyle=:box, ylims=ylims, yticks=yticks, dpi=dpi, guidefontsize=9)

        col_u = "brown3"
        col_d = "seagreen"
        col_np = "deepskyblue2"
        col_nn = "lightslateblue"
        col_c = "goldenrod2"
        col_tp = "black"
        col_tn = "grey18"

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

        # Net radiative fluxes
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

            plot!(plt, absnet[  posnet], arr_P[  posnet], label="RAD NET"*L">0", lw=w*2.0, lc=col_np, ls=:solid)
            plot!(plt, absnet[.!posnet], arr_P[.!posnet], label="RAD NET"*L"<0", lw=w    , lc=col_nn, ls=:solid)
        end 

        # Convective flux (MLT)
        if any(x->x!=0.0, atmos.flux_c)
            plot!(plt, atmos.flux_c, arr_P, label="CONVECT", lw=w*1.2, lc=col_c, ls=:solid)
            max_fl = max(max_fl, maximum(atmos.flux_c))
        end 

        # Sensible heat
        if atmos.flux_sens != 0.0
            if atmos.flux_sens > 0.0
                scatter!(plt, [atmos.flux_sens],      [arr_P[end]], markershape=:utriangle, markercolor=col_u, label="SENSIBLE")
            else
                scatter!(plt, [abs(atmos.flux_sens)], [arr_P[end]], markershape=:dtriangle, markercolor=col_d, label="SENSIBLE")
            end
        end

        # Total flux
        abstot = zeros(Float64, atmos.nlev_l)
        postot = trues(atmos.nlev_l)
        for i in 1:atmos.nlev_l
            abstot[i] = abs(atmos.flux_tot[i])
            postot[i] = (atmos.flux_tot[i] >= 0)
        end 
        max_fl = max(max_fl, maximum(abstot))
        plot!(plt, abstot[  postot], arr_P[  postot], label="TOTAL"*L">0", lw=w*0.7, lc=col_tp, ls=:solid)
        plot!(plt, abstot[.!postot], arr_P[.!postot], label="TOTAL"*L"<0", lw=w*0.4, lc=col_tn, ls=:solid)

        # Set limits
        xlims  = (1e-1, max_fl * 1.5)
        xticks = 10.0 .^ round.(Int,range( log10(xlims[1]), stop=log10(xlims[2]), step=1))

        # Overplot convection and condensation mask
        for i in 1:atmos.nlev_c
            if atmos.mask_c[i] > 0
                plot!(plt, [xlims[1],xlims[2]], [atmos.p[i]/1.0e5, atmos.p[i]/1e5], opacity=0.2, linewidth=3.5, color="goldenrod2", label="")
            end 
            if atmos.mask_p[i] > 0
                plot!(plt, [xlims[1],xlims[2]], [atmos.p[i]/1.0e5, atmos.p[i]/1e5], opacity=0.2, linewidth=3.5, color="dodgerblue", label="")
            end 
        end

        xlabel!(plt, "Unsigned flux [W m-2]")
        ylabel!(plt, "Pressure [bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)
        xaxis!(plt, xscale=:log10, xlims=xlims, xticks=xticks)

        savefig(plt, fname)

        return nothing
    end

    """
    Plot emission spectrum at the TOA
    """
    function plot_emission(atmos, fname; dpi::Int=250, planck_tmp::Float64=0.0)

        # Check that we have data 
        if !(atmos.is_out_lw && atmos.is_out_sw)
            error("Cannot plot emission spectrum because radiances have not been calculated")
        end

        # Get emission spectrum data
        xe = zeros(Float64, atmos.nbands)
        ye = zeros(Float64, atmos.nbands)
        for ba in 1:atmos.nbands
            # x value - band centres [nm]
            xe[ba] = 0.5 * (atmos.bands_min[ba] + atmos.bands_max[ba]) * 1.0e9
            
            # y value - spectral flux [erg s-1 cm-2 nm-1]
            w  = (atmos.bands_max[ba] - atmos.bands_min[ba]) * 1.0e9 # band width in nm
            f  = atmos.band_u_lw[1, ba] + atmos.band_u_sw[1, ba] # raw flux in W m-2
            ff = f / w * 1000.0 # converted to erg s-1 cm-2 nm-1
            ye[ba] = ff
        end

        # Get planck function values 
        plot_planck = false
        nsamps = 300
        if planck_tmp > 1.0 
            plot_planck = true
            xp = 10 .^ range( log10(xe[1]), stop=log10(xe[end]), length=nsamps)
            yp = zeros(Float64, nsamps)
            for i in 1:nsamps 
                lambda = xp[i] * 1.0e-9 # metres

                # Calculate planck function value [W m-2 sr-1 m-1]
                # http://spiff.rit.edu/classes/phys317/lectures/planck.html
                yp[i] = 2.0 * phys.h_pl * phys.c_vac^2 / lambda^5.0   *   1.0 / ( exp(phys.h_pl * phys.c_vac / (lambda * phys.k_B * planck_tmp)) - 1.0)

                # Integrate solid angle (hemisphere), scale by albedo, convert units
                yp[i] = yp[i] * pi * 1.0e-9 # [W m-2 nm-1]
                yp[i] = yp[i] * (1.0 - atmos.albedo_s)
                yp[i] = yp[i] * 1000.0 # [erg s-1 cm-2 nm-1]
            end 
        end 

        # Make plot
        plt = plot(framestyle=:box, dpi=dpi, guidefontsize=9)

        if plot_planck
            plot!(plt, xp, yp, label="Surface",  color="brown3") # surface planck function
        end
        plot!(plt, xe, ye, label="Outgoing spectrum", color="black")  # emission spectrum

        xlims  = (minimum(xe), min(maximum(xe), 50000.0))
        xticks = 10.0 .^ round.(Int,range( log10(xlims[1]), stop=log10(xlims[2]), step=1))

        ylims  = (minimum(ye) / 2, maximum(ye) * 2)
        yticks = 10.0 .^ round.(Int,range( log10(ylims[1]), stop=log10(ylims[2]), step=1))

        xlabel!(plt, "Wavelength [nm]")
        ylabel!(plt, "Spectral flux density [erg s-1 cm-2 nm-1]")
        yaxis!(plt, yscale=:log10, ylims=ylims, yticks=yticks)
        xaxis!(plt, xscale=:log10, xlims=xlims, xticks=xticks, minorgrid=true)

        savefig(plt, fname)

        return nothing 
    end 

    """
    Plot contribution function (per band)
    """
    function plot_contfunc(atmos, fname; dpi::Int=250)

        # Check that we have data 
        if !atmos.is_out_lw
            error("Cannot plot contribution function because radiances have not been calculated")
        end

        # Get data
        x = zeros(Float64, atmos.nbands)    # band centres (reverse order)
        y = zeros(Float64, atmos.nlev_c)    # pressure levels
        z = zeros(Float64, (atmos.nlev_c,atmos.nbands))

        # x value - band centres [nm]
        for ba in 1:atmos.nbands
            br = atmos.nbands - ba + 1
            x[br] = 0.5 * (atmos.bands_min[ba] + atmos.bands_max[ba]) * 1.0e9
        end

        # y value - pressures [bar]
        for i in 1:atmos.nlev_c 
            y[i] = atmos.p[i] * 1.0e-5
        end 

        # z value - log'd and normalised contribution function 
        cf_min = 1.0e-9
        for ba in 1:atmos.nbands
            for i in 1:atmos.nlev_c 
                br = atmos.nbands - ba + 1
                z[i,br] = log10(max(atmos.contfunc_norm[i,ba],cf_min))
            end 
        end 

        # Make plot
        plt = plot(framestyle=:box, dpi=dpi, guidefontsize=9, colorbar_title="log " * L"\widehat {cf}(\lambda, p)")

        heatmap!(plt, x,y,z, c=:devon)

        xlims  = (minimum(x), maximum(x))
        xticks = 10.0 .^ round.(Int,range( log10(xlims[1]), stop=log10(xlims[2]), step=1))
        xlabel!(plt, "Wavelength [nm]")
        xaxis!(plt, xscale=:log10, xlims=xlims, xticks=xticks, minorgrid=true)

        ylims  = (y[1], y[end])
        yticks = 10.0 .^ round.(Int,range( log10(ylims[1]), stop=log10(ylims[2]), step=1))
        ylabel!(plt, "Pressure [bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10, yticks=yticks, ylims=ylims, minorgrid=true)

        savefig(plt, fname)

        return nothing 
    end 

    """
    Plot spectral albedo (ratio of LW_UP to SW_DN)
    """
    function plot_albedo(atmos, fname; dpi::Int=250)

        # Check that we have data 
        if !(atmos.is_out_lw && atmos.is_out_sw)
            error("Cannot plot contribution function because radiances have not been calculated")
        end

        # Get data
        x = zeros(Float64, atmos.nbands)
        y = zeros(Float64, atmos.nbands)

        for ba in 1:atmos.nbands
            # x value - band centres [nm]
            x[ba] = 0.5 * (atmos.bands_min[ba] + atmos.bands_max[ba]) * 1.0e9
            
            # y value - spectral albedo [dimensionless]
            y[ba] = atmos.band_u_lw[1, ba]/atmos.band_d_sw[1, ba]
        end

        # Make plot
        plt = plot(framestyle=:box, dpi=dpi, guidefontsize=9)

        plot!(plt, x, y, label="", color="black")

        xlabel!(plt, "Wavelength [nm]")
        ylabel!(plt, "Spectral albedo")

        savefig(plt, fname)

        return nothing 
    end 

    """
    Animate output frames from solver, using ffmpeg
    """
    function anim_solver(atmos)

        # Command line format:
        # bash> ffmpeg -framerate 16 -i out/zzframe_%04d.png -pix_fmt yuv420p -y out/anim.mp4

        runtime = 15.0 # seconds

        # Find output files
        out = atmos.OUT_DIR
        frames = glob("zzframe_*.png",out)
        nframes = length(frames)

        # Animating fluxes?
        frames2 = glob("zyframe_*.png",out)
        fluxes = length(frames2) > 0

        # Create animation
        if nframes < 1
            println("WARNING: Cannot animate solver because no output frames were found")
        else 
            fps = max(nframes/runtime, 5)
            run(`ffmpeg -loglevel quiet -framerate $fps -pattern_type glob -i "$out/zzframe_*.png" -pix_fmt yuv420p -y $out/anim_tmp.mp4`)
            if fluxes
                run(`ffmpeg -loglevel quiet -framerate $fps -pattern_type glob -i "$out/zyframe_*.png" -pix_fmt yuv420p -y $out/anim_flx.mp4`)
            end
        end

        return nothing
    end


end

