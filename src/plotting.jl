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
    using Printf

    import atmosphere
    import phys

    # Symmetric log
    function _symlog(v::Float64)::Float64
        if abs(v) < 1.0
            return 0.0
        end 
        return sign(v)*max(log10(abs(v)), 0.0)
    end 

    # Int to string 
    function _intstr(v::Int)::String
        return @sprintf("%d",v)
    end 

    """
    Plot the temperature-pressure profile.
    """
    function plot_pt(atmos::atmosphere.Atmos_t, fname::String; dpi::Int=250, incl_magma::Bool=false, condensates::Array=[])
        
        # Interleave cell-centre and cell-edge arrays
        arr_P, arr_T = atmosphere.get_interleaved_pt(atmos)
        arr_P *= 1e-5

        ylims  = (arr_P[1]/1.5, arr_P[end]*1.5)
        yticks = 10.0 .^ round.(Int,range( log10(ylims[1]), stop=log10(ylims[2]), step=1))

        # Create plot
        plt = plot(framestyle=:box, ylims=ylims, yticks=yticks, legend=:outertopright, dpi=dpi, size=(500,400), guidefontsize=9)

        # Plot condensation curves 
        if length(condensates) > 0
            sat_t = zeros(Float64, atmos.nlev_l)
            for c in condensates
                for i in 1:atmos.nlev_l
                    sat_t[i] = phys.calc_Tdew(c, atmos.pl[i])
                end 
                plot!(plt, sat_t, atmos.pl*1e-5, lc=phys.lookup_color[c], ls=:dot, label=phys.lookup_pretty[c])
            end 
        end

        # Plot tmp_magma 
        if incl_magma
            scatter!(plt, [atmos.tmp_magma], [atmos.pl[end]*1e-5], color="cornflowerblue", label=L"T_m") 
        end

        # Plot tmp_surf 
        scatter!(plt, [atmos.tmp_surf], [atmos.pl[end]*1e-5], color="brown3", label=L"T_s")

        # Plot profile 
        plot!(plt, arr_T, arr_P, lc="black", lw=2, label=L"T(p)")

        # Decorate
        xlabel!(plt, "Temperature [K]")
        ylabel!(plt, "Pressure [bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)

        if !isempty(fname)
            savefig(plt, fname)
        end
        return plt 
    end

    """
    Plot the VMRs of the atmosphere at each cell-centre location.
    """
    function plot_x(atmos::atmosphere.Atmos_t, fname; dpi::Int=250)
        
        arr_P = atmos.p .* 1.0e-5 # Convert Pa to bar
        ylims  = (arr_P[1]/1.5, arr_P[end]*1.5)
        yticks = 10.0 .^ round.(Int,range( log10(ylims[1]), stop=log10(ylims[2]), step=1))

        
        # Create plot
        plt = plot(framestyle=:box, ylims=ylims, yticks=yticks, dpi=dpi, legend=:outertopright, size=(500,400), guidefontsize=9)

        # Plot mole fractions for each gas
        xmin::Float64 = 0.1
        this_min::Float64 = 0.0
        for (i_gas,gas) in enumerate(atmos.gases)  
            if haskey(phys.lookup_color, gas)
                c = phys.lookup_color[gas]
            else 
                c = "#"*bytes2hex(rand(UInt8, 3))
            end 
            plot!(atmos.layer_x[1:end,i_gas], arr_P, label=phys.lookup_pretty[gas], lw=2.5, linealpha=0.7, color=c)
            this_min = minimum(atmos.layer_x[1:end,i_gas])
            if this_min > 1.0e-90
                xmin = min(xmin, this_min)
            end
        end

        xlims  = (max(xmin, 1.0e-20)*0.5, 1.2)
        xticks = 10.0 .^ round.(Int,range( log10(xlims[1]), stop=0, step=1))

        # Set figure properties
        xlabel!(plt, "Mole fraction")
        xaxis!(plt, xscale=:log10, xlims=xlims, xticks=xticks)

        ylabel!(plt, "Pressure [bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)
                
        if !isempty(fname)
            savefig(plt, fname)
        end
        return plt 
    end

    """
    Plot the current temperature-pressure profile, current heating rates, and
    optionally the previous states that the atmosphere has taken.
    """
    function plot_solver(atmos::atmosphere.Atmos_t, fname::String; hist_tmpl::Array=[], incl_magma::Bool=false, dpi::Int=250, step::Int=-1)

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
        scatter!(plt1, [atmos.tmp_surf],     [atmos.pl[end]*1e-5], color="brown3",         label=L"T_s") 

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
    function plot_fluxes(atmos::atmosphere.Atmos_t, fname::String; dpi::Int=250, 
                            incl_eff::Bool=false, incl_mlt::Bool=true, incl_cdct::Bool=true, incl_phase::Bool=true
                        )

        arr_P = atmos.pl .* 1.0e-5 # Convert Pa to bar
        ylims  = (arr_P[1]*0.95, arr_P[end]*2.0)
        yticks = 10.0 .^ round.(Int,range( log10(ylims[1]), stop=log10(ylims[2]), step=1))

        max_fl = log10(max(100.0, maximum(abs.(atmos.flux_tot)), maximum(atmos.flux_u), maximum(atmos.flux_d)))
        xticks_pos = unique(ceil.(Int, range( start=1, stop=max_fl+1, step=1)))
        xticks = unique(vcat(-1.0.*reverse(xticks_pos), 0.0, xticks_pos))
        xlims = (-xticks_pos[end], xticks_pos[end])
        xticklabels = _intstr.(round.(Int, abs.(xticks)))

        w = 2.0
        plt = plot(legend=:outertopright, framestyle=:box, ylims=ylims, yticks=yticks, xticks=(xticks, xticklabels), xlims=xlims, dpi=dpi, guidefontsize=9)

        col_r = "#c0c0c0"
        col_n = "#000000"
        col_c = "#6495ed"
        col_t = "#ff4400"
        col_o = "#66CD00"
        col_p = "#ecb000"

        alpha = 0.7

        # Legend dummy plots
        plot!(plt, [-9e99, -8e99], [-9e99, -8e99], ls=:dot,   lw=w, lc=col_r, label="SW")
        plot!(plt, [-9e99, -8e99], [-9e99, -8e99], ls=:dash,  lw=w, lc=col_r, label="LW")
        plot!(plt, [-9e99, -8e99], [-9e99, -8e99], ls=:solid, lw=w, lc=col_r, label="LW+SW")
        plot!(plt, [-9e99, -8e99], [-9e99, -8e99], ls=:solid, lw=w, lc=col_n, label="UP+DN")

        # Zero line 
        plot!(plt, [0.0, 0.0], [arr_P[1], arr_P[end]], lw=0.4, lc="black", label="")

        # Effective flux
        if incl_eff
            plot!(plt, [_symlog(atmos.flux_eff)], [arr_P[1], arr_P[end]], ls=:dashdot, lw=0.4, lc="black", label="EFF")
        end

        # LW component
        if atmos.is_out_lw
            plot!(plt, _symlog.(-1.0*atmos.flux_d_lw), arr_P, label="", lw=w, lc=col_r, ls=:dash, linealpha=alpha)
            plot!(plt, _symlog.(     atmos.flux_u_lw), arr_P, label="", lw=w, lc=col_r, ls=:dash, linealpha=alpha)
        end

        # SW component
        if atmos.is_out_sw
            plot!(plt, _symlog.(-1.0.*atmos.flux_d_sw), arr_P, label="", lw=w, lc=col_r, ls=:dot, linealpha=alpha)
            plot!(plt, _symlog.(      atmos.flux_u_sw), arr_P, label="", lw=w, lc=col_r, ls=:dot, linealpha=alpha)
        end 

        # Net radiative fluxes
        if atmos.is_out_lw && atmos.is_out_sw
            plot!(plt, _symlog.(      atmos.flux_u), arr_P, label="", lw=w, lc=col_r, ls=:solid, linealpha=alpha)
            plot!(plt, _symlog.(-1.0.*atmos.flux_d), arr_P, label="", lw=w, lc=col_r, ls=:solid, linealpha=alpha)
            plot!(plt, _symlog.(      atmos.flux_n), arr_P, label="", lw=w, lc=col_n, ls=:solid, linealpha=alpha)
        end 

        # Convective flux (MLT)
        if incl_mlt
            plot!(plt, _symlog.(atmos.flux_cdry), arr_P, label="C_DRY", lw=w*1.2, lc=col_c, ls=:solid, linealpha=alpha)
        end 

        # Conduction 
        if incl_cdct
            plot!(plt, _symlog.(atmos.flux_cdct), arr_P, label="CNDCT", lw=w*1.2, lc=col_o, ls=:solid, linealpha=alpha)
        end 

        # Condensation 
        if incl_phase
            plot!(plt, _symlog.(atmos.flux_p), arr_P, label="PHASE", lw=w*1.2, lc=col_p, ls=:solid, linealpha=alpha)
        end

        # Sensible heat
        if atmos.flux_sens != 0.0
            scatter!(plt, [_symlog(atmos.flux_sens)], [arr_P[end]], markershape=:utriangle, markercolor=col_r, label="SENS")
        end

        # Total flux
        plot!(plt, _symlog.(atmos.flux_tot), arr_P, label="TOTAL", lw=w, lc=col_t, ls=:solid, linealpha=alpha)

        # Overplot convection and condensation mask
        for i in 1:atmos.nlev_c
            if atmos.mask_c[i] > 0
                scatter!(plt, [0.0], [arr_P[i]], opacity=0.9, markersize=2, color=col_c, label="")
            end 
            if atmos.mask_p[i] > 0
                scatter!(plt, [0.0], [arr_P[i]], opacity=0.9, markersize=2, color=col_p, label="")
            end 
        end

        # Labels 
        annotate!(plt, xlims[1]/2.0, arr_P[end], text("Downward", :black, :center, 9))
        annotate!(plt, xlims[2]/2.0, arr_P[end], text("Upward"  , :black, :center, 9))

        # Finalise + save
        xlabel!(plt, "log Unsigned flux [W m-2]")
        ylabel!(plt, "Pressure [bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)

        if !isempty(fname)
            savefig(plt, fname)
        end
        return plt 
    end

    """
    Plot emission spectrum at the TOA
    """
    function plot_emission(atmos::atmosphere.Atmos_t, fname::String; dpi::Int=250, incl_surf::Bool=true)

        # Check that we have data 
        if !(atmos.is_out_lw && atmos.is_out_sw)
            error("Cannot plot emission spectrum because radiances have not been calculated")
        end

        # Get emission spectrum data
        xe = zeros(Float64, atmos.nbands)
        yt = zeros(Float64, atmos.nbands)
        yl = zeros(Float64, atmos.nbands)
        ys = zeros(Float64, atmos.nbands)
        for ba in 1:atmos.nbands
            # x value - band centres [nm]
            xe[ba] = 0.5 * (atmos.bands_min[ba] + atmos.bands_max[ba]) * 1.0e9
            
            # y value - spectral flux [erg s-1 cm-2 nm-1]
            w  = (atmos.bands_max[ba] - atmos.bands_min[ba]) * 1.0e9 # band width in nm
            yl[ba] = atmos.band_u_lw[1, ba] / w * 1000.0 # converted to erg s-1 cm-2 nm-1
            ys[ba] = atmos.band_u_sw[1, ba] / w * 1000.0 # converted to erg s-1 cm-2 nm-1
            yt[ba] = yl[ba] + ys[ba]

        end

        # Get planck function values 
        if incl_surf
            planck_tmp::Float64 = atmos.tmp_surf
            nsamps::Int = 300
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

        if incl_surf
            plot!(plt, xp, yp, label="Surface",  color="green") # surface planck function
        end

        plot!(plt, xe, ys, label="SW spectrum", color="blue")  
        plot!(plt, xe, yl, label="LW spectrum", color="red" ) 
        plot!(plt, xe, yt, label="Total spectrum", color="black")  # emission spectrum 

        xlims  = ( max(1.0e-10,minimum(xe)), min(maximum(xe), 70000.0))
        xticks = 10.0 .^ round.(Int,range( log10(xlims[1]), stop=log10(xlims[2]), step=1))

        ylims  = (max(1.0e-10,minimum(yt)) / 2, max(maximum(yt),maximum(yp)) * 2)
        yticks = 10.0 .^ round.(Int,range( log10(ylims[1]), stop=log10(ylims[2]), step=1))

        xlabel!(plt, "Wavelength [nm]")
        ylabel!(plt, "Spectral flux density [erg s⁻¹ cm⁻² nm⁻¹]")
        yaxis!(plt, yscale=:log10, ylims=ylims, yticks=yticks)
        xaxis!(plt, xscale=:log10, xlims=xlims, xticks=xticks, minorgrid=true)

        if !isempty(fname)
            savefig(plt, fname)
        end
        return plt 
    end 

    """
    Plot contribution function (per band)
    """
    function plot_contfunc(atmos::atmosphere.Atmos_t, fname::String; dpi::Int=250)

        # Check that we have data 
        if !atmos.is_out_lw
            error("Cannot plot contribution function because radiances have not been calculated")
        end

        # Get data
        x = zeros(Float64, atmos.nbands)    # band centres (reverse order)
        y = zeros(Float64, atmos.nlev_c)    # pressure levels
        z = zeros(Float64, (atmos.nlev_c,atmos.nbands))

        # Reversed?
        reversed::Bool = (atmos.bands_min[1] > atmos.bands_min[end])

        # x value - band centres [nm]
        for ba in 1:atmos.nbands
            if reversed
                br = atmos.nbands - ba + 1
            else 
                br = ba 
            end
            x[br] = 0.5 * (atmos.bands_min[ba] + atmos.bands_max[ba]) * 1.0e9
        end

        # y value - pressures [bar]
        for i in 1:atmos.nlev_c 
            y[i] = atmos.p[i] * 1.0e-5
        end 

        # z value - log'd and normalised contribution function 
        cf_min = 1.0e-9
        for ba in 1:atmos.nbands
            if reversed
                br = atmos.nbands - ba + 1
            else 
                br = ba 
            end
            for i in 1:atmos.nlev_c 
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

        if !isempty(fname)
            savefig(plt, fname)
        end
        return plt 
    end 

    """
    Plot spectral albedo (ratio of SW_UP to SW_DN)
    """
    function plot_albedo(atmos::atmosphere.Atmos_t, fname::String; dpi::Int=250)

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
            y[ba] = 100.0 * atmos.band_u_sw[1, ba]/atmos.band_d_sw[1, ba]
        end

        # Make plot
        plt = plot(framestyle=:box, dpi=dpi, guidefontsize=9)

        plot!(plt, x, y, label="", color="black")

        xlims  = (200.0, 1000.0)
        xticks = range( xlims[1], xlims[2], step=100.0)
        xaxis!(plt, xlims=xlims, xticks=xticks, minorgrid=true)

        xlabel!(plt, "Wavelength [nm]")
        ylabel!(plt, "Spectral albedo [%]")

        if !isempty(fname)
            savefig(plt, fname)
        end
        return plt 
    end 

    """
    Animate output frames from solver, using ffmpeg
    """
    function anim_solver(atmos::atmosphere.Atmos_t)

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
            @warn "Cannot animate solver because no output frames were found"
        else 
            fps = max(nframes/runtime, 5)
            run(`ffmpeg -loglevel quiet -framerate $fps -pattern_type glob -i "$out/zzframe_*.png" -pix_fmt yuv420p -y $out/anim_tmp.mp4`)
            if fluxes
                run(`ffmpeg -loglevel quiet -framerate $fps -pattern_type glob -i "$out/zyframe_*.png" -pix_fmt yuv420p -y $out/anim_flx.mp4`)
            end
        end

        return nothing
    end


end # end module plotting

