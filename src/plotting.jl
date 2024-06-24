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
    using Printf
    using FFMPEG

    import ..atmosphere
    import ..phys

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
    function plot_pt(atmos::atmosphere.Atmos_t, fname::String; 
                            dpi::Int=250,
                            size_x::Int=500, size_y::Int=400,
                            incl_magma::Bool=false, 
                            title::String="")
        
        # Interleave cell-centre and cell-edge arrays
        ylims  = (1e-5*atmos.pl[1]/1.5, 1e-5*atmos.pl[end]*1.5)
        yticks = 10.0 .^ round.(Int,range( log10(ylims[1]), stop=log10(ylims[2]), step=1))

        # Create plot
        plt = plot(framestyle=:box, ylims=ylims, yticks=yticks, legend=:outertopright, dpi=dpi, size=(size_x,size_y), guidefontsize=9, titlefontsize=9)

        # Plot phase boundary 
        if length(atmos.condensates) > 0
            sat_t::Array{Float64,1} = zeros(Float64, atmos.nlev_c)
            crt_i::Int = atmos.nlev_c  # index at which criticality occurs
            for c in atmos.condensates

                # for each level...
                i = 1
                crt_i = atmos.nlev_c
                while (i < crt_i) && (i <= atmos.nlev_c)
                    # get dew point temperature from this partial pressure 
                    sat_t[i] = phys.get_Tdew(atmos.gas_dat[c], atmos.pl[i]) 

                    # check if supercritical
                    if sat_t[i] >= atmos.gas_dat[c].T_crit
                        crt_i = i
                    end 

                    i += 1
                end 

                # plot phase boundary for this condensate, stopping at critical point
                plot!(plt, sat_t[1:crt_i], atmos.p[1:crt_i]*1e-5, lc=atmos.gas_dat[c].plot_color, ls=:dot, label=atmos.gas_dat[c].plot_label)

                # plot critical temperature 
                # scatter!(plt, [atmos.gas_dat[c].T_crit], [atmos.p[crt_i]*1e-5], mc=atmos.gas_dat[c].plot_color, label="")
            end 
        end

        # Plot tmp_magma 
        if incl_magma
            scatter!(plt, [atmos.tmp_magma], [atmos.pl[end]*1e-5], color="cornflowerblue", label=L"T_m") 
        end

        # Plot tmp_surf 
        scatter!(plt, [atmos.tmp_surf], [atmos.pl[end]*1e-5], color="brown3", label=L"T_s")

        # Plot profile 
        plot!(plt, atmos.tmpl, atmos.pl*1e-5, lc="black", lw=2, label=L"T(p)")

        # Decorate
        xlabel!(plt, "Temperature [K]")
        ylabel!(plt, "Pressure [bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)
        if !isempty(title)
            title!(plt, title)
        end 

        if !isempty(fname)
            savefig(plt, fname)
        end
        return plt 
    end

    """
    Plot the VMRs of the atmosphere at each cell-centre location.
    """
    function plot_vmr(atmos::atmosphere.Atmos_t, fname::String; 
                            dpi::Int=250,
                            size_x::Int=500, size_y::Int=400)
        
        arr_P = atmos.p .* 1.0e-5 # Convert Pa to bar
        ylims  = (arr_P[1]/1.5, arr_P[end]*1.5)
        yticks = 10.0 .^ round.(Int,range( log10(ylims[1]), stop=log10(ylims[2]), step=1))
        
        # Create plot
        plt = plot(framestyle=:box, ylims=ylims, yticks=yticks, dpi=dpi, legend=:outertopright, size=(size_x,size_y), guidefontsize=9)

        # Plot log10 mole fractions for each gas
        xmin::Float64 = -2
        this_min::Float64 = 0.0

        for gas in atmos.gas_names
            # get VMR
            x_arr = log10.(clamp.(atmos.gas_vmr[gas][:],1e-100, 1e100))
            this_min = minimum(x_arr)
            if this_min > -90
                xmin = min(xmin, this_min)
            end

            # plot gas
            plot!(x_arr, arr_P, label=phys.pretty_name(gas), lw=2.5, linealpha=0.7, color=phys.pretty_color(gas))
        end

        xlims  = (max(xmin-0.1, -12), 0.1)
        xticks = round.(Int,range( xlims[1], stop=0, step=1))

        # Set figure properties
        xlabel!(plt, "log(Mole fraction)")
        xaxis!(plt, xlims=xlims, xticks=xticks)

        ylabel!(plt, "Pressure [bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)
                
        if !isempty(fname)
            savefig(plt, fname)
        end
        return plt 
    end

    """
    Plot the fluxes at each pressure level
    """
    function plot_fluxes(atmos::atmosphere.Atmos_t, fname::String; dpi::Int=250, 
                            size_x::Int=550, size_y::Int=400,
                            incl_eff::Bool=false, incl_mlt::Bool=true, 
                            incl_cdct::Bool=true, incl_latent::Bool=true,
                            title::String=""
                        )

        arr_P = atmos.pl .* 1.0e-5 # Convert Pa to bar
        ylims  = (arr_P[1]*0.95, arr_P[end]*2.0)
        yticks = 10.0 .^ round.(Int,range( log10(ylims[1]), stop=log10(ylims[2]), step=1))

        max_fl = log10(max(100.0, maximum(abs.(atmos.flux_tot)), maximum(atmos.flux_u), maximum(atmos.flux_d)))
        xticks_pos = unique(ceil.(Int, range( start=1, stop=max_fl+1, step=1)))
        xticks = unique(vcat(-1.0.*reverse(xticks_pos), 0.0, xticks_pos))
        xlims = (-xticks_pos[end], xticks_pos[end])
        xticklabels = _intstr.(round.(Int, abs.(xticks)))

        plt = plot(legend=:outertopright, framestyle=:box, ylims=ylims, yticks=yticks, xticks=(xticks, xticklabels), xlims=xlims, dpi=dpi, guidefontsize=9, size=(size_x,size_y), titlefontsize=9)

        col_r::String = "#c0c0c0"
        col_n::String = "#000000"
        col_c::String = "#6495ed"
        col_t::String = "#ff4400"
        col_o::String = "#66CD00"
        col_p::String = "#ecb000"

        alpha = 0.7
        w = 2.0

        # Legend dummy plots
        plot!(plt, [-9e99, -8e99], [-9e99, -8e99], ls=:dot,   lw=w, lc=col_r, label="SW")
        plot!(plt, [-9e99, -8e99], [-9e99, -8e99], ls=:dash,  lw=w, lc=col_r, label="LW")
        plot!(plt, [-9e99, -8e99], [-9e99, -8e99], ls=:solid, lw=w, lc=col_r, label="LW+SW")
        plot!(plt, [-9e99, -8e99], [-9e99, -8e99], ls=:solid, lw=w, lc=col_n, label="UP-DN")

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
            plot!(plt, _symlog.(atmos.flux_cdry), arr_P, label="Convect", lw=w*1.2, lc=col_c, ls=:solid, linealpha=alpha)
        end 

        # Conduction 
        if incl_cdct
            plot!(plt, _symlog.(atmos.flux_cdct), arr_P, label="Conduct", lw=w*1.2, lc=col_o, ls=:solid, linealpha=alpha)
        end 

        # Condensation 
        if incl_latent
            plot!(plt, _symlog.(atmos.flux_l), arr_P, label="Latent", lw=w*1.2, lc=col_p, ls=:solid, linealpha=alpha)
        end

        # Sensible heat
        if atmos.flux_sens != 0.0
            scatter!(plt, [_symlog(atmos.flux_sens)], [arr_P[end]], markershape=:utriangle, markercolor=col_r, label="Sensible")
        end

        # Total flux
        plot!(plt, _symlog.(atmos.flux_tot), arr_P, label="Total", lw=w, lc=col_t, ls=:solid, linealpha=alpha)

        # Overplot convection and condensation mask
        #    by indicating it with scatter points of the corresponding colour
        for i in 1:atmos.nlev_c
            if atmos.mask_c[i] > 0
                scatter!(plt, [-0.2], [arr_P[i]], opacity=0.9, markersize=2, msw=0.5, color=col_c, label="")
            end 
            if atmos.mask_l[i] > 0
                scatter!(plt, [0.2], [arr_P[i]], opacity=0.9, markersize=2, msw=0.5, color=col_p, label="")
            end 
        end

        # Labels 
        annotate!(plt, xlims[1]/2.0, arr_P[end], text("Downward", :black, :center, 9))
        annotate!(plt, xlims[2]/2.0, arr_P[end], text("Upward"  , :black, :center, 9))

        # Finalise + save
        xlabel!(plt, "log Unsigned flux [W m⁻²]")
        ylabel!(plt, "Pressure [bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)
        if !isempty(title)
            title!(plt, title)
        end

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
        xe::Array{Float64, 1} = zeros(Float64, atmos.nbands)
        yt::Array{Float64, 1} = zeros(Float64, atmos.nbands)
        yl::Array{Float64, 1} = zeros(Float64, atmos.nbands)
        ys::Array{Float64, 1} = zeros(Float64, atmos.nbands)
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

        plot!(plt, xe, ys, lw=0.9, label="SW spectrum", color="blue")  
        plot!(plt, xe, yl, lw=0.9, label="LW spectrum", color="red" ) 
        plot!(plt, xe, yt, lw=0.5, label="Total spectrum", color="black")  # emission spectrum 

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
        x::Array{Float64, 1} = zeros(Float64, atmos.nbands)    # band centres (reverse order)
        y::Array{Float64, 1} = zeros(Float64, atmos.nlev_c)    # pressure levels
        z::Array{Float64, 2} = zeros(Float64, (atmos.nlev_c,atmos.nbands))

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
        x::Array{Float64, 1} = zeros(Float64, atmos.nbands)
        y::Array{Float64, 1} = zeros(Float64, atmos.nbands)

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
    Combined plot 
    """
    function combined(plt_pt, plt_fl, plt_mr, info::String, fname::String; dpi::Int=180, size_x::Int=800, size_y::Int=650)

        plt_info = plot(legend=false, showaxis=false, grid=false)
        annotate!(plt_info, (0.02, 0.7, text(info, family="Courier", :black, :left, 10)))

        plt = plot(plt_pt, plt_fl, plt_mr, plt_info, dpi=dpi, layout=(2,2), size=(size_x, size_y))

        if !isempty(fname)
            savefig(plt, fname)
        end
        return plt 
    end 

    """
    Animate output frames from solver, using ffmpeg
    """
    function animate(atmos::atmosphere.Atmos_t, duration::Float64=12.0)

        # Find output files
        out::String = atmos.OUT_DIR
        frames = glob("frames/*.png",out)
        nframes::Int = length(frames)

        # Create animation
        if nframes < 1
            @warn "Cannot animate solver because no output frames were found"
        else 
            fps::Float64 = nframes/duration*1.0
            @ffmpeg_env run(`$(FFMPEG.ffmpeg) -loglevel quiet -framerate $fps -pattern_type glob -i "$out/frames/*.png" -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" -y $out/animation.mp4`)
        end

        return nothing
    end

    """
    Plot jacobian matrix 
    """
    function jacobian(b::Array{Float64,2}, fname::String; 
                            perturb::Array{Bool,1}=Bool[],
                            dpi::Int=200, size_x::Int=600, size_y::Int=500)
        
        lim::Float64 = maximum(abs.(b))     # colourbar limits
        l::Int = length(perturb)            # show perturbed levels?

        plt = plot(framestyle=:box, dpi=dpi, 
                    guidefontsize=9, size=(size_x, size_y),
                    title="∂r/∂x [W m⁻² K⁻¹]",
                    clim=(-lim,lim), grid=false, yflip=true)

        # show jacobian 
        heatmap!(plt, b, color=:RdBu)

        # show perturbed levels 
        if l > 0
            scatter!(plt, collect(1:l)[perturb], ones(Float64, l)[perturb]*(l+1),
                        color=:black,label="",markershape=:utriangle)
        end 

        xlabel!(plt, "Cell-centre index")
        ylabel!(plt, "Cell-centre index")

        if !isempty(fname)
            savefig(plt, fname)
        end

        return plt 
    end 

end # end module plotting

