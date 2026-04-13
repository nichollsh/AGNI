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
    using Printf
    using FFMPEG
    using Statistics
    import Glob:glob

    import ..atmosphere
    import ..phys
    import ..multicol

    # Allowed plot file extensions
    const ALLOWED_EXTS::Set{String} = Set(["png", "pdf", "svg"])

    # Colors
    const col_r::String = "#c0c0c0"
    const col_n::String = "#000000"
    const col_c::String = "#6495ed"
    const col_t::String = "#ff4400"
    const col_o::String = "#66CD00"
    const col_p::String = "#ecb000"
    const col_d::String = "#8B008B"

    # Default plotting configuration
    const la::Float64 = 0.7
    const lw::Float64 = 1.9
    const size_x_default::Int64 = 500
    const size_y_default::Int64 = 400
    const plt_default = Dict(:fontfamily => "sans-serif",
                             :framestyle => :box,
                             :grid       => true,
                             :guidefontsize => 9,
                             :titlefontsize => 9,
                             :dpi => 240)

    """
    **Apply a signed symmetric log10 transform, returning zero for |v| < thresh.**

    Arguments:
    - `v::Float64`                  value to transform
    - `thresh::Float64`             threshold below which to return zero (default: 1.0)

    Returns:
    - `Float64`                     log10 transformed value
    """
    function _symlog(v::Float64; thresh::Float64=1.0)::Float64
        if abs(v) < thresh
            return 0.0
        end
        return sign(v)*max(log10(abs(v)), 0.0)
    end

    """
    Format an integer as a plain decimal string.
    """
    function _intstr(v::Int64)::String
        return @sprintf("%d",v)
    end

    """
    Return `(p_top, p_bot)` pressure axis limits [bar] for a standard profile plot.
    """
    function _get_ylims(atmos::atmosphere.Atmos_t)::Tuple
        return (1e-5*atmos.pl[1]/1.5, 1e-5*max(atmos.p_oboa, atmos.p_boa)*1.5 )
    end

    """
    Return pressure tick mark values [bar] on a log10 scale for a standard profile plot.
    """
    function _get_yticks(atmos::atmosphere.Atmos_t)::Array
        ylims = _get_ylims(atmos)
        return 10.0 .^ round.(Int,range( log10(ylims[1]), stop=log10(ylims[2]), step=1))
    end

    # Reusable plotting snippets
    macro _plt_pboa()
        return esc(:(hline!(plt, [atmos.p_boa/1e5 ], label="", color="black",  ls=:solid)))
    end
    macro _plt_poboa()
        return esc(:(hline!(plt, [atmos.p_oboa/1e5], label="", color="black", ls=:dot)))
    end

    """
    **Plot the temperature-pressure and Kzz profile.**

    Arguments:
    - `atmos::atmosphere.Atmos_t`   atmosphere object
    - `fname::String`               filename to save the plot (if empty, does not save)
    - `size_x::Int64`               width of the plot in pixels
    - `size_y::Int64`               height of the plot in pixels
    - `incl_magma::Bool`            include the magma temperature as a scatter point
    - `title::String`               title for the plot
    """
    function plot_pt(atmos::atmosphere.Atmos_t, fname::String;
                            size_x::Int64=size_x_default, size_y::Int64=size_y_default,
                            incl_magma::Bool=false,
                            title::String="")

        y = atmos.pl ./ 1e5 # pressure -> bar

        # Create plot
        plt = plot(ylims=_get_ylims(atmos), yticks=_get_yticks(atmos),
                        legend=:outerbottomright,
                        size=(size_x,size_y); plt_default...)

        # Plot phase boundary and demixing binodal
        if atmos.condense_any
            sat_t::Array{Float64,1} = zeros(Float64, atmos.nlev_c)
            dmx_t::Array{Float64,1} = zeros(Float64, atmos.nlev_c)
            for c in atmos.condensates
                # demixing binodal for this condensate
                if atmos.demixing
                    for i in 1:atmos.nlev_c
                        dmx_t[i] = phys.get_Tdemix(atmos.gas_dat[c], atmos.p[i], atmos.gas_vmr[c][i])
                        if dmx_t[i] < 0.0
                            dmx_t[i] = NaN
                        end
                    end
                    if !all(isnan.(dmx_t))
                        plot!(plt, dmx_t, atmos.p*1e-5, lc=atmos.gas_dat[c].plot_color, ls=:dash,
                                label=atmos.gas_dat[c].plot_label*" dmx")
                    end
                end

                # saturation curve
                if atmos.gas_dat[c].no_sat
                    continue
                end
                for i in 1:atmos.nlev_c
                    sat_t[i] = phys.get_Tdew(atmos.gas_dat[c], atmos.p[i]*atmos.gas_vmr[c][i])
                    if sat_t[i] > atmos.gas_dat[c].T_crit-0.1
                        sat_t[i] = NaN
                    end
                end
                plot!(plt, sat_t, atmos.p*1e-5, lc=atmos.gas_dat[c].plot_color, ls=:dot,
                            label=atmos.gas_dat[c].plot_label*" sat")
            end
        end

        # Plot tmp_magma
        if incl_magma
            scatter!(plt, [atmos.tmp_magma], [atmos.p_boa/1e5],
                        color="cornflowerblue", label=L"T_m")
        end

        # Plot tmp_surf
        scatter!(plt, [atmos.tmp_surf], [atmos.p_boa/1e5], color="brown3", label=L"T_s")

        # Plot profile
        plot!(plt, atmos.tmpl, y, lc="black", lw=lw, label=L"T(p)")

        # Plot current surface pressure and original
        @_plt_pboa
        @_plt_poboa

        # Decorate
        xlims!(plt, (0.0, maximum(atmos.tmpl)+15.0))
        xlabel!(plt, "Temperature [K]")
        ylabel!(plt, "Pressure [bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)
        if !isempty(title)
            title!(plt, title)
        end

        # Add secondary x-axis for Kzz profile
        plt2 = twiny(plt)

        x = atmos.Kzz .* 1e4 # convert from m²/s to cm²/s for plotting
        xmin = 5.0
        if any(x .> 0.0)
            xmin = min(xmin, log10(minimum(x[x.>0])))
        end
        x = log10.(clamp.(x, 10.0 ^ xmin, Inf64))
        x_con = copy(x)
        x_rad = copy(x)
        for i in 2:atmos.nlev_l-1
            # convective regions, and adjacent layers
            x_con[i] = any(atmos.flux_cdry[i-1:i+1] .> 0.0) ? x[i] : NaN
            # radiative regions, strictly
            x_rad[i] = atmos.flux_cdry[i] <= 0.0 ? x[i] : NaN
        end

        plot!(plt2, [xmin, xmin], [1.0, 1.0], lc="darkgreen", lw=lw, ls=:solid, label=L"K_{zz}")
        plot!(plt2, x_con, y, lc="darkgreen", label="Con.", ls=:solid)
        plot!(plt2, x_rad, y, lc="darkgreen", label="Rad.", ls=:dot)

        xlabel!(plt2, "log₁₀ Kzz [cm²/s]")
        ylims!(plt2, _get_ylims(atmos))
        yticks!(plt2, _get_yticks(atmos))
        yflip!(plt2)
        yaxis!(plt2, yscale=:log10)
        xaxis!(plt2, xlims=(xmin, maximum(x)+1), grid=false, gridalpha=0.0)

        # ensures that axes spines are same size
        plot!(plt2, legend=:outertopright, tick_direction=:out)

        if !isempty(fname)
            savefig(plt, fname)
        end
        return plt
    end

    """
    **Plot the radius vs pressure profile.**

    Arguments:
    - `atmos::atmosphere.Atmos_t`   atmosphere object
    - `fname::String`               filename to save the plot (if empty, does not save)
    - `size_x::Int64`               width of the plot in pixels
    - `size_y::Int64`               height of the plot in pixels
    - `title::String`               title for the plot
    """
    function plot_radius(atmos::atmosphere.Atmos_t, fname::String;
                                size_x::Int64=size_x_default, size_y::Int64=size_y_default,
                                title::String="")

        # Create plot
        plt = plot(ylims=_get_ylims(atmos), yticks=_get_yticks(atmos),
                        legend=:outertopright,
                        size=(size_x,size_y); plt_default...)

        # Plot surface
        scatter!(plt, [atmos.rp*1e-3], [atmos.pl[end]*1e-5], color="brown3", label=L"P_s")

        # Plot cell-centres and cell-edges
        scatter!(plt, atmos.r*1e-3,  atmos.p*1e-5,  msa=0.0, msw=0, ms=1.2, shape=:diamond, label="Centres")
        scatter!(plt, atmos.rl*1e-3, atmos.pl*1e-5, msa=0.0, msw=0, ms=1.2, shape=:diamond, label="Edges")

        # Plot current surface pressure and original
        @_plt_pboa
        @_plt_poboa

        # Decorate
        xlabel!(plt, "Radius [km]")
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
    **Plot the cloud and aerosol mass mixing ratios.**

    Arguments:
    - `atmos::atmosphere.Atmos_t`   atmosphere object
    - `fname::String`               filename to save the plot (if empty, does not save)
    - `size_x::Int64`               width of the plot in pixels
    - `size_y::Int64`               height of the plot in pixels
    - `title::String`               title for the plot
    """
    function plot_cloud(atmos::atmosphere.Atmos_t, fname::String;
                            size_x::Int64=size_x_default, size_y::Int64=size_y_default,
                            title::String="")

        xlims = (-8.0, 0.0)
        xticks = collect(range(start=xlims[1], stop=xlims[2], step=1))

        y = atmos.p * 1e-5 # pressure -> bar

        # Create plot
        plt = plot( xlims=xlims, xticks=xticks,
                    ylims=_get_ylims(atmos), yticks=_get_yticks(atmos),
                    legend=:outertopright, size=(size_x,size_y); plt_default...)

        # Temperature profile for reference
        tmp_nrm = (atmos.tmp .- minimum(atmos.tmp))./(maximum(atmos.tmp)-minimum(atmos.tmp))
        @. tmp_nrm = xlims[1] + (xlims[2]-xlims[1])*tmp_nrm
        plot!(plt, tmp_nrm, y, lc="black",
                        linealpha=0.3, lw=lw, label=L"\hat{T}(p)")

        # Plot cloud profiles
        ls = atmos.control.l_cloud ? :solid : :dot
        plot!(plt, log10.(clamp.(atmos.cloud_arr_l,10^xlims[1],10^xlims[2])), y,
                    lw=lw, ls=ls, label="Cloud", linealpha=la)

        # Plot aerosol profiles
        ls = atmos.control.l_aerosol ? :solid : :dot
        for k_aer in keys(atmos.aerosol_arr_l)
            plot!(plt, log10.(clamp.(atmos.aerosol_arr_l[k_aer], 10^xlims[1], 10^xlims[2])), y,
                    lw=lw, ls=ls, label=k_aer, linealpha=la)
        end

        # Plot current surface pressure and original
        @_plt_pboa
        @_plt_poboa

        # Decorate
        xlabel!(plt, "log₁₀ Mass mixing ratio)")
        ylabel!(plt, "Pressure [bar]")
        yflip!(plt)
        yaxis!(plt, yscale=:log10)
        xaxis!(plt, xlims=xlims)
        if !isempty(title)
            title!(plt, title)
        end

        if !isempty(fname)
            savefig(plt, fname)
        end
        return plt
    end

    """
    **Plot the gas phase volume mixing ratios at each cell-centre location.**

    Arguments:
    - `atmos::atmosphere.Atmos_t`   atmosphere object
    - `fname::String`               filename to save the plot (if empty, does not save)
    - `size_x::Int64`               width of the plot in pixels
    - `size_y::Int64`               height of the plot in pixels
    """
    function plot_vmr(atmos::atmosphere.Atmos_t, fname::String;
                            size_x::Int64=size_x_default, size_y::Int64=size_y_default)

        # X-axis minimum allowed left-hand-side limit (log units)
        minmin_x::Float64 = -10

        arr_P = atmos.p .* 1.0e-5 # Convert Pa to bar

        # Create plot
        plt = plot(ylims=_get_ylims(atmos), yticks=_get_yticks(atmos),
                        legend=:outertopright,
                        size=(size_x,size_y); plt_default...)


        # Plot log10 VMRs for each gas
        gas_xsurf::Array = zeros(Float64, atmos.gas_num)
        gas::String = ""
        for i in 1:atmos.gas_num
            gas = atmos.gas_names[i]

            # store surface value
            gas_xsurf[i] = log10(clamp(atmos.gas_vmr[gas][end], eps(0.0), phys.BIGFLOAT))
        end

        num_plotted::Int64 = 0
        arr_x::Array{Float64, 1} = zeros(Float64, atmos.nlev_c)
        min_x::Float64 = -3
        for i in reverse(sortperm(gas_xsurf))
            gas = atmos.gas_names[i]
            col = atmos.gas_dat[gas].plot_color
            # Plot gases in order of descending abundance, so that the legend
            #    shows the most interesting gases at the top of the list.

            # Too many gases makes legend and plot unreadable
            if num_plotted > 20
                break
            end
            num_plotted += 1

            # Plot post-chemistry values, before rainout
            if atmos.condense_any
                @. arr_x = atmos.gas_cvmr[gas]
                if minimum(arr_x) < eps(0.0)
                    continue
                end
                @. arr_x = log10(arr_x)
                min_x = min(min_x, minimum(arr_x))
                plot!(arr_x, arr_P, label=nothing, linestyle=:dot,
                        lw=lw, linealpha=la, color=col)
            end


            # Plot runtime values, after rainout
            @. arr_x = atmos.gas_vmr[gas]
            if minimum(arr_x) < eps(0.0)
                continue
            end
            @. arr_x = log10(arr_x)
            min_x = min(min_x, minimum(arr_x))
            plot!(arr_x, arr_P,  label=atmos.gas_dat[gas].plot_label, linestyle=:solid,
                    lw=lw, linealpha=la, color=col)



            # Plot original value as surface scatter point
            scatter!([log10(atmos.gas_ovmr[gas][end])], [arr_P[end]],
                        opacity=0.9, markersize=2, msw=0.5,
                        color=atmos.gas_dat[gas].plot_color, label="")
        end

        # Plot current surface pressure and original
        @_plt_pboa
        @_plt_poboa

        xlims  = (max(min_x, minmin_x)-0.1, 0.1)
        xticks = round.(Int,range( xlims[1], stop=0, step=1))

        # Set figure properties
        xlabel!(plt, "log₁₀ Volume Mixing Ratio")
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
    **Plot the fluxes at each pressure level**

    Arguments:
    - `atmos::atmosphere.Atmos_t`  atmosphere object
    - `fname::String`              filename to save the plot (if empty, does not save)
    - `size_x::Int64`              width of the plot in pixels
    - `size_y::Int64`              height of the plot in pixels
    - `incl_eff::Bool`             whether to include the intrinsic (or interior) heat flux as a dashed line
    - `incl_mlt::Bool`             whether to include the convective flux as a solid line
    - `incl_cdct::Bool`            whether to include the conductive flux as a solid line
    - `incl_latent::Bool`          whether to include the latent heating flux as a solid line
    - `incl_deep::Bool`            whether to include the deep heating flux as a solid line
    - `title::String`              title for the plot
    """
    function plot_fluxes(atmos::atmosphere.Atmos_t, fname::String;
                            size_x::Int64=size_x_default, size_y::Int64=size_y_default,
                            incl_eff::Bool=false, incl_mlt::Bool=true,
                            incl_cdct::Bool=true, incl_latent::Bool=true,
                            incl_deep::Bool=true,
                            title::String=""
                        )

        arr_P = atmos.pl .* 1.0e-5 # Convert Pa to bar

        max_fl = log10(max(100.0, maximum(abs.(atmos.flux_tot)), maximum(atmos.flux_u), maximum(atmos.flux_d)))
        xticks_pos = unique(ceil.(Int, range( start=1, stop=max_fl+1, step=1)))
        xticks = unique(vcat(-1.0.*reverse(xticks_pos), 0.0, xticks_pos))
        xlims = (-xticks_pos[end], xticks_pos[end])
        xticklabels = _intstr.(round.(Int, abs.(xticks)))
        ylims = _get_ylims(atmos)

        plt = plot(legend=:outertopright,
                    ylims=ylims, yticks=_get_yticks(atmos),
                    xticks=(xticks, xticklabels), xlims=xlims,
                    size=(size_x,size_y); plt_default...)

        # Legend dummy plots
        plot!(plt, [-9e99, -8e99], [-9e99, -8e99], ls=:dot,   lw=lw, lc=col_r, label="SW")
        plot!(plt, [-9e99, -8e99], [-9e99, -8e99], ls=:dash,  lw=lw, lc=col_r, label="LW")
        plot!(plt, [-9e99, -8e99], [-9e99, -8e99], ls=:solid, lw=lw, lc=col_r, label="LW+SW")
        plot!(plt, [-9e99, -8e99], [-9e99, -8e99], ls=:solid, lw=lw, lc=col_n, label="UP-DN")

        # Zero line
        vline!(plt, [0.0], lw=0.4, lc="black", label="")

        # Indicate the target intrinsic (or interior) heat flux
        if incl_eff
            plot!(plt, [_symlog(atmos.flux_int)], [arr_P[1], arr_P[end]], ls=:dashdot, lw=0.4, lc="black", label="INT")
        end

        # LW component
        if atmos.is_out_lw
            plot!(plt, _symlog.(-1.0*atmos.flux_d_lw), arr_P, lw=lw, lc=col_r, ls=:dash, linealpha=la, label="")
            plot!(plt, _symlog.(     atmos.flux_u_lw), arr_P, lw=lw, lc=col_r, ls=:dash, linealpha=la, label="")
        end

        # SW component
        if atmos.is_out_sw
            plot!(plt, _symlog.(-1.0*atmos.flux_d_sw),  arr_P, lw=lw, lc=col_r, ls=:dot, linealpha=la, label="")
            plot!(plt, _symlog.(      atmos.flux_u_sw), arr_P, lw=lw, lc=col_r, ls=:dot, linealpha=la, label="")
        end

        # Net radiative fluxes
        if atmos.is_out_lw && atmos.is_out_sw
            plot!(plt, _symlog.(      atmos.flux_u), arr_P, lw=lw, lc=col_r, ls=:solid, linealpha=la, label="")
            plot!(plt, _symlog.(-1.0*atmos.flux_d),  arr_P, lw=lw, lc=col_r, ls=:solid, linealpha=la, label="")
            plot!(plt, _symlog.(      atmos.flux_n), arr_P, lw=lw, lc=col_n, ls=:solid, linealpha=la, label="")
        end

        # Convective flux (MLT)
        if incl_mlt
            plot!(plt, _symlog.(atmos.flux_cdry), arr_P, label="Convect", lw=lw*1.2, lc=col_c, ls=:solid, linealpha=la)
        end

        # Conduction
        if incl_cdct
            plot!(plt, _symlog.(atmos.flux_cdct), arr_P, label="Conduct", lw=lw*1.2, lc=col_o, ls=:solid, linealpha=la)
        end

        # Latent heating
        if incl_latent
            plot!(plt, _symlog.(atmos.flux_l), arr_P, label="Latent", lw=lw*1.2, lc=col_p, ls=:solid, linealpha=la)
        end

        # Deep heating
        if incl_deep
            plot!(plt, _symlog.(atmos.flux_deep), arr_P, label="Deep", lw=lw*1.2, lc=col_d, ls=:solid, linealpha=la)
        end

        # Sensible heat
        scatter!(plt, [_symlog(atmos.flux_sens)], [arr_P[end]], markershape=:utriangle, markercolor=col_r, label="Sensible")

        # Total flux
        plot!(plt, _symlog.(atmos.flux_tot), arr_P, label="Total", lw=lw, lc=col_t, ls=:solid, linealpha=la)

        # Overplot convection and condensation mask
        #    by indicating it with scatter points of the corresponding colour
        for i in 1:atmos.nlev_c
            if atmos.mask_c[i]
                scatter!(plt, [-0.2], [arr_P[i]], opacity=0.9, markersize=2, msw=0.5, color=col_c, label="")
            end
            if atmos.mask_l[i]
                scatter!(plt, [0.2], [arr_P[i]], opacity=0.9, markersize=2, msw=0.5, color=col_p, label="")
            end
        end

        # Plot current surface pressure and original
        @_plt_pboa
        @_plt_poboa

        # Labels
        annotate!(plt, xlims[1]/2.0, ylims[1]*1.4, text("Downward", :black, :center, 9))
        annotate!(plt, xlims[2]/2.0, ylims[1]*1.4, text("Upward"  , :black, :center, 9))

        # Finalise + save
        xlabel!(plt, "log Unsigned Flux [W m⁻²]")
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
    **Plot emission spectrum at the TOA**

    Arguments:
    - `atmos::atmosphere.Atmos_t`    atmosphere object
    - `fname::String`               filename to save the plot (if empty, does not save)
    """
    function plot_emission(atmos::atmosphere.Atmos_t, fname::String)

        # Check that we have data
        if !(atmos.is_out_lw && atmos.is_out_sw)
            @warn "Cannot plot emission spectrum because radiances have not been calculated"
            return
        end

        # Get emission spectrum data
        xe::Array{Float64, 1} = zeros(Float64, atmos.nbands)
        yt::Array{Float64, 1} = zeros(Float64, atmos.nbands)
        yl::Array{Float64, 1} = zeros(Float64, atmos.nbands)
        ys::Array{Float64, 1} = zeros(Float64, atmos.nbands)
        ye::Array{Float64, 1} = zeros(Float64, atmos.nbands)
        wd::Array{Float64, 1} = zeros(Float64, atmos.nbands)
        yp::Array{Float64, 1} = zeros(Float64, atmos.nbands)

        # band widths
        wd = atmos.bands_wid * 1e9 # convert to nm

        # band centres
        xe = atmos.bands_cen * 1e9 # convert to nm

        # TOA upward spectral flux [erg s-1 cm-2 nm-1]
        @. yl = atmos.band_u_lw[1, :] / wd * 1000.0 # converted to erg s-1 cm-2 nm-1
        @. ys = atmos.band_u_sw[1, :] / wd * 1000.0 # converted to erg s-1 cm-2 nm-1
        @. yt = yl + ys

        # surface upward spectral flux
        @. ye = (atmos.band_u_lw[end, :] + atmos.band_u_sw[end, :]) / wd * 1000.0

        # Get planck function values
        @. yp = phys.evaluate_planck(xe, atmos.tmp_surf) * 1000.0

        # Make plot
        plt = plot(size=(600,400); plt_default...)

        plot!(plt, xe, yp, label=L"Planck @ $T_s$",  color="green")
        plot!(plt, xe, ye, label="Surface LW+SW",    color="green", ls=:dash)

        plot!(plt, xe, ys, lw=0.9, label="Planetary SW",    color="blue")
        plot!(plt, xe, yl, lw=0.9, label="Planetary LW",    color="red" )
        plot!(plt, xe, yt, lw=0.5, label="Planetary LW+SW", color="black")

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
    **Plot contribution function at different bands.**

    The contribution function is plotted with one line (vs pressure) per spectral band.

    Arguments:
    - `atmos::atmosphere.Atmos_t`    atmosphere object
    - `fname::String`               filename to save the plot (if empty, does not save)
    - `size_x::Int64`              width of the plot in pixels
    - `size_y::Int64`              height of the plot in pixels
    - `cf_min::Float64`            minimum contribution function value to plot (log10 units)
    """
    function plot_contfunc1(atmos::atmosphere.Atmos_t, fname::String;
                                    size_x::Int64=size_x_default, size_y::Int64=size_y_default,
                                    cf_min::Float64=1e-6)

        # Check that we have data
        if !atmos.is_out_lw
            @warn "Cannot plot contrib func because radiances have not been calculated"
            return
        end

        # Make plot
        plt = plot(legend=:bottomleft, size=(size_x, size_y); plt_default...)
        x_min::Float64 = log10(cf_min)
        x_max::Float64 = x_min + 0.5

        # Define arrays
        cff::Array{Float64, 1} = zeros(Float64, atmos.nlev_c) # log10 cont func
        prs::Array{Float64, 1} = zeros(Float64, atmos.nlev_c) # pressure [bar]
        @. prs = atmos.p * 1.0e-5

        # Band limits
        wl_min  = 0.1 * 1e-6 # 100 nm
        wl_imin = findmin(abs.(atmos.bands_cen .- wl_min))[2]
        wl_max  = 150 * 1e-6 # 150 um
        wl_imax = findmin(abs.(atmos.bands_cen .- wl_max))[2]

        # reversed?
        if wl_imin > wl_imax
            wl_imin, wl_imax = wl_imax, wl_imin
        end

        # plot statistical contributions at each level, from bands in given range
        for i in 1:atmos.nlev_c
            cff[i] = log10(max(maximum(atmos.contfunc_band[i,wl_imin:wl_imax]), cf_min))
        end
        x_max = max(x_max, maximum(cff))
        plot!(plt, cff, prs, c=:black, label="Maximum", ls=:solid)

        for i in 1:atmos.nlev_c
            cff[i] = log10(max(mean(atmos.contfunc_band[i,wl_imin:wl_imax]), cf_min))
        end
        plot!(plt, cff, prs, c=:black, label="Mean", ls=:dash)

        for i in 1:atmos.nlev_c
            cff[i] = log10(max(median(atmos.contfunc_band[i,wl_imin:wl_imax]), cf_min))
        end
        plot!(plt, cff, prs, c=:black, label="Median", ls=:dot)

        # plot per-band contributions [um] as their own lines
        for wl_tgt in Float64[1.0, 5.0, 10.0, 15.0]
            # find nearest band
            wl_tgt *= 1e-6
            iband = findmin(abs.(atmos.bands_cen .- wl_tgt))[2]
            wl_i = atmos.bands_cen[iband] * 1e6

            # get contribution function
            for i in 1:atmos.nlev_c
                cff[i] = log10(max(atmos.contfunc_band[i,iband], cf_min))
            end

            # plot
            plot!(plt, cff, prs, label=@sprintf("%.1f μm",wl_i))
        end

        xlabel!(plt, "log₁₀ Contribution function")
        xaxis!(plt, xlims=(x_min+0.05, x_max+0.1))

        ylims  = (prs[1], prs[end])
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
    **Plot normalised contribution function (per band)**

    The data displayed in this plot are fine, but the x-axis ticks are labelled
    incorrectly by the plotting library. I don't know why this is.

    Arguments:
    - `atmos::atmosphere.Atmos_t`    atmosphere object
    - `fname::String`               filename to save the plot (if empty, does not save)
    """
    function plot_contfunc2(atmos::atmosphere.Atmos_t, fname::String)

        # Check that we have data
        if !atmos.is_out_lw
            @warn "Cannot plot contribution func because radiances have not been calculated"
            return
        end

        @warn "Contribution func 2D colormesh x-axis is incorrect!"

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

        # z value - log contribution function, normalised
        cf_min = 1.0e-9
        for ba in 1:atmos.nbands
            if reversed
                br = atmos.nbands - ba + 1
            else
                br = ba
            end
            for i in 1:atmos.nlev_c
                z[i,br] = max(atmos.contfunc_band[i,ba],cf_min)
            end
        end
        z /= maximum(z)
        z[:] = log10.(z[:])

        # Make plot
        plt = plot(colorbar_title="log " * L"\widehat {cf}(\lambda, p)"; plt_default...)

        heatmap!(plt, x,y,z, c=:devon, label="")

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
    **Plot spectral albedo (ratio of SW_UP to SW_DN)**

    Arguments:
    - `atmos::atmosphere.Atmos_t`    atmosphere object
    - `fname::String`               filename to save the plot (if empty, does not save)
    - `size_x::Int64`              width of the plot in pixels
    - `size_y::Int64`              height of the plot in pixels
    """
    function plot_albedo(atmos::atmosphere.Atmos_t, fname::String;
                            size_x::Int64=size_x_default, size_y::Int64=size_y_default)

        # Check that we have data
        if !(atmos.is_out_lw && atmos.is_out_sw)
            @warn "Cannot plot spectral albedo because radiances have not been calculated"
            return
        end

        # spectral albedo [percentage]
        y::Array{Float64, 1} = zeros(Float64, atmos.nbands)
        @. y = 100.0 * atmos.band_u_sw[1, :]/atmos.band_d_sw[1, :]

        # Make plot
        ylims  = (0.0, 100.0)
        plt = plot(ylims=ylims, size=(size_x, size_y); plt_default...)

        plot!(plt, atmos.bands_cen*1e9, y, color="black", label="")

        xlims  = (200.0, 1500.0)
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
    **Combined multi-panel plot used for tracking behaviour of the solver at runtime.**
    """
    function combined(plt_pt, plt_fl, plt_mr, plt_ra, info::String, fname::String;
                        size_x::Int64=800, size_y::Int64=700)

        # plt_info = plot(legend=false, showaxis=false, grid=false)
        # annotate!(plt_info, (0.02, 0.7, text(info, family="Courier", :black, :left, 10)))

        plt = plot(plt_pt, plt_fl, plt_mr, plt_ra,
                        plot_title=info,
                        layout=(2,2), size=(size_x, size_y); plt_default...)

        if !isempty(fname)
            savefig(plt, fname)
        end
        return plt
    end

    """
    **Globe plot of multi-column atmosphere climate and energy balance.**

    Left panel: temperature profiles, coloured by column longitude.
    Right panel: heating rate profiles, coloured by column longitude.

    Arguments:
    - `globe::multicol.Globe_t`     globe object
    - `fname::String`               filename to save the plot (if empty, does not save)
    - `size_x::Int64`              width of the plot in pixels
    - `size_y::Int64`              height of the plot in pixels
    """
    function plot_globe(globe::multicol.Globe_t, fname::String;
        size_x::Int64=size_x_default, size_y::Int64=size_y_default)

        atmos = globe.atmos_wrk

        # Plot config
        tlim = Float64[0.0, 200.0]
        hlim = Float64[-10.0, 10.0]
        ylims  = _get_ylims(atmos)
        yticks = _get_yticks(atmos)

        # Create plot
        plt1 = plot(size=(size_x/2, size_y), ylims=ylims, yticks=yticks,
                        legend=:topright; plt_default...)
        plt2 = plot(size=(size_x/2, size_y), ylims=ylims, yticks=(yticks,[]),
                        legend=false; plt_default...)

        # Zero line
        vline!(plt2, [0.0], lw=0.4, lc="black", label="")

        # Plot profiles
        cmap = cgrad(:batlow, globe.ncol, categorical=true, rev=true)
        for i in 1:globe.ncol

            c = cmap[i]
            y = globe.atmos_arr[i].p * 1e-5 # pressure -> bar

            # temperature
            plot!(plt1, globe.atmos_arr[i].tmp, y, color=c, lw=lw, linealpha=la,
                    label=@sprintf("Lo=%4.1f°", globe.lons_arr[i]))

            # store min, max temperatures
            tlim[1] = min(tlim[1], minimum(globe.atmos_arr[i].tmpl)+15.0)
            tlim[2] = max(tlim[2], maximum(globe.atmos_arr[i].tmpl)+15.0)

            # heating rate
            hr = _symlog.(globe.atmos_arr[i].heating_rate) # transformed
            plot!(plt2, hr, y, color=c, lw=lw, linealpha=la, label="")

            # store max heating rate
            hlim[2] = max(hlim[2], maximum(abs.(hr)) * 1.2)
        end
        hlim[1] = -hlim[2]

        # decorate axes
        xaxis!(plt1, xlims=tlim, xlabel="Temperature [K]")
        xaxis!(plt2, xlims=hlim, xlabel="log Unsigned Heating [K day⁻¹]")
        ylabel!(plt1, "Pressure [bar]")

        annotate!(plt2, hlim[1]/2.0, ylims[1]*1.4, text("Cooling", :black, :center, 9))
        annotate!(plt2, hlim[2]/2.0, ylims[1]*1.4, text("Heating", :black, :center, 9))

        # combine into multi-panel plot
        plt = plot(plt1, plt2, layout=(1,2), size=(size_x, size_y),
                        left_margin=[1*Plots.mm -4*Plots.mm]; plt_default...)
        yflip!(plt)
        yaxis!(plt, yscale=:log10)

        # surface pressure
        @_plt_pboa
        @_plt_poboa

        if !isempty(fname)
            savefig(plt, fname)
        end
        return plt
    end

    """
    **Convert still-frame images into an animation**

    Uses FFMPEG provided by Julia library.

    Arguments:
    - `output_dir::String`          folder in which to save the animation
    - `frames_dir::String`          folder containing the frames

    Optional arguments:
    - `output_fmt::String`          output file format: mp4, mov, ...
    - `frames_fmt::String`          input frame format: png, jpg, ...
    - `duration::Float64`           animation duration [seconds]
    """
    function animate(output_dir::String, frames_dir::String;
                            output_fmt::String="mp4",
                            frames_fmt::String="png",
                            duration::Float64=12.0)

        # Find output files
        frames = glob("*.$frames_fmt",frames_dir)
        nframes::Int64 = length(frames)

        # Create animation
        if nframes < 1
            @warn "Cannot create animation; no frames found in $frames_dir"
        else
            fps = Float64(nframes)/duration
            @ffmpeg_env run(`$(FFMPEG.ffmpeg) -loglevel quiet -framerate $fps -pattern_type glob -i "$frames_dir/*.$frames_fmt" -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" -y $output_dir/animation.$output_fmt`)
        end

        return nothing
    end

    """
    Plot jacobian matrix
    """
    function jacobian(b::Array{Float64,2}, fname::String;
                            perturb::Array{Bool,1}=Bool[],
                            size_x::Int64=size_x_default, size_y::Int64=size_y_default)

        lim::Float64 = maximum(abs.(b))     # colourbar limits
        l::Int64 = length(perturb)            # show perturbed levels?

        plt = plot(size=(size_x, size_y),
                    title="∂r/∂x [W m⁻² K⁻¹]",
                    clim=(-lim,lim), yflip=true; plt_default...)

        # show jacobian
        heatmap!(plt, b, color=:RdBu, label="")

        # show perturbed levels
        if l > 0
            scatter!(plt, collect(1:l)[perturb], ones(Float64, l)[perturb]*(l+1),
                        color=:black,markershape=:utriangle, label="")
        end

        xlabel!(plt, "Cell-centre index")
        ylabel!(plt, "Cell-centre index")

        if !isempty(fname)
            savefig(plt, fname)
        end

        return plt
    end

end # end module plotting

