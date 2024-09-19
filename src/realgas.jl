# Contains routines for loading and working with real-gas equations of state

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module realgas

    import ..phys

    using ScatteredInterpolation  # 2D interpolation
    using LoggingExtras
    using DelimitedFiles

    mutable struct Aqua_t

        # Axis
        npts_t::Int
        npts_p::Int

        # Array limits
        lim_prs::Array{Float64}
        lim_tmp::Array{Float64}

        # Interpolators
        itp_rho       # Interpolator for density [kg/m^3]
        itp_gad       # Interpolator for adiabatic gradient (dlog(T)/dlog(P))_S
        itp_mmw       # Interpolator for mmw [kg/mol]

        # Struct
        Aqua_t() = new()
    end

    function read_aqua(fpath::String)::Aqua_t

        # Read the file
        @debug "Reading AQUA table at '$fpath'"
        data::Array{Float64,2} = readdlm(fpath, Float64; header=false, skipstart=19)

        # Parse the columns
        #    input
        arr_prs::Array{Float64} = data[:,1]
        arr_tmp::Array{Float64} = data[:,2]
        #    output
        arr_rho::Array{Float64} = data[:,3]
        arr_gad::Array{Float64} = data[:,4]
        arr_mmw::Array{Float64} = data[:,8]

        # Create struct
        aqua = Aqua_t()

        # Dimensions are hardcoded inside the file
        aqua.npts_p = 1093
        aqua.npts_t = 301

        # Check dimensions vs data
        if length(arr_prs) != aqua.npts_p*aqua.npts_t
            @error "AQUA lookup table has invalid dimensions"
        end

        # Limits
        aqua.lim_prs = [minimum(arr_prs), maximum(arr_prs)]
        aqua.lim_tmp = [minimum(arr_tmp), maximum(arr_tmp)]

        # Create interpolators
        eval_pts = vcat(arr_prs', arr_tmp')
        aqua.itp_rho = interpolate(NearestNeighbor(), eval_pts, arr_rho)
        aqua.itp_gad = interpolate(NearestNeighbor(), eval_pts, arr_gad)
        aqua.itp_mmw = interpolate(NearestNeighbor(), eval_pts, arr_mmw)

        return aqua
    end

    # Evaluate density from temperature and pressure
    function eval_aqua_rho(aqua::Aqua_t, prs::Float64, tmp::Float64)::Float64
        prs = max(min(prs, aqua.lim_prs[2]), aqua.lim_prs[1])
        tmp = max(min(tmp, aqua.lim_tmp[2]), aqua.lim_tmp[1])
        return aqua.itp_rho(tmp,prs)
    end

end
