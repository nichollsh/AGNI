# Contains routines for loading and working with real-gas equations of state

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module realgas

    import ..phys

    using Dierckx           # 2D spline interpolation
    using LoggingExtras
    using DelimitedFiles

    mutable struct Aqua_t

        # Axis
        npts_t::Int
        npts_p::Int

        # Array limits
        lim_prs::Array{Float64,2}
        lim_tmp::Array{Float64,2}

        # Interpolators
        itp_rho::Spline2D       # Interpolator for density [kg/m^3]
        itp_gad::Spline2D       # Interpolator for adiabatic gradient (dlog(T)/dlog(P))_S
        itp_mmw::Spline2D       # Interpolator for mmw [kg/mol]

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

        # Reshape arrays
        evl_prs::Array{Float64}     = arr_prs[1:aqua.npts_t:end]
        evl_tmp::Array{Float64}     = arr_tmp[1:1:aqua.npts_t]
        evl_rho::Array{Float64,2}   = reshape(arr_rho, (aqua.npts_p,aqua.npts_t))
        evl_gad::Array{Float64,2}   = reshape(arr_gad, (aqua.npts_p,aqua.npts_t))
        evl_mmw::Array{Float64,2}   = reshape(arr_mmw, (aqua.npts_p,aqua.npts_t))

        # Limits
        aqua.lim_prs = [minimum(evl_prs), maximum(evl_prs)]
        aqua.lim_tmp = [minimum(evl_tmp), maximum(evl_tmp)]

        # Create interpolators
        aqua.itp_rho = Spline2D(evl_prs, evl_tmp, evl_rho)
        aqua.itp_gad = Spline2D(evl_prs, evl_tmp, evl_gad)
        aqua.itp_mmw = Spline2D(evl_prs, evl_tmp, evl_mmw)

        return aqua
    end

    # Evaluate density from temperature and pressure
    function eval_aqua_rho(aqua::Aqua_t, prs::Float64, tmp::Float64)::Float64
        prs = max(min(prs, aqua.lim_prs[2]), aqua.lim_prs[1])
        tmp = max(min(tmp, aqua.lim_tmp[2]), aqua.lim_tmp[1])
        return aqua.itp_rho(tmp,prs)
    end

end
