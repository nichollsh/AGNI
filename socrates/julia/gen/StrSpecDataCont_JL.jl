# Autogenerated from src/radiance_core/def_spectrum.F90
# svn revision 1226


    mutable struct StrSpecDataCont{P}
        parent_handle::P
    end


    # this is used to show values in the REPL and when using IJulia
    function Base.show(io::IO, m::MIME"text/plain", handle::StrSpecDataCont)
        println(io, handle)
        dump_properties(io, handle)                
    end


    function Base.propertynames(handle::StrSpecDataCont, private::Bool=false)
        names = [
            :n_band_continuum,
            :index_continuum,
            :index_water,
            :i_scale_fnc_cont,
            :k_cont,
            :scale_cont,
            :p_ref_cont,
            :t_ref_cont,
            :k_cont_ses,
            :k_h2oc,
        ]

        return names
    end


    function Base.getproperty(handle::StrSpecDataCont, field::Symbol)


        cptr = getfield(getfield(handle, :parent_handle), :cptr)
        cptr != Ptr{Cvoid}() || error("invalid handle (null cptr)")            

        if field in (
            :index_water,
        )
            val = Ref{Cint}()
            field_ok = ccall(
                (:PS_StrSpecDataCont_get_integer, libSOCRATES_C),
                Cuchar,
                (Ptr{Cvoid}, Cstring, Ref{Cint}),
                cptr, String(field), val
            )
            Bool(field_ok) || error("StrSpecDataCont integer field $field not present - coding error")
            return val[]
        elseif field in (
            :n_band_continuum,
            :index_continuum,
            :i_scale_fnc_cont,
        )
            loc = Ref{Ptr{Cint}}()
            ndim = Ref{Cint}(10)
            dims = zeros(Cint, ndim[])
            lbounds = zeros(Cint, ndim[])
            field_ok = ccall(
                (:PS_StrSpecDataCont_get_integer_array, libSOCRATES_C),
                Cuchar, 
                (Ptr{Cvoid}, Cstring, Ref{Ptr{Cint}}, Ref{Cint}, Ref{Cint}, Ref{Cint}),
                cptr, String(field), loc, dims, ndim, lbounds
            )
            Bool(field_ok) || error("StrSpecDataCont Cint array field $field not present - coding error")
            if loc[] == C_NULL
                return nothing
            else
                a = unsafe_wrap(Array, loc[], Tuple(dims[1:ndim[]]), own=false)
                oa = OffsetArray(a, (lbounds[1:ndim[]] .- 1)...)
                return oa
            end
        elseif field in (
            :k_cont,
            :scale_cont,
            :p_ref_cont,
            :t_ref_cont,
            :k_cont_ses,
            :k_h2oc,
        )
            loc = Ref{Ptr{Cdouble}}()
            ndim = Ref{Cint}(10)
            dims = zeros(Cint, ndim[])
            lbounds = zeros(Cint, ndim[])
            field_ok = ccall(
                (:PS_StrSpecDataCont_get_real_array, libSOCRATES_C),
                Cuchar, 
                (Ptr{Cvoid}, Cstring, Ref{Ptr{Cdouble}}, Ref{Cint}, Ref{Cint}, Ref{Cint}),
                cptr, String(field), loc, dims, ndim, lbounds
            )
            Bool(field_ok) || error("StrSpecDataCont Cdouble array field $field not present - coding error")
            if loc[] == C_NULL
                return nothing
            else
                a = unsafe_wrap(Array, loc[], Tuple(dims[1:ndim[]]), own=false)
                oa = OffsetArray(a, (lbounds[1:ndim[]] .- 1)...)
                return oa
            end
        else
            return getfield(handle, field)
        end    
           
    end


    function Base.setproperty!(handle::StrSpecDataCont, field::Symbol, val)
              

        cptr = getfield(getfield(handle, :parent_handle), :cptr)
        cptr != Ptr{Cvoid}() || error("invalid handle (null cptr)")            

        if field in (
            :index_water,
        )
                    
            field_ok = ccall(
                (:PS_StrSpecDataCont_set_integer, libSOCRATES_C),
                Cuchar,
                (Ptr{Cvoid}, Cstring, Cint),
                cptr, String(field), val
            )
            Bool(field_ok) || error("StrSpecDataCont integer field $field not present - coding error")
            return val
        else
            error("type StrSpecDataCont has no writeable field $field")    
        end    
           
    end
