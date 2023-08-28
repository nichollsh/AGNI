# Contains functions for setting-up the temperature profile analytically

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module setup_pt

    include("phys.jl")

    function isothermal!(atmos, set_tmp)
        if !atmos.is_param
            error("Atmosphere parameters not set")
        end 

        atmos.tmpl[:] .= set_tmp 
        atmos.tmp[:] .= set_tmp

    end # end isothermal

    function dry_adiabat!(atmos)
        if !(atmos.is_alloc && atmos.is_param) 
            error("Atmosphere is not setup")
        end 

        # Calculate cell-centre values
        for idx in 1:atmos.nlev_c
            atmos.tmp[idx] = atmos.tstar * ( atmos.p[idx] / atmos.pl[end] ) ^ ( phys.R_gas / atmos.layer_cp[idx] )
        end

        # Interpolate to cell-edge values 
        atmos.tmpl[1]   = atmos.tmp[1]
        atmos.tmpl[end] = atmos.tmp[end]
        for idx in 2:atmos.nlev_l-1
            atmos.tmpl[idx] = 0.5 * (atmos.tmp[idx-1] + atmos.tmp[idx])
        end 
    end # end dry_adiabat
     
end 
