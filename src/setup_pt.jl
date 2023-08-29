# Contains functions for setting-up the temperature profile analytically

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module setup_pt

    include("phys.jl")
    import DifferentialEquations

    # Set atmosphere to be isothermal
    function isothermal!(atmos, set_tmp)
        if !atmos.is_param
            error("Atmosphere parameters not set")
        end 
        atmos.tmpl[:] .= set_tmp 
        atmos.tmp[:] .= set_tmp
    end 

    # Set atmosphere to dry adiabat
    function dry_adiabat!(atmos)
        # Validate input
        if !(atmos.is_alloc && atmos.is_param) 
            error("Atmosphere is not setup")
        end 

        # Calculate cell-centre values
        for idx in 1:atmos.nlev_c
            atmos.tmp[idx] = atmos.tstar * ( atmos.p[idx] / atmos.pl[end] ) ^ ( phys.R_gas / atmos.layer_cp[idx] )
        end
        
        # Calculate cell-edge values
        for idx in 2:atmos.nlev_l-1
            cp = 0.5 * ( atmos.layer_cp[idx-1] + atmos.layer_cp[idx])
            atmos.tmpl[idx] = atmos.tstar * ( atmos.pl[idx] / atmos.pl[end] ) ^ ( phys.R_gas / cp )
        end

        # Calculate top boundary
        dt = atmos.tmp[1]-atmos.tmpl[2]
        dp = atmos.p[1]-atmos.pl[2]
        atmos.tmpl[1] = atmos.tmp[1] + dt/dp * (atmos.pl[1] - atmos.p[1])
    end 

    # Set atmosphere to phase curve of 'con' when it enters condensible region
    function condensing!(atmos, con::String)

        # Validate input
        if !(atmos.is_alloc && atmos.is_param) 
            error("Atmosphere is not setup")
        end 
        if !(con in keys(phys.lookup_L_vap))
            error("Invalid condensible $con")
        end 

        L = phys.lookup_L_vap[con]
        R = phys.R_gas / phys.lookup_mmw[con]
        p0 = phys.lookup_P_trip[con]
        T0 = phys.lookup_T_trip[con]

        # Check if each level is condensing. If it is, place on phase curve
        for i in 1:atmos.nlev_c
            # Cell centre
            p = atmos.p[i]
            Tsat = 1.0/(  1/T0 - R/L * log(p/p0)  )
            if atmos.tmp[i] < Tsat
                atmos.tmp[i] = Tsat
            end
                
            # Cell edge
            p = atmos.pl[i]
            Tsat = 1.0/(  1/T0 - R/L * log(p/p0)  )
            if atmos.tmpl[i] < Tsat
                atmos.tmpl[i] = Tsat
            end
        end

    end 

end 
