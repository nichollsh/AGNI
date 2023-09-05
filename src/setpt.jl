# Contains functions for setting-up the temperature profile analytically

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module setpt

    include("phys.jl")

    using PCHIPInterpolation

    # Read atmosphere T(p) from a file 
    function csv!(atmos, fpath)

        # Check file exists
        if !isfile(fpath)
            error("File '$fpath' does not exist")
        end 

        # Read file 
        content = readlines(fpath)

        # Parse file once, to get number of level edges 
        nlev_e = 0  
        for l in content
            if isempty(l) || (l[1] == '#') || (l[1] == '\n')
                continue 
            end
            nlev_e += 1
        end

        # Validate 
        if nlev_e < 3
            error("Csv file contains too few levels (contains $nlev_e edge values)")
        end
        nlev_c = nlev_e - 1

        # Allocate T and P arrays 
        tmpl = zeros(Float64,nlev_e)
        tmp  = zeros(Float64,nlev_c)
        pl   = zeros(Float64,nlev_e)
        p    = zeros(Float64,nlev_c)

        # Parse file again, storing data this time
        idx = 1
        for l in content
            # Skip
            if isempty(l) || (l[1] == '#') || (l[1] == '\n')
                continue 
            end

            # Read
            lsplit = split(strip(l,[' ','#','\t','\n']),',')
            p_val = parse(Float64,lsplit[1])  # Pressure [Pa]
            t_val = parse(Float64,lsplit[2])  # Temperature [K]

            # Validate 
            if p_val <= 0.0
                error("Negative pressure(s) in csv file")
            end 
            if t_val <= 0.0
                error("Negative temperature(s) in csv file")
            end 

            # Store
            pl[idx] = p_val 
            tmpl[idx] = t_val

            # Iterate 
            idx += 1
        end 

        # Check if arrays are flipped 
        if pl[1] > pl[2]  
            pl = reverse(pl)
            tmpl = reverse(tmpl)
        end 

        # Check that pressure is monotonic
        for i in 1:nlev_c
            if pl[i] > pl[i+1]
                error("Pressure is not monotonic in csv file")
            end
        end

        # Set pressure grid
        atmos.p_boa = pl[end]
        atmos.p_toa = pl[1]
        atmosphere.generate_pgrid!(atmos)

        # Set temperatures
        atmos.tstar = tmpl[end]

        itp = Interpolator(pl, tmpl) # Cell edges 
        atmos.tmpl[:] .= itp.(atmos.pl)

        itp = Interpolator(p, tmp) # Cell centres 
        atmos.tmp[:] .= itp.(atmos.p)
    end

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

    # Set pure 'con' atmosphere to phase curve of 'con' when it enters condensible region
    function condensing!(atmos, con::String)
        if !(atmos.is_alloc && atmos.is_param) 
            error("Atmosphere is not setup")
        end 
        if !(con in keys(phys.lookup_L_vap))
            error("Invalid condensible $con")
        end 

        # Get properties
        L = phys.lookup_L_vap[con]
        R = phys.R_gas / phys.lookup_mmw[con]
        p0 = phys.lookup_P_trip[con]
        T0 = phys.lookup_T_trip[con]

        # Check surface pressure (should not be supersaturated)
        tsurf = atmos.tmpl[end]
        psat = p0 * exp( L/R * (1/T0 - 1/tsurf)   )
        if atmos.pl[end] > psat 
            atmos.p_boa = psat 
            atmosphere.generate_pgrid!(atmos)
        end

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
