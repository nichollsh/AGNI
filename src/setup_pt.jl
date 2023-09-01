# Contains functions for setting-up the temperature profile analytically

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module setup_pt

    include("phys.jl")

    # Read atmosphere T(p) from a file (DOES NOT SET ATMOS STRUCT DIRECTLY) 
    function readcsv(fpath)

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

        # Interpolate to cell-centres
        p[:]   .= 0.5 .* (pl[2:end]   + pl[1:end-1]  )
        tmp[:] .= 0.5 .* (tmpl[2:end] + tmpl[1:end-1])

        # Return dict 
        output = Dict([
                        ("p_surf", pl[end]),
                        ("T_surf", tmpl[end]), 
                        ("p", p),
                        ("pl",pl),
                        ("tmp",tmp),
                        ("tmpl",tmpl)
                      ])
        return output
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
