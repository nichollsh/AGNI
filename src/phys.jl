# Contains physical data

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module phys 

    # Sources:
    # - Python files that accompany Ray's book (=> phys.py)
    # - SOCRATES source code
    # - NIST

    # Universal gas constant, J K-1 mol-1
    const R_gas = 8.31446261815324 

    # Stefan-boltzmann constant, W m−2 K−4
    const sigma = 5.670367e-8 

    # Molecule mean molecular weight, kg mol-1
    const lookup_mmw = Dict([
        ("H2O", 1.801530E-02 ), 
        ("CO2", 4.401000E-02 ), 
        ("O3", 4.799820E-02 ), 
        ("N2O", 4.401280E-02 ), 
        ("CO", 2.801060E-02 ), 
        ("CH4", 1.604300E-02 ), 
        ("O2", 3.199880E-02 ), 
        ("NO", 3.000610E-02 ), 
        ("SO2", 6.406280E-02 ), 
        ("NO2", 4.600550E-02 ), 
        ("NH3", 1.703060E-02 ), 
        ("HNO3", 6.301290E-02 ), 
        ("N2", 2.801340E-02 ), 
        ("CFC11", 1.373686E-01 ), 
        ("CFC12", 1.209140E-01 ), 
        ("CFC113", 1.873765E-01 ), 
        ("HCFC22", 8.646892E-02 ), 
        ("HFC125", 1.200223E-01 ), 
        ("HFC134A", 1.020318E-01 ), 
        ("CFC114", 1.709210E-01 ), 
        ("TiO", 6.386600E-02 ), 
        ("VO", 6.694090E-02 ), 
        ("H2", 2.015880E-03 ), 
        ("He", 4.002602E-03 ), 
        ("OCS", 6.007500E-02 ), 
        ("Na", 2.298977E-02 ), 
        ("K", 3.909830E-02 ), 
        ("FeH", 5.685300E-02 ), 
        ("CrH", 5.300400E-02 ), 
        ("Li", 6.941000E-03 ), 
        ("Rb", 8.546780E-02 ), 
        ("Cs", 1.329055E-01 ), 
        ("PH3", 3.399758E-02 ), 
        ("C2H2", 2.603730E-02 ), 
        ("HCN", 2.702530E-02 ), 
        ("H2S", 3.408100E-02 ), 
        ("Ar", 3.994800E-02 ), 
        ("O", 1.599940E-02 ), 
        ("N", 1.400674E-02 ), 
        ("NO3", 6.301280E-02 ), 
        ("N2O5", 1.080104E-01 ), 
        ("HONO", 4.701340E-02 ), 
        ("HO2NO2", 7.901220E-02 ), 
        ("H2O2", 3.401470E-02 ), 
        ("C2H6", 3.006900E-02 ), 
        ("CH3", 1.503450E-02 ), 
        ("H2CO", 3.002600E-02 ), 
        ("HO2", 3.300670E-02 ), 
        ("HDO", 1.902140E-02 ), 
        ("HCl", 3.646100E-02 ), 
        ("HF", 2.000689E-02 )
    ])

    # Molecule heat capacity at constant pressure, J K-1 kg-1
    const lookup_cp = Dict([
        ("H2O", 1.847000E+03), 
        ("CO2", 8.200000E+02), 
        ("CO", 1.040000E+03), 
        ("CH4", 2.195000E+03), 
        ("O2", 9.160000E+02), 
        ("NH3", 2.060000E+03), 
        ("N2", 1.037000E+03), 
        ("H2", 1.423000E+04), 
        ("He", 5.196000E+03)
    ])

    # Critical point temperature, K
    const lookup_T_crit = Dict([
        ("H2O", 6.471000e+02), 
        ("CO2", 3.042000e+02), 
        ("CO", 1.134450e+02), 
        ("CH4", 1.904400e+02), 
        ("O2", 1.545400e+02), 
        ("NH3", 4.055000e+02), 
        ("N2", 1.262000e+02), 
        ("H2", 3.320000e+01), 
        ("He", 5.100000e+00)
    ])

    # Triple point temperature, K
    const lookup_T_trip = Dict([
        ("H2O", 2.731500E+02), 
        ("CO2", 2.165400E+02), 
        ("CO", 6.795000E+01), 
        ("CH4", 9.067000E+01), 
        ("O2", 5.430000E+01), 
        ("NH3", 1.954000E+02), 
        ("N2", 6.314000E+01), 
        ("H2", 1.395000E+01), 
        ("He", 2.170000E+00)
    ])

    # Triple point pressure, Pa
    const lookup_P_trip = Dict([
        ("H2O", 6.110000E+02),
        ("CO2", 5.185000E+05),
        ("CO", 1.530000E+04),
        ("CH4", 1.170000E+04),
        ("O2", 1.500000E+02),
        ("NH3", 6.100000E+03),
        ("N2", 1.253000E+04),
        ("H2", 7.200000E+03),
        ("He", 5.070000E+03)
    ])

    # Latent heat of vapourisation, J kg-1
    const lookup_L_vap= Dict([
        ("H2O", 2.493000E+06), 
        ("CO2", 3.970000E+05), 
        ("CO", 2.142857E+05), 
        ("CH4", 5.360000E+05), 
        ("O2", 2.420000E+05), 
        ("NH3", 1.658000E+06), 
        ("N2", 2.180000E+05), 
        ("H2", 4.540000E+05), 
        ("He", 2.030000E+04)
    ])

    # Get values from thermodynamic property lookup tables 
    function lookup_safe(prop::String,gas::String)

        prop = lowercase(prop)

        # Find table 
        if prop == "l_vap"
            table = lookup_L_vap
        elseif prop == "p_trip"
            table = lookup_P_trip
        elseif prop == "t_trip"
            table = lookup_T_trip
        elseif prop == "t_crit"
            table = lookup_T_crit
        elseif prop == "mmw"
            table = lookup_mmw
        elseif prop == "cp"
            table = lookup_cp
        else 
            error("Invalid thermodynamic property '$thermo'")
        end 

        # Find value 
        if gas in keys(table)
            return table[gas]
        end 

        # Default case
        return 0.0
    end

    # Get saturation pressure at a given temperature
    function calc_Psat(gas::String, T_eval::Float64)

        # Get properties
        L = phys.lookup_safe("L_vap",gas)
        if L < 1e-20
            return 0.0
        end 

        R = R_gas /  phys.lookup_safe("mmw",gas)
        p0 =  phys.lookup_safe("P_trip",gas)
        T0 =  phys.lookup_safe("T_trip",gas)

        # Calculate Psat
        return p0*exp(-(L/R)*(1.0/T_eval - 1.0/T0))
    end

    # Get dew point temperature at a given pressure (derived from AEOLUS)
    function calc_Tdew(gas::String, p_eval::Float64)

        # Get properties
        L = phys.lookup_safe("L_vap",gas)
        if L < 1e-20
            return 0.0
        end 

        R = R_gas /  phys.lookup_safe("mmw",gas)

        # Avoid math error for p = 0
        p = max(p_eval, 1.0e-100)

        if gas == "H2O"
            Tref = 373.15 # K, boiling point of H2O at 1 atm 
            pref = 1e5 # esat('H2O',Tref) returns 121806.3 Pa, should return 1e5    

        elseif gas == "CH4"
            Tref = 148.15 # K, arbitrary point (148.15K,esat(148.15K)=9.66bar) on the L/G coexistence curve of methane 
            pref = calc_Psat(gas,Tref)

        elseif gas == "CO2"
            Tref = 253.0 # K, arbitrary point (253K,esat(253K)=20.9bar) on the coexistence curve of CO2 
            pref = calc_Psat(gas,Tref)

        elseif gas == "CO"
            Tref = 100. # K, arbitrary point (100K,esat(100K)=4.6bar) on the coexistence curve of CO 
            pref = calc_Psat(gas,Tref)

        elseif gas == "N2"
            Tref = 98.15 # K, arbitrary point (98.15K,esat(98.15K)=7.9bar) on the coexistence curve of N2 
            pref = calc_Psat(gas,Tref)

        elseif gas == "O2"
            Tref = 123.15 # K, arbitrary point (123.15K,esat(123.15K)=21.9bar) on the coexistence curve of O2 
            pref = calc_Psat(gas,Tref)

        elseif gas == "H2"
            Tref = 23.15 # K, arbitrary point (23.15K,esat(23.15K)=1.7bar) on the coexistence curve of H2 
            pref = calc_Psat(gas,Tref)

        elseif gas == "He"
            Tref = 4.22 # K, boiling point of He at 1 atm 
            pref = 1e5 # esat('He',Tref) returns 45196 Pa, should return 1e5

        elseif gas == "NH3"
            Tref = 273.15 # K, arbitrary point (273.15K,esat(273.15K)=8.6bar) on the coexistence curve of NH3 
            pref = calc_Psat(gas,Tref)
        
        else 
            return 0.0
        end

        Tsat = Tref/(1.0-(Tref*R/L)*log(p/pref))
        return Tsat 
    end 

end 
