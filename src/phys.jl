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
    const R_gas::Float64 = 8.314462618 # NIST CODATA

    # Stefan-boltzmann constant, W m−2 K−4
    const sigma::Float64 =  5.670374419e-8 # NIST CODATA

    # Von Karman constant, dimensionless
    const k_vk::Float64 = 0.40  # Hogstrom 1988

    # Planck's constant, J s
    const h_pl::Float64 = 6.62607015e-34  # NIST CODATA

    # Boltzmann constant, J K-1 
    const k_B::Float64 = 1.380649e-23 # NIST CODATA

    # Speed of light, m s-1
    const c_vac::Float64 = 299792458.0 # NIST CODATA

    # Pretty-formatted names
    const lookup_pretty = Dict{String, String}([
        ("H2O",    "H₂O"      ), 
        ("CO2",    "CO₂"      ), 
        ("O3",     "O₃"       ), 
        ("N2O",    "N₂O"      ), 
        ("CO",     "CO"       ), 
        ("CH4",    "CH₄"      ), 
        ("O2",     "O₂"       ), 
        ("NO",     "NO"       ), 
        ("SO2",    "SO₂"      ), 
        ("NO2",    "NO₂"      ), 
        ("NH3",    "NH₃"      ), 
        ("HNO3",   "HNO₃"     ), 
        ("N2",     "N₂"       ), 
        ("CFC11",  "CFC11"    ), 
        ("CFC12",  "CFC12"    ), 
        ("CFC113", "CFC113"   ), 
        ("HCFC22", "HCFC22"   ), 
        ("HFC125", "HFC125"   ), 
        ("HFC134A","HFC134A"  ), 
        ("CFC114", "CFC114"   ), 
        ("TiO",    "TiO"      ), 
        ("VO",     "VO"       ), 
        ("H2",     "H₂"       ), 
        ("He",     "He"       ), 
        ("OCS",    "OCS"      ), 
        ("Na",     "Na"       ), 
        ("K",      "K"        ), 
        ("FeH",    "FeH"      ), 
        ("CrH",    "CrH"      ), 
        ("Li",     "Li"       ), 
        ("Rb",     "Rb"       ), 
        ("Cs",     "Cs"       ), 
        ("PH3",    "PH₃"      ), 
        ("C2H2",   "C₂H₂"     ), 
        ("HCN",    "HCN"      ), 
        ("H2S",    "H₂S"      ), 
        ("Ar",     "Ar"       ), 
        ("O",      "O"        ), 
        ("N",      "N"        ), 
        ("NO3",    "NO₃"      ), 
        ("N2O5",   "N₂O₅"     ), 
        ("HONO",   "HONO"     ), 
        ("HO2NO2", "HO₂NO₂"   ), 
        ("H2O2",   "H₂O₂"     ), 
        ("C2H6",   "C₂H₆"     ), 
        ("CH3",    "CH3"      ), 
        ("H2CO",   "H₂CO"     ), 
        ("HO2",    "HO₂"      ), 
        ("HDO",    "HDO"      ), 
        ("HCl",    "HCl"      ), 
        ("HF",     "HF" )
    ])

    # Pre-defined absorber colors 
    const lookup_color = Dict{String, String}([
        ("H2O", "#C720DD" ),
        ("CO2", "#D24901" ),
        ("H2" , "#008C01" ),
        ("CH4", "#027FB1" ),
        ("CO" , "#D1AC02" ),
        ("N2" , "#870036" ),
        ("S"  , "#FF8FA1" ),
        ("O2" , "#00008B" ),
        ("He" , "#30FF71" ),
        ("NH3", "#675200" ),
    ])

    # Molecule mean molecular weight, kg mol-1
    const lookup_mmw = Dict{String, Float64}([
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

    # Molecule heat capacity at constant pressure and temperature, J K-1 kg-1
    const lookup_cp = Dict{String, Float64}([
        ("H2O", 1.847000E+03), 
        ("CO2", 8.200000E+02), 
        ("CO",  1.040000E+03), 
        ("CH4", 2.195000E+03), 
        ("O2",  9.160000E+02), 
        ("NH3", 2.060000E+03), 
        ("N2",  1.037000E+03), 
        ("H2",  1.423000E+04), 
        ("He",  5.196000E+03),
        ("O3",  819.37)
    ])

    # Molecule thermal conductivity at constant pressure and temperature, J K-1 kg-1
    const lookup_kc = Dict{String, Float64}([
        ("H2O", 0.11967),   # at 1173 K   | all values in this dict sourced 
        ("CO2", 0.07396),   # at 1050 K   | from the engineering toolbox website 
        ("CH4", 0.17670),   # at 1000 K   | https://www.engineeringtoolbox.com
    ])

    # Critical point temperature, K
    const lookup_T_crit = Dict{String, Float64}([
        ("H2O", 6.471000e+02), 
        ("CO2", 3.042000e+02), 
        ("CO", 1.134450e+02), 
        ("CH4", 1.904400e+02), 
        ("O2", 1.545400e+02), 
        ("NH3", 4.055000e+02), 
        ("N2", 1.262000e+02), 
        ("H2", 3.320000e+01), 
        ("He", 5.100000e+00),
        ("O3", 261.15)
    ])

    # Triple point temperature, K
    const lookup_T_trip = Dict{String, Float64}([
        ("H2O", 2.731500E+02), 
        ("CO2", 2.165400E+02), 
        ("CO", 6.795000E+01), 
        ("CH4", 9.067000E+01), 
        ("O2", 5.430000E+01), 
        ("NH3", 1.954000E+02), 
        ("N2", 6.314000E+01), 
        ("H2", 1.395000E+01), 
        ("He", 2.170000E+00),
        ("O3", 80.0),

    ])

    # Triple point pressure, Pa
    const lookup_P_trip = Dict{String, Float64}([
        ("H2O", 6.110000E+02),
        ("CO2", 5.185000E+05),
        ("CO", 1.530000E+04),
        ("CH4", 1.170000E+04),
        ("O2", 1.500000E+02),
        ("NH3", 6.100000E+03),
        ("N2", 1.253000E+04),
        ("H2", 7.200000E+03),
        ("He", 5.070000E+03),
        ("O3", 7.346E-1)
    ])

    # Latent heat of vapourisation, J kg-1
    const lookup_L_vap= Dict{String, Float64}([
        ("H2O", 2.493000E+06), 
        ("CO2", 3.970000E+05), 
        ("CO", 2.142857E+05), 
        ("CH4", 5.360000E+05), 
        ("O2", 2.420000E+05), 
        ("NH3", 1.658000E+06), 
        ("N2", 2.180000E+05), 
        ("H2", 4.540000E+05), 
        ("He", 2.030000E+04),
        ("O3", 2.8849E+5)
    ])

    # Shomate formulation for heat capacity, J K-1 kg-1
    function shomate_cp(gas::String, tmp::Float64)::Float64
        # Get molar mass for this molecule 
        mmw::Float64 = 0.0
        if gas in keys(lookup_mmw)
            mmw = lookup_mmw[gas]  # kg mol-1
        else 
            return 0.0
        end 

        # Coefficient base values 
        cA::Float64 = 0.0
        cB::Float64 = 0.0
        cC::Float64 = 0.0
        cD::Float64 = 0.0
        cE::Float64 = 0.0

        # Set coefficients
        if gas == "H2O"
            # NIST Tables | Chase, 1998
            tmp = clamp(tmp, 500.0, 6000.0)
            if tmp < 1700.0 
                cA = 30.09200	
                cB = 6.832514	
                cC = 6.793435	
                cD = -2.534480
                cE = 0.082139	
            else 
                cA = 41.96426
                cB = 8.622053
                cC = -1.499780
                cD = 0.098119
                cE = -11.15764
            end 
        elseif gas == "CO2"
            # NIST Tables | Chase, 1998
            tmp = clamp(tmp, 298.0, 6000.0)
            if tmp < 1200.0 
                cA = 24.99735	
                cB = 55.18696	
                cC = -33.69137
                cD = 7.948387	
                cE = -0.136638
            else 
                cA = 58.16639
                cB = 2.720074
                cC = -0.492289
                cD = 0.038844
                cE = -6.447293
            end 
        elseif gas == "CO"
            # NIST Tables | Chase, 1998
            tmp = clamp(tmp, 298.0, 6000.0)
            if tmp < 1300.0 
                cA = 25.56759	
                cB = 6.096130	
                cC = 4.054656	
                cD = -2.671301
                cE = 0.131021	
            else 
                cA = 35.15070
                cB = 1.300095
                cC = -0.205921
                cD = 0.013550
                cE = -3.282780
            end 
        elseif gas == "CH4"
            # NIST Tables | Chase, 1998
            tmp = clamp(tmp, 298.0, 6000.0)
            if tmp < 1300.0 
                cA = -0.703029
                cB = 108.4773
                cC = -42.52157	
                cD = 5.862788
                cE = 0.678565	
            else 
                cA = 85.81217
                cB = 11.26467
                cC = -2.114146
                cD = 0.138190
                cE = -26.42221
            end 
        elseif gas == "O2"
            # NIST Tables | Chase, 1998
            tmp = clamp(tmp, 100.0, 6000.0)
            if tmp < 700.0 
                cA = 31.32234	
                cB = -20.23531
                cC = 57.86644	
                cD = -36.50624
                cE = -0.007374
            elseif tmp < 2000.0
                cA = 30.03235	
                cB = 8.772972	
                cC = -3.988133
                cD = 0.788313	
                cE = -0.741599
            else 
                cA = 20.91111
                cB = 10.72071
                cC = -2.020498
                cD = 0.146449
                cE = 9.245722
            end 
        elseif gas == "NH3"
            # NIST Tables | Chase, 1998
            tmp = clamp(tmp, 298.0, 6000.0)
            if tmp < 1400.0 
                cA = 19.99563	
                cB = 49.77119	
                cC = -15.37599
                cD = 1.921168	
                cE = 0.189174	
            else 
                cA = 52.02427
                cB = 18.48801
                cC = -3.765128
                cD = 0.248541
                cE = -12.45799
            end 
        elseif gas == "N2"
            # NIST Tables | Chase, 1998
            tmp = clamp(tmp, 100.0, 6000.0)
            if tmp < 500.0 
                cA = 28.98641	
                cB = 1.853978	
                cC = -9.647459
                cD = 16.63537	
                cE = 0.000117	
            elseif tmp < 2000.0
                cA = 19.50583	
                cB = 19.88705	
                cC = -8.598535
                cD = 1.369784	
                cE = 0.527601	
            else 
                cA = 35.51872
                cB = 1.128728
                cC = -0.196103
                cD = 0.014662
                cE = -4.553760
            end
        elseif gas == "H2"
            # NIST Tables | Chase, 1998
            tmp = clamp(tmp, 298.0, 6000.0)
            if tmp < 1000.0 
                cA = 33.066178	
                cB = -11.363417
                cC = 11.432816	
                cD = -2.772874	
                cE = -0.158558	
            elseif tmp < 2500.0
                cA = 18.563083
                cB = 12.257357
                cC = -2.859786
                cD = 0.268238	
                cE = 1.977990	
            else 
                cA = 43.413560
                cB = -4.293079
                cC = 1.272428
                cD = -0.096876
                cE = -20.533862
            end
        elseif gas == "He"
            # NIST Tables | Chase, 1998
            tmp = clamp(tmp, 298.0, 6000.0)
            cA = 20.78603
            cB = 4.850638e-10
            cC = -1.582916e-10
            cD = 1.525102e-11
            cE = 3.196347e-11
        elseif gas == "O3"
            # NIST Tables | Chase, 1998
            tmp = clamp(tmp, 298.0, 6000.0)
            if tmp < 1200.0 
                cA = 21.66157	
                cB = 79.86001	
                cC = -66.02603
                cD = 19.58363	
                cE = -0.079251
            else
                cA = 57.81409
                cB = 0.730941
                cC = -0.039253
                cD = 0.002610
                cE = -3.560367
            end
        else 
            return 0.0
        end 

        # Evalulate heat capacity and convert units to J K-1 kg-1
        t1::Float64 = tmp/1000.0
        t2::Float64 = t1*t1
        return (cA + cB*t1 + cC*t2 + cD*t2*t1 + cE/t2)/mmw 
    end 

    # Temperature dependent latent heat interpolation, J kg-1 
    function interp_Lv(gas::String, tmp::Float64)
        out = 0.0

        # Get molar mass for this molecule 
        mmw = 0.0
        if gas in keys(lookup_mmw)
            mmw = lookup_mmw[gas]  # kg mol-1
        else 
            return 0.0
        end 

        error("Temperature-dependent latent heat not yet implemented")

    end 

    # Get values from thermodynamic property lookup tables 
    function lookup_safe(prop::String,gas::String;tmp::Float64=-1.0)::Float64

        prop = lowercase(prop)

        _table::Dict{String, Float64} = Dict()
        _funct = nothing

        # Find table 
        if prop == "l_vap"
            _table = lookup_L_vap
            # funct = interp_Lv
        elseif prop == "p_trip"
            _table = lookup_P_trip
        elseif prop == "t_trip"
            _table = lookup_T_trip
        elseif prop == "t_crit"
            _table = lookup_T_crit
        elseif prop == "mmw"
            _table = lookup_mmw
        elseif prop == "cp"
            _table = lookup_cp
            _funct = shomate_cp
        elseif prop == "kc"
            _table = lookup_kc
        else 
            error("Invalid thermodynamic property '$prop'")
        end 

        # Try temperature-dependent cases 
        if !isnothing(_funct) && (tmp > 0.1)
            return _funct(gas,tmp)

        # Otherwise, try lookup tables
        else
            if gas in keys(_table)
                return _table[gas]
            end
        end 
         
        # Default case
        return 0.0
    end

    # Get saturation pressure at a given temperature
    function calc_Psat(gas::String, T_eval::Float64)::Float64

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

    # Get dew point temperature at a given pressure
    function calc_Tdew(gas::String, p_eval::Float64)::Float64

        # Get properties
        L::Float64 = phys.lookup_safe("L_vap",gas)
        if L < 1e-20
            return 0.0
        end 

        R::Float64 = R_gas /  phys.lookup_safe("mmw",gas)
        Tref::Float64 = 0.0
        pref::Float64 = 0.0

        # Avoid math error for p = 0
        p::Float64 = max(p_eval, 1.0e-50)

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

        return Tref/(1.0-(Tref*R/L)*log(p/pref)) 
    end 

end 
