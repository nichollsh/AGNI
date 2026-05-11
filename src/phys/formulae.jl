module formulae

    # Atom counts in some standard species
    const _lookup_count_atoms::Dict{String, Dict} = Dict([
        "H2O" => Dict("H"=>2, "O"=>1),
        "H2"  => Dict("H"=>2),
        "O2"  => Dict("O"=>2),
        "O3"  => Dict("O"=>3),
        "CO"  => Dict("C"=>1, "O"=>1),
        "CO2" => Dict("C"=>1, "O"=>2),
        "CH4" => Dict("C"=>1, "H"=>4),
        "N2"  => Dict("N"=>2),
        "NH3" => Dict("N"=>1, "H"=>3),
        "SO2" => Dict("S"=>1, "O"=>2),
        "H2S" => Dict("S"=>1, "H"=>2),
        "S2"  => Dict("S"=>2),
        "S8"  => Dict("S"=>8),
    ]
    )


    # Table of gas molecule mean molecular weights, kg mol-1
    const _lookup_mmw::Dict{String, Float64} = Dict([
        # molecules
        ("H2O",     1.801530E-02 ),
        ("HDO",     2.002760E-02 ),
        ("CO2",     4.401000E-02 ),
        ("O3",      4.799820E-02 ),
        ("N2O",     4.401280E-02 ),
        ("CO",      2.801060E-02 ),
        ("CH4",     1.604300E-02 ),
        ("O2",      3.199880E-02 ),
        ("NO",      3.000610E-02 ),
        ("SO2",     6.406280E-02 ),
        ("NO2",     4.600550E-02 ),
        ("NH3",     1.703060E-02 ),
        ("HNO3",    6.301290E-02 ),
        ("N2",      2.801340E-02 ),
        ("TiO",     6.386600E-02 ),
        ("VO",      6.694090E-02 ),
        ("H2",      2.015880E-03 ),
        ("OCS",     6.007500E-02 ),
        ("FeH",     5.685300E-02 ),
        ("CrH",     5.300400E-02 ),
        ("PH3",     3.399758E-02 ),
        ("C2H2",    2.603730E-02 ),
        ("HCN",     2.702530E-02 ),
        ("H2S",     3.408100E-02 ),
        ("NO3",     6.301280E-02 ),
        ("N2O5",    1.080104E-01 ),
        ("HONO",    4.701340E-02 ),
        ("HO2NO2",  7.901220E-02 ),
        ("H2O2",    3.401470E-02 ),
        ("C2H6",    3.006900E-02 ),
        ("CH3",     1.503450E-02 ),
        ("H2CO",    3.002600E-02 ),
        ("HO2",     3.300670E-02 ),
        ("HCl",     3.646100E-02 ),
        ("HF",      2.000689E-02 ),

        # elements from https://iupac.qmul.ac.uk/AtWt/
        ("H",   1.008000000e-03 ),
        ("D",   2.014100000e-03 ), # https://chemlin.org/isotope/hydrogen-2
        ("He",  4.002000000e-03 ),
        ("Li",  6.940000000e-03 ),
        ("Be",  9.012000000e-03 ),
        ("B",   1.081000000e-02 ),
        ("C",   1.201100000e-02 ),
        ("N",   1.400700000e-02 ),
        ("O",   1.599900000e-02 ),
        ("F",   1.899800000e-02 ),
        ("Ne",  2.017970000e-02 ),
        ("Na",  2.298900000e-02 ),
        ("Mg",  2.430500000e-02 ),
        ("Al",  2.698100000e-02 ),
        ("Si",  2.808500000e-02 ),
        ("P",   3.097300000e-02 ),
        ("S",   3.206000000e-02 ),
        ("Cl",  3.545000000e-02 ),
        ("Ar",  3.995000000e-02 ),
        ("K",   3.909830000e-02 ),
        ("Ca",  4.007800000e-02 ),
        ("Sc",  4.495500000e-02 ),
        ("Ti",  4.786700000e-02 ),
        ("V",   5.094150000e-02 ),
        ("Cr",  5.199610000e-02 ),
        ("Mn",  5.493800000e-02 ),
        ("Fe",  5.584500000e-02 ),
        ("Co",  5.893300000e-02 ),
        ("Ni",  5.869340000e-02 ),
        ("Cu",  6.354600000e-02 ),
        ("Zn",  6.538000000e-02 ),
        ("Ga",  6.972300000e-02 ),
        ("Ge",  7.263000000e-02 ),
        ("As",  7.492100000e-02 ),
        ("Se",  7.897100000e-02 ),
        ("Br",  7.990400000e-02 ),
        ("Kr",  8.379800000e-02 ),
        ("Rb",  8.546780000e-02 ),
        ("Sr",  8.762000000e-02 ),
        ("Y",   8.890500000e-02 ),
        ("Zr",  9.122400000e-02 ),
        ("Nb",  9.290600000e-02 ),
        ("Mo",  9.595000000e-02 ),
        ("Ru",  1.010700000e-01 ),
        ("Rh",  1.029050000e-01 ),
        ("Pd",  1.064200000e-01 ),
        ("Ag",  1.078682000e-01 ),
        ("Cd",  1.124140000e-01 ),
        ("In",  1.148180000e-01 ),
        ("Sn",  1.187100000e-01 ),
        ("Sb",  1.217600000e-01 ),
        ("Te",  1.276000000e-01 ),
        ("I",   1.269040000e-01 ),
        ("Xe",  1.312930000e-01 ),
        ("Cs",  1.329050000e-01 ),
        ("Ba",  1.373270000e-01 ),
        ("La",  1.389050000e-01 ),
        ("Ce",  1.401160000e-01 ),
        ("Pr",  1.409070000e-01 ),
        ("Nd",  1.442420000e-01 ),
        ("Sm",  1.503600000e-01 ),
        ("Eu",  1.519640000e-01 ),
        ("Gd",  1.572500000e-01 ),
        ("Tb",  1.589250000e-01 ),
        ("Dy",  1.625000000e-01 ),
        ("Ho",  1.649300000e-01 ),
        ("Er",  1.672590000e-01 ),
        ("Tm",  1.689340000e-01 ),
        ("Yb",  1.730450000e-01 ),
        ("Lu",  1.749668000e-01 ),
        ("Hf",  1.784860000e-01 ),
        ("Ta",  1.809470000e-01 ),
        ("W",   1.838400000e-01 ),
        ("Re",  1.862070000e-01 ),
        ("Os",  1.902300000e-01 ),
        ("Ir",  1.922170000e-01 ),
        ("Pt",  1.950840000e-01 ),
        ("Au",  1.969660000e-01 ),
        ("Hg",  2.005920000e-01 ),
        ("Tl",  2.043800000e-01 ),
        ("Pb",  2.072000000e-01 ),
        ("Bi",  2.089800000e-01 ),
        ("Th",  2.320377000e-01 ),
        ("Pa",  2.310350000e-01 ),
        ("U",   2.380280000e-01),
    ])

    """
    Get number of atoms from formula, returning a dictionary
    """
    function count_atoms(molec::String)::Dict{String,Int64}

        # Pre-defined molecules
        if haskey(_lookup_count_atoms, molec)
            return _lookup_count_atoms[molec]
        end

        # Remove unsafe chars
        m = String(molec)
        for c in ['(',')','[',']','{','}','-','+',' ']
            m = replace(m, c => "")
        end

        # Setup
        out::Dict{String,Int64} = Dict{String,Int64}()
        nchar::Int64 = length(m)
        i::Int64 = 1
        elem::String = ""
        count::Int64=-1
        last::Bool=false

        # Loop through string
        while i <= nchar
            last = (i == nchar)

            # new element
            if isuppercase(m[i])
                count = 0
                elem = string(m[i])
                if !last && islowercase(m[i+1])  # two letter element name
                    elem = elem*string(m[i+1])
                    i += 1
                end
            end

            last = (i == nchar)

            # get count
            if count == 0   # expecting number
                # number of atoms
                if last || isletter(m[i+1]) # got letter => count=1
                    count = 1
                else
                    count = parse(Int, m[i+1])
                end
                # repeated element
                if elem in keys(out)
                    out[elem] += count
                else
                    out[elem] = count
                end
                # reset
                elem = ""
                count = -1
            end
            i += 1
        end

        return out
    end
    export count_atoms

    """
    Check if two gas atom dicts are equivalent
    """
    function same_atoms(d1::Dict, d2::Dict)::Bool

        # check if have same atoms at all
        for k in keys(d1)
            if !(k in keys(d2))
                return false
            end
        end

        # ^^ reverse combination
        for k in keys(d2)
            if !(k in keys(d1))
                return false
            end
        end

        # check counts
        for k in keys(d1)
            if d1[k] != d2[k]
                return false
            end
        end

        # if we haven't returned false so far, then it must be true
        return true
    end
    export same_atoms

    """
    Calculate species mean molecular weight [kg mol-1] from formula or use known value
    """
    function get_mmw(m::String)::Float64

        # already defined?
        if m in keys(_lookup_mmw)
           return _lookup_mmw[m]
        end

        # get atoms
        atoms::Dict{String, Int64} = count_atoms(m)

        # add up atoms
        mmw::Float64 = 0.0
        for k in keys(atoms)
            mmw += _lookup_mmw[k]*atoms[k]
        end

        return mmw
    end
    export get_mmw

end
