module style

    import ..formulae: count_atoms

    # Table of pre-defined colors for plotting
    const _lookup_colour::Dict{String, String} = Dict([
        # basic gases
        ("H2O", "#027FB1" ),
        ("CO2", "#D24901" ),
        ("H2" , "#008C01" ),
        ("CH4", "#C720DD" ),
        ("CO" , "#D1AC02" ),
        ("N2" , "#870036" ),
        ("NH3", "#675200" ),
        ("S2" , "#FF8FA1" ),
        ("SO2", "#00008B" ),

        # volatile elements
        ("H"  , "#0000aa"),
        ("D"  , "#0000ff"),
        ("C"  , "#ff0000"),
        ("O"  , "#00cc00"),
        ("N"  , "#ffaa00"),
        ("S" ,  "#ff22ff"),
        ("P" ,  "#33ccff"),
        ("He" , "#30FF71" ),
        ("Ar",  "#FF007F"),
        ("F" ,  "#00FFFF"),
        ("Cl",  "#00FF00"),
        ("Br",  "#FF00FF"),

        # (semi)refractory elements
        ("Na",   "#1144ff"),
        ("Si",   "#aa2277"),
        ("Ti",   "#779922"),
        ("V",    "#555555"),
        ("Mg",   "#996633"),
        ("K",    "#bbbbee"),
        ("Fe",   "#aa8888"),
        ("Li",   "#ffaaaa"),
        ("Rb",   "#ddaacc"),
        ("Cs",   "#ccaaee"),
        ("Ca",   "#ccffdd"),
        ("Al",   "#ff7711"),
        ("Cr",   "#77aa77"),
    ])

    """
    Convert formula to pretty unicode string
    """
    function pretty_name(gas::String)::String
        out::String = ""
        for c in gas
            if isnumeric(c)
                d = parse(Int, string(c))
                out *= Char(parse(Int,"208$d", base=16))
            else
                out *= c
            end
        end
        return out
    end
    export pretty_name

    """
    Generate a colour hex code from a molecular formula
    """
    function pretty_colour(gas::String)::String
        # Defined
        if gas in keys(_lookup_colour)
            return _lookup_colour[gas]
        end

        # Else, generate colour from atoms
        atoms = count_atoms(gas)
        r::Float64 = 0.0
        g::Float64 = 0.0
        b::Float64 = 0.0
        for e in keys(atoms)
            r += parse(Int,_lookup_colour[e][2:3],base=16)*atoms[e]
            g += parse(Int,_lookup_colour[e][4:5],base=16)*atoms[e]
            b += parse(Int,_lookup_colour[e][6:7],base=16)*atoms[e]
        end
        m::Float64 = max(r,g,b)

        # prevents the colour getting too close to white
        if r+g+b > 705
            m *= 255.0/235.0
        end

        # convert to hex code
        out::String = "#"
        out *= string(floor(Int,255 * r/m),base=16,pad=2)
        out *= string(floor(Int,255 * g/m),base=16,pad=2)
        out *= string(floor(Int,255 * b/m),base=16,pad=2)
        return out
    end
    export pretty_colour
end
