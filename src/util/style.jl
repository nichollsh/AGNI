# This file is part of AGNI. License is GPL-3.0: https://www.gnu.org/licenses

"""
**Module for defining styles for plotting and printing.**
"""
module style

    import ..formulae: count_atoms

    # Colors for plotting
    const col_r::String = "#c0c0c0"; export col_r # radiation
    const col_n::String = "#000000"; export col_n # net
    const col_c::String = "#6495ed"; export col_c # convection
    const col_t::String = "#ff4400"; export col_t # temperature
    const col_o::String = "#66CD00"; export col_o # conduction
    const col_p::String = "#ecb000"; export col_p # phase change
    const col_d::String = "#8B008B"; export col_d # deepheating
    const col_b::String = "#eeeeee"; export col_b # bandpass

    # Allowed plot file extensions
    const ALLOWED_EXTS::Set{String} = Set(["png", "pdf", "svg"])
    export ALLOWED_EXTS

    # Telescope bandpasses [micron], from PROTEUS plot.py
    const OBSERVER_BANDS = Dict{String, Dict{String, NTuple{2, Float64}}}(
        # https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-instrumentation/miri-filters-and-dispersers
        "MIRI" => Dict(
            "F560W" => (5.054, 6.171),
            "F770W" => (6.581, 8.687),
            "F1000W" => (9.023, 10.891),
            "F1130W" => (10.953, 11.667),
            "F1280W" => (11.588, 14.115),
            "F1500W" => (13.527, 16.64),
            "F1800W" => (16.519, 19.502),
            "F2100W" => (18.477, 23.159),
            "F2550W" => (23.301, 26.733),
            "LRS" => (5.0, 14.0),
        ),
        # https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-instrumentation/nirspec-dispersers-and-filters#
        "NIRSpec" => Dict(
            "F070LP" => (0.70, 1.27),
            "F100LP" => (0.97, 1.84),
            "F170LP" => (1.66, 3.07),
            "F290LP" => (2.87, 5.10),
            "PRISM" => (0.60, 5.30),
        ),
        # https://jwst-docs.stsci.edu/jwst-near-infrared-imager-and-slitless-spectrograph#gsc.tab=0
        "NIRISS" => Dict(
            "SOSS" => (0.6, 2.8),
            "WFSS" => (0.8, 2.2),
            "AMI" => (2.8, 4.8),
        ),
        # https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters
        "NIRCam" => Dict(
            "Short" => (0.6, 2.3),
            "Long" => (2.4, 5.0),
        ),
        # https://www.esa.int/Science_Exploration/Space_Science/Ariel/Ariel_s_instruments
        "ARIEL" => Dict(
            "AIRS0" => (1.95, 3.9),
            "AIRS1" => (3.9, 7.8),
        ),
        # https://en.wikipedia.org/wiki/Infrared_astronomy
        # https://www.eso.org/sci/facilities/paranal/instruments/crires/inst.html
        # https://www.eso.org/sci/facilities/paranal/instruments/gravity/overview.html
        "IR" => Dict(
            "R" => (0.65, 1.0),
            "J" => (1.115, 1.362),
            "H" => (1.423, 1.769),
            "K" => (1.972, 2.624),
            "L" => (2.869, 4.188),
            "M" => (4.6, 5.0),
            "N" => (7.5, 14.5),
            "Q" => (17.0, 25.0),
            "Z" => (28.0, 40.0),
        ),
        # https://link.springer.com/article/10.1007/s10686-020-09660-1
        "PLATO" => Dict(
            "blue" => (0.500, 0.675),
            "red" => (0.675, 1.125),
        ),
        # https://doi.org/10.1051/0004-6361/202140366
        "LIFE" => Dict(
            "LIFE" => (4.0, 18.5),
        ),
        # https://ntrs.nasa.gov/api/citations/20240006497/downloads/HWO%20Engineering%20View%20Status%20Plans%20Opportunities.pdf
        "HWO" => Dict(
            "Coronograph" => (0.4, 1.8),
            "Highres imager" => (0.2, 2.5),
            "Spectrograph" => (0.1, 1.0),
        ),
        # https://www.gemini.edu/instrumentation/maroon-x
        # https://www.gemini.edu/instrumentation/igrins-2
        "GEMINI-N" => Dict(
            "MAROON-X" => (0.5, 0.92),
            "IGRINS-2" => (1.49, 2.46),
        ),
        # https://www.eso.org/sci/facilities/paranal/instruments/sphere.html
        # https://www.eso.org/sci/facilities/paranal/instruments/espresso/overview.html
        "VLT" => Dict(
            "ESPRESSO" => (0.38, 0.788),
            "SPHERE" => (0.95, 2.32),
        ),
        # https://carmenes.caha.es/ext/instrument/index.html
        "CARMENES" => Dict(
            "CARMENES" => (0.520, 1.710),
        ),
        # https://www.tng.iac.es/instruments/harps/
        "HARPS" => Dict(
            "HARPS-N" => (0.383, 0.690),
        ),
        # https://noirlab.edu/public/programs/kitt-peak-national-observatory/wiyn-35m-telescope/neid/
        "NEID" => Dict(
            "NEID" => (0.38, 0.93),
        ),
        # https://elt.eso.org/instrument/
        "ELT" => Dict(
            "HARMONI" => (0.47, 2.45),
            "MICADO" => (0.8, 2.4),
            "ANDES" => (0.40, 1.80),
        ),
    )
    export OBSERVER_BANDS

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
