using Documenter
using DocumenterCitations
using DocumenterPages
using DocumenterTools: Themes
using AGNI

# dirs
ASSETS_DIR = joinpath(@__DIR__, "src", "assets")

# metadata
header::String = ""
footer::String = "Copyright © 2023-Present Harrison Nicholls. AGNI source code is available under [Apache-2.0](https://apache.org/licenses/LICENSE-2.0) license. Documentation and assets are available under the [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/) license."
description::String = "AGNI is an open-source model for simulating extreme atmospheres on rocky planets and exoplanets. It simulates 1D radiative-convective-chemical atmosphere profiles and is designed to be fast, flexible, and user-friendly. AGNI is developed by Harrison Nicholls and is a part of the PROTEUS framework."
authors::String = "Harrison Nicholls"
sitename::String = "AGNI"

# following https://github.com/JuliaMusic/JuliaMusic_documentation.jl/blob/master/docs/make.jl
# combine style and defs files into single scss files for compilation
for w in ("light",) # "dark")
    style = read(joinpath(ASSETS_DIR, "style.scss"), String)
    theme = read(joinpath(ASSETS_DIR, "$(w)defs.scss"), String)
    write(joinpath(ASSETS_DIR, "$(w).scss"), style*"\n"*theme)
end

# dark theme is duplicate of light theme
cp(joinpath(ASSETS_DIR, "light.scss"), joinpath(ASSETS_DIR, "dark.scss"), force=true)

# compile styles into scss files
Themes.compile(joinpath(ASSETS_DIR, "light.scss"),
                joinpath(ASSETS_DIR, "themes/documenter-light.css"))
Themes.compile(joinpath(ASSETS_DIR, "dark.scss"),
                joinpath(ASSETS_DIR, "themes/documenter-dark.css"))

format = Documenter.HTML(   edit_link = nothing,
                            collapselevel = 1,
                            size_threshold=Int(1e6),
                            description = description,
                            footer = footer,
                            prettyurls = get(ENV, "CI", nothing) == "true",
                            assets = [
                                # HTML content
                                Documenter.HTMLWriter.RawHTMLHeadContent(header),

                                # local assets
                                "assets/style.css",
                                "assets/logo.ico",

                                # remote assets
                                asset("https://fonts.googleapis.com/css?family=Inter:400&family=JetBrains+Mono:400&family=Lato", class=:css),
                            ]
    )

bib = CitationBibliography(
    joinpath(ASSETS_DIR, "refs.bib"),
    style=:numeric  # :authoryear
)

makedocs(
    sitename = sitename,
    authors = authors,
    format=format,
    pages = [
        "Home" => "index.md",

        PageNode("How-to guides" => "howto/index.md", [
            "Getting started"  => "howto/getting_started.md",
            "Obtaining input data"        => "howto/data.md",
            "Configuring AGNI"            => "howto/configure.md",
            "Using FastChem"              => "howto/fastchem.md",
            "Accessing AGNI from Python"  => "howto/python.md",
            "Running a grid of models"    => "howto/grid.md",
            "Line-by-line radiative transfer" => "howto/lbl.md",
            "Troubleshooting"             => "howto/troubleshoot.md",
            "Contributing to AGNI"        => "howto/contribute.md",
            ],
        ),

        PageNode("Tutorials" => "tutorials/index.md", [
            "Your first calculation"  => "tutorials/01_nosolve.md",
            "Radiative-convective solution"    => "tutorials/02_rce.md",
            "Aerosol radiative properties"       => "tutorials/03_aerosol.md",
            "Steam runaway greenhouse"      => "tutorials/04_runaway.md",
            ],
        ),

        PageNode("Explanation" => "explanation/index.md", [
            "Model description" => "explanation/model.md",
            "Convection" => "explanation/model_convection.md",
            "Radiation" => "explanation/model_radiation.md",
            "Height and gravity" => "explanation/model_height.md",
            "Thermodynamics" => "explanation/model_thermodynamics.md",
            "Sensible heating" => "explanation/model_sensible.md",
            "Advective terms" => "explanation/model_advection.md",
            "Bibliography"       => "explanation/references.md",
            ],
        ),

        PageNode("Reference" => "reference/index.md", [
            "Configuration reference"   => "reference/configuration.md",
            "Solver and output flags"   => "reference/solver_flags.md",
            "API reference"             => "reference/api.md",
            ],
        ),

        "Community" => "thanks.md",

        "Other PROTEUS modules" => "ecosystem.md",

    ],
    plugins=[bib],
)

deploydocs(
    repo = "github.com/nichollsh/AGNI.git",
    push_preview=true
)

