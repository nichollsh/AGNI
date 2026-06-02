using Documenter
using DocumenterPages
using AGNI

# following https://github.com/JuliaMusic/JuliaMusic_documentation.jl/blob/master/docs/make.jl
# combine style and defs files into single scss files for compilation
for w in ("light",) # "dark")
    header = read(joinpath(@__DIR__, "style.scss"), String)
    theme = read(joinpath(@__DIR__, "$(w)defs.scss"), String)
    write(joinpath(@__DIR__, "$(w).scss"), header*"\n"*theme)
end

# dark theme is duplicate of light theme
cp(joinpath(@__DIR__, "light.scss"), joinpath(@__DIR__, "dark.scss"), force=true)

# compile styles into scss files
using DocumenterTools: Themes
Themes.compile(joinpath(@__DIR__, "light.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-light.css"))
Themes.compile(joinpath(@__DIR__, "dark.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"))

format = Documenter.HTML(edit_link = "main",
                         prettyurls = get(ENV, "CI", nothing) == "true",
                         assets = [
                            # local assets
                            "assets/style.css",
                            "assets/logo.ico",

                            # remote assets
                            asset("https://fonts.googleapis.com/css?family=Inter:400&family=JetBrains+Mono:400&family=Lato", class=:css),
                        ]
)

makedocs(
    sitename="AGNI",
    authors = "Harrison Nicholls",
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
            "Model description"  => "explanation/model.md",
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

    ]
)

deploydocs(
    repo = "github.com/nichollsh/AGNI.git",
    push_preview=true
)

