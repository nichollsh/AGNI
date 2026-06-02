using Documenter
using DocumenterPages
using AGNI

format = Documenter.HTML(edit_link = "main",
                         prettyurls = get(ENV, "CI", nothing) == "true",
                         assets = [
                            # local assets
                            "assets/style.css",
                            "assets/logo.ico",

                            # remote assets
                            asset("https://fonts.googleapis.com/css?family=Space+Grotesk:400&family=JetBrains+Mono:400&family=Lato", class=:css),
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
            "Running a full model"    => "tutorials/02_rce.md",
            "Aerosol formation"       => "tutorials/03_aerosol.md",
            "Runaway greenhouse"      => "tutorials/04_runaway.md",
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

