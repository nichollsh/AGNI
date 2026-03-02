using Documenter
using DocumenterPages
using AGNI

format = Documenter.HTML(edit_link = "main",
                         prettyurls = get(ENV, "CI", nothing) == "true",
                         assets = [
                             "assets/style.css",
                             "assets/logo.ico",
                        ]
)

makedocs(
    sitename="AGNI",
    format=format,
    pages = [
        "Home" => "index.md",

        PageNode("Tutorials" => "tutorials/index.md", [
            "Getting started"  => "tutorials/getting_started.md",
            "Example outputs"  => "tutorials/examples.md",
            ]
        ),

        PageNode("How-to guides" => "howto/index.md", [
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

        PageNode("Explanation" => "explanation/index.md", [
            "Model description"  => "explanation/model.md",
            "Related codes"      => "explanation/ecosystem.md",
            "Bibliography"       => "explanation/references.md",
            ],
        ),

        PageNode("Reference" => "reference/index.md", [
            "Configuration reference"   => "reference/configuration.md",
            "Solver and output flags"   => "reference/solver_flags.md",
            "API reference"             => "reference/api.md",
            ],
        ),
        "Acknowledgements" => "thanks.md",
    ]
)

deploydocs(
    repo = "github.com/nichollsh/AGNI.git",
    push_preview=true
)

