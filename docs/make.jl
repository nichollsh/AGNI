using Documenter
using AGNI

format = Documenter.HTML(edit_link = "main",
                         prettyurls = get(ENV, "CI", nothing) == "true",
                         assets = [
                             joinpath("assets", "style.css"),
                             joinpath("assets", "logo.ico")
                        ]
)

makedocs(
    sitename="AGNI",
    format=format,
    pages = [
        "Home" => "index.md",
        "model/index.md",
        "setup.md" ,
        "usage.md",
        "examples/index.md",
        "troubleshooting.md",
        "Software manual" => "manual.md",
        "Related codes" => "ecosystem.md",
        "Contributors" => "contributors.md",
        "Contributing" => "contributing.md",
    ]
)

deploydocs(
    repo = "github.com/nichollsh/AGNI.git",
    push_preview=true
)

