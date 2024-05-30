using Documenter
# using AGNI

format = Documenter.HTML(edit_link = "main",
                         prettyurls = get(ENV, "CI", nothing) == "true",
                         assets = [
                             joinpath("assets", "style.css"),
                        ]
)

makedocs(
    sitename="AGNI",
    format=format,
    pages = [
        "Home" => "index.md",
        "model.md",
        "setup.md" ,
        "usage.md",
        "examples/index.md",
        "manual/index.md"
    ]
)

deploydocs(
    repo = "github.com/nichollsh/AGNI.git",
    push_preview=true
)

