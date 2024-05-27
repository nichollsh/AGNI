using Documenter

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
        "Introduction"  => ["index.md"],
        "Examples"      => ["examples.md"],
        "Model"         => ["model.md"],
        "Usage"         => ["usage.md"],
        "Examples"      => ["examples.md"],
        # "Subsection" => [
        #     ...
        # ]
    ]
)

deploydocs(
    repo = "github.com/nichollsh/AGNI.git",
    push_preview=true,
    devbranch="main"
)

