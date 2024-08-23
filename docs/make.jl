using Documenter, SymBoltz

makedocs(
    sitename = "SymBoltz",
    pages = [
        "index.md",
        "Tutorials" => [
            "getting_started.md",
            "automatic_differentiation.md",
            "extended_models.md",
            "parameter_fitting.md"
        ],
        "components.md",
        "conventions.md",
        "comparison.md",
    ],
    pagesonly = true, # only process files listed in pages (useful for testing)
)
deploydocs(repo = "github.com/hersle/SymBoltz.jl.git")
