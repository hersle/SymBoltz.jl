using Documenter, SymBoltz

ENV["JULIA_DEBUG"] = "Documenter" # make Documenter.jl output more verbose

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
        "Components and models" => [
            "components.md",
            "models.md"
        ],
        "conventions.md",
        "comparison.md",
    ],
    pagesonly = true, # only process files listed in pages (useful for testing)
)
deploydocs(repo = "github.com/hersle/SymBoltz.jl.git")
