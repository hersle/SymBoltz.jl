using Documenter, SymBoltz

makedocs(
    sitename = "SymBoltz",
    pages = [
        "index.md",
        "getting_started.md",
        "automatic_differentiation.md",
        "extended_models.md",
        "parameter_fitting.md",
        "components.md",
        "conventions.md"
    ]
)
deploydocs(repo = "github.com/hersle/SymBoltz.jl.git")
