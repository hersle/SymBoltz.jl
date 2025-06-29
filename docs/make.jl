using Documenter, SymBoltz

ENV["JULIA_DEBUG"] = "Documenter" # make Documenter.jl output more verbose

using Latexify
set_default(starred = true) # equation numbering looks bad

using Base.Threads
println("Julia is running with $(nthreads()) threads")

using LinearAlgebra: BLAS
BLAS.set_num_threads(1)

makedocs(
    sitename = "SymBoltz",
    pages = [
        "index.md",
        "Tutorials" => [
            "getting_started.md",
            "automatic_differentiation.md",
            "extended_models.md",
            "parameter_fitting.md",
            "forecasting.md"
        ],
        "Building models" => [
            "components.md",
            "models.md"
        ],
        "solve.md",
        "plot.md",
        "observables.md",
        "comparison.md",
        "benchmarks.md",
    ],
    pagesonly = true, # only process files listed in pages (useful for testing)
    linkcheck = true, warnonly = [:linkcheck, :cross_references], # warn if documentation has dead links
)
deploydocs(
    repo = "github.com/hersle/SymBoltz.jl.git",
    forcepush = true # don't keep gh-pages history
)
