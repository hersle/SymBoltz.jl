# Introduction

SymBoltz.jl is a package for solving the Einstein–Boltzmann equations in cosmology.
It can solve the background and linear perturbations, and compute derived observables such as supernova distances and matter and CMB power spectra.
It is easy to use and extend, and integrates with other Julia packages for scientific analysis.

Compared to traditional codes like [CAMB](https://camb.info/) and [CLASS](http://class-code.net/), SymBoltz offers three main features:

- **Symbolic-numeric:** Represents the Einstein–Boltzmann equations symbolically and solves them numerically, and makes it easier to implement modified cosmological models (using [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/)).
- **Approximation-free:** Integrates the full stiff equations with implicit solvers (using [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/)). Remains as fast as approximation-based codes (e.g. the [TCA, UFA and RSA](https://arxiv.org/abs/1104.2933)) due to analytical and sparse Jacobians.
- **Differentiable:** Computes accurate gradients of any output with respect to any input parameters with automatic differentiation (using [ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl/)).

## Installation and usage

[Install Julia](https://julialang.org/downloads/) and launch it with `julia -tauto`.
Then install and load SymBoltz:

```julia
using SymBoltz
```

This prompts to install SymBoltz if it is uninstalled.
The option `-tauto` enables multithreading to parallelize internal computations.
See examples in the sidebar on how to use SymBoltz, in particular the [Getting started](@ref "Getting started") example.

!!! danger "Avoid restarting Julia!"
    Do not use SymBoltz by restarting `julia -tauto myscript.jl` repeatedly!
    This unnecessarily recompiles code and is slow with Julia's just-in-time compilation.
    Instead start `julia -tauto` once and run commands interactively, with `include("myscript.jl")` or with [IJulia.jl](https://github.com/JuliaLang/IJulia.jl) notebooks.
    See [Julia workflow tips](https://docs.julialang.org/en/v1/manual/workflow-tips/) for more.

## Contributing

Contributions to SymBoltz are welcome!
To modify the package source code, install SymBoltz in development mode:
```julia
using Pkg; Pkg.develop("SymBoltz")
```
Please feel free to open [issues](https://github.com/hersle/SymBoltz.jl/issues) and [pull requests](https://github.com/hersle/SymBoltz.jl/pulls) in SymBoltz repository.

!!! tip "Smoother development workflow with Revise"
    Use [Revise.jl](https://timholy.github.io/Revise.jl/) to automatically reload changes to source files without restarting Julia.
