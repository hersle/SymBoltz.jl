# Introduction

SymBoltz.jl is a package for solving the linearized Einstein-Boltzmann system of equations. It is designed to be simple to use, extend and integrate with the wider Julia and SciML ecosystems for scientific analysis.

## Features

- **Symbolic-numeric:** Represents the Einstein-Boltzmann system symbolically and solves it numerically.
- **Extendable:** Facilitates development of extended models by representing each component (e.g. gravitational theory and every particle species) as modular blocks that are compiled into a full model.
- **Approximation-free:**  Does not rely on tight coupling, ultra-relativistic fluid and radiation streaming approximations, but implicitly integrates the full stiff equations with automatically generated Jacobians.
- **Differentiable:** Enables sensitivity analysis with automatic differentiation of any output quantity.
- **Convenient post-processing:** Compute and plot any derived quantity by its expression with no extra code.
- **Spectra:** Compute linear and non-linear matter and CMB power spectra.
- **GPU-accelerated:** Optionally accelerates the solution over GPUs (TODO).

This is made possible by the packages
[ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit),
[DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs),
[ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl),
[MatterPower.jl](https://github.com/komatsu5147/MatterPower.jl),
[DiffEqGPU.jl](https://docs.sciml.ai/DiffEqGPU)
and more.

## Installation

[Install](https://julialang.org/downloads/) and launch `julia` and install SymBoltz.jl in one of these two ways:

### For usage

```julia
using Pkg; Pkg.add("SymBoltz")
```

This installs the latest release of SymBoltz for use *without* permission to modify internals.
You can still build and modify cosmological models as intended from the symbolic interface.
The installation is tracked by Julia's package manager, and you can easily `] update` to new releases.
**This is the recommended method for most users.**

### For usage and development

```julia
using Pkg; Pkg.develop("SymBoltz")
```

This clones the development repository of SymBoltz for use *with* permission to modify any parts of the code.
The installation is not tracked by Julia's package manager, and you are responsible for managing the local repository and pulling remote updates.
See [the package manager documentation](https://pkgdocs.julialang.org/v1/managing-packages/#developing) to learn more.
**This is the recommended method for contributing users.**

!!! tip
    For a smooth development workflow, it is strongly recommended to [install Revise.jl](https://timholy.github.io/Revise.jl/stable/) and load it *before* Symboltz.jl:
    ```julia
    using Pkg; Pkg.add("Revise"); Pkg.dev("SymBoltz")
    using Revise, SymBoltz # Revise first!
    # 1) Modify SymBoltz.jl code
    # 2) Interactively try out new changes
    # 3) Rinse and repeat
    ```
    Revise automatically tracks changes made to the source files and reloads them in the active Julia session, so Julia does not have to be restarted after every change.

## Usage workflow

Like most Julia packages, SymBoltz is intended for interactive use within a running `julia` REPL session, **not by repeatedly executing scripts** like `julia script.jl` from the shell.
See [Julia workflow tips](https://docs.julialang.org/en/v1/manual/workflow-tips/) to learn more.

Make sure to run `julia --threads=auto` with multi-threading to take advantage of internal parallellizations.
