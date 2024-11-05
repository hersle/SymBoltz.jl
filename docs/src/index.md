# Introduction

SymBoltz.jl is a package for solving the linearized Einstein-Boltzmann system of equations. It is designed to be simple to use, extend and integrate with the wider Julia and SciML ecosystems for scientific analysis.

## Features

- **Symbolic-numeric:** Represents the Einstein-Boltzmann system symbolically and solves it numerically.
- **Extendable:** Facilitates development of extended models by representing each component (e.g. gravitational theory and every particle species) as modular blocks that are joined into a full model.
- **Approximation-free:**  Does not rely on tight coupling, ultra-relativistic fluid and radiation streaming approximations, but implicitly integrates the full stiff equations with automatically generated Jacobians.
- **Differentiable:** Enables sensitivity analysis with automatic differentiation of any output quantity.
- **Convenient post-processing:** Compute and plot any derived quantity by its expression with no extra code.
- **GPU-accelerated:** Optionally accelerates the solution over GPUs (TODO).

This is made possible by the packages
[ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit),
[DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs),
[ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl) and
[DiffEqGPU.jl](https://docs.sciml.ai/DiffEqGPU).

## Installation

[Install](https://julialang.org/downloads/) and launch `julia`, enter the package manager by hitting `]`, and execute one of the commands below.

### For use

```
add https://github.com/hersle/SymBoltz.jl
```

### For use and development

```
dev https://github.com/hersle/SymBoltz.jl
```

## Usage

See the tutorials in the sidebar.
