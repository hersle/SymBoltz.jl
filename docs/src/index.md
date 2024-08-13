# SymBoltz.jl: linear Einstein-Boltzmann equation solver

SymBoltz.jl is a package for solving the linearized Einstein-Boltzmann system of equations. It is designed to be simple to use, extend and integrate with the wider Julia and SciML ecosystems for scientific analysis.

## Features

- Represents the Einstein-Boltzmann system symbolically and solves it numerically (using [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit)).
- Facilitates development of extended models by representing each component (e.g. theory of gravity and different particle species) as modular blocks of symbolic expressions that are composed to build a custom model.
- Does not rely on approximations (like the tight coupling, ultra-relativistic fluid and radiation streaming approximations), but solves the resulting stiff equations with implicit integrators (using [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs)).
- Enables automatic differentiation of any output quantity (using [ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl)).
- ~~Accelerates the solution over GPUs (using [DiffEqGPU.jl](https://docs.sciml.ai/DiffEqGPU)).~~
- ~~Uses the variable conventions of the [seminal paper by Ma and Bertschinger (1995)](https://arxiv.org/abs/astro-ph/9506072) in most places.~~
- Interoperates with other packages from the wider Julia and SciML ecosystem for parameter fitting, sensitivity analysis, and so on.

## Installation

Launch `julia`, enter the package manager by hitting `]` and enter (one of)
```
add https://github.com/hersle/SymBoltz.jl # for usage
dev https://github.com/hersle/SymBoltz.jl # for development
```

## Usage

TODO: tutorials
