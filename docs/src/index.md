# Introduction

SymBoltz.jl is a package for solving the linearized Einstein-Boltzmann system of equations. It is designed to be simple to use, extend and integrate with the wider Julia and SciML ecosystems for scientific analysis.

## Features

- Represents the Einstein-Boltzmann system symbolically and solves it numerically (using [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit)).
- Facilitates development of extended models by representing each component (e.g. theory of gravity and particle species) as modular blocks that are joined into a full model.
- Does not rely on approximations (like the tight coupling, ultra-relativistic fluid and radiation streaming approximations), but solves the resulting stiff equations with implicit integrators (using [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs)).
- Enables sensitivity analysis with automatic differentiation of any output quantity (using [ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl)).
- TODO: Accelerates the solution over GPUs (using [DiffEqGPU.jl](https://docs.sciml.ai/DiffEqGPU)).
- Also interoperates with other packages from the wider Julia and SciML ecosystem.

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
