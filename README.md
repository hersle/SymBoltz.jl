# SymBoltz.jl

SymBoltz.jl is a Julia package for solving the linear Einstein-Boltzmann equations in cosmology. It is:

- **Symbolic-numeric:** Models are specified with symbolic equations in a simple, convenient and extensible high-level format, then automatically compiled to fast numerical functions that are solved.
- **Approximation-free:** Full equations are solved at all times without tight-coupling, ultrarelativistic fluid and radiation-streaming approximations (TCA, UFA and RSA) using efficient implicit ODE integrators.
- **Differentiable**: Get derivatives of any output (e.g. power spectra) with respect to any input (e.g. cosmological parameters) using automatic differentiation (exact to numerical precision and does not need step size tuning).

## Installation

[Install Julia](https://julialang.org/install/), launch `julia` and run:
```julia
using Pkg; Pkg.add("SymBoltz")
```

## Documentation

Visit [**the documentation pages**](https://hersle.github.io/SymBoltz.jl/) to get started with tutorials and examples.

## Gallery

<img width="1677" height="1295" alt="bilde" src="https://github.com/user-attachments/assets/72982f26-4a9d-4962-8075-b94462fc711a" />
