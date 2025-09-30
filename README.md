# SymBoltz.jl

[![Documentation](https://img.shields.io/badge/documentation%20-%20stable%20-%20%232e63b8)](https://hersle.github.io/SymBoltz.jl/stable/)
[![Code paper](https://img.shields.io/badge/paper%20-%20arXiv%20-%20%23b31b1b)](https://arxiv.org/abs/2509.24740)
[![Build status](https://img.shields.io/github/actions/workflow/status/hersle/SymBoltz.jl/ContinuousIntegration.yml)](https://github.com/hersle/SymBoltz.jl/actions/workflows/ContinuousIntegration.yml)
[![License](https://img.shields.io/github/license/hersle/SymBoltz.jl)](https://github.com/hersle/SymBoltz.jl/blob/main/LICENSE)

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

## Citation

If you use SymBoltz in your work, please cite [this paper](https://arxiv.org/abs/2509.24740):
```
@article{SymBoltz,
  title = {{SymBoltz.jl}: a symbolic-numeric, approximation-free and differentiable linear {Einstein-Boltzmann} solver},
  author = {Herman Sletmoen},
  year = {2025},
  journal = {arXiv},
  eprint = {2509.24740},
  archiveprefix = {arXiv},
  primaryclass = {astro-ph.CO},
  doi = {10.48550/arXiv.2509.24740},
  url = {http://arxiv.org/abs/2509.24740},
}
```
