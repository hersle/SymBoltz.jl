# SymBoltz.jl

[![Version](https://img.shields.io/github/v/release/hersle/SymBoltz.jl?color=000000)](https://github.com/hersle/SymBoltz.jl/releases)
[![Documentation](https://img.shields.io/badge/documentation%20-%20stable%20-%20%232e63b8)](https://hersle.github.io/SymBoltz.jl/stable/)
[![Journal paper](https://img.shields.io/badge/paper%20-%20Astronomy%20&%20Astrophysics%20-%20%23004e7c)](http://dx.doi.org/10.1051/0004-6361/202557450)
[![ArXiv paper](https://img.shields.io/badge/preprint%20-%20arXiv%20-%20%23b31b1b)](http://dx.doi.org/10.1051/0004-6361/202557450)
[![Build status](https://img.shields.io/github/actions/workflow/status/hersle/SymBoltz.jl/ContinuousIntegration.yml)](https://github.com/hersle/SymBoltz.jl/actions/workflows/ContinuousIntegration.yml)
[![License](https://img.shields.io/github/license/hersle/SymBoltz.jl)](https://github.com/hersle/SymBoltz.jl/blob/main/LICENSE)

SymBoltz.jl is a Julia package for solving the linear Einstein-Boltzmann equations in cosmology. It is:

- **Symbolic-numeric:** Models are specified with symbolic equations in a simple, convenient and extensible high-level format, then automatically compiled to fast numerical functions that are solved.
- **Approximation-free:** Full equations are solved at all times without tight-coupling, ultrarelativistic fluid and radiation-streaming approximations (TCA, UFA and RSA) using efficient implicit ODE integrators.
- **Differentiable**: Get derivatives of any output (e.g. power spectra) with respect to any input (e.g. cosmological parameters) using automatic differentiation (exact to numerical precision and does not need step size tuning).

## Installation

[Install Julia](https://julialang.org/install/), launch `julia -tauto` (for multithreading) and run:
```julia
using Pkg; Pkg.add("SymBoltz")
```

## Documentation

Visit the [**documentation**](https://hersle.github.io/SymBoltz.jl/) to get started with tutorials and examples.

## Gallery

<img width="1677" height="1295" alt="bilde" src="https://github.com/user-attachments/assets/72982f26-4a9d-4962-8075-b94462fc711a" />

## Citation

If you use SymBoltz in your work, please cite [this peer-reviewed paper](http://dx.doi.org/10.1051/0004-6361/202557450) (a fully equivalent [arXiv preprint](https://arxiv.org/abs/2509.24740) is also available):
```
@article{SymBoltz,
  title = {{{SymBoltz.jl}}: {{A}} symbolic-numeric, approximation-free, and differentiable linear {{Einstein}}--{{Boltzmann}} solver},
  author = {Herman Sletmoen},
  year = 2026,
  month = mar,
  journal = {Astronomy \& Astrophysics},
  volume = {707},
  pages = {A128},
  doi = {10.1051/0004-6361/202557450},
  url = {http://dx.doi.org/10.1051/0004-6361/202557450},
}
```
