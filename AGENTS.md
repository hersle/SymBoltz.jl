# SymBoltz.jl: agent guide

SymBoltz.jl is a Julia package that symbolically builds (using ModelingToolkit) and numerically solves (using OrdinaryDiffEq) the linear Einstein-Boltzmann equations in cosmology (like CLASS/CAMB). It is approximation-free (no TCA/UFA/RSA) and differentiable via ForwardDiff.

## Code execution

- Read and run Julia code through an active Julia MCP server, so a single Julia REPL stays alive and code isn't recompiled from scratch every time.
- The REPL preserves state between calls (variables, loaded packages, etc.), so check existing state before redoing setup.
- Precompilation can be necessary and takes a while. When running e.g. `using SymBoltz` (or anything else that triggers compilation), use a timeout of ~10 minutes.
- Run one command at a time. Do not batch them together like `cmd1; cmd2; cmd3; ...` as this may time out due to precompilation.
- If Revise is loaded, file changes are picked up automatically and you do not have to restart the REPL after editing files. Do not manually `include()` source files in the same session that desyncs Revise.
- If the REPL gets into a bad state, restart it rather than fighting it.

## Standard workflow

```julia
using SymBoltz, Plots
M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)

ks = 10.0 .^ range(-3, 2, length=50)
sol = solve(prob, ks)

Ps = spectrum_matter(sol, ks)
plot(log10.(ks), log10.(Ps))

ls = 25:25:2500
jl = SphericalBesselCache(ls)
Dls = spectrum_cmb(:TT, prob, jl; normalization = :Dl)
plot(log10.(jl.l), log10.(Dls))
```

## Repository layout

- `src/SymBoltz.jl`: module entry point, `include` all files and  `export` public functionality.
- `src/solve.jl`: `CosmologyProblem`/`CosmologySolution` types and the `solve`/`solvebg`/`solvept` pipeline (solving background/thermodynamics, perturbations; shooting method, etc.).
- `src/parameters.jl`: canonical parameter sets, e.g. `parameters_Planck18`.
- `src/models/`: symbolic model components (`System` builders), one file per physical sector:
    - `metric.jl`: submodels for the spacetime metric (conformal Newtonian gauge metric) and related variables/equations.
    - `gravity.jl`: submodels for gravitational theories (e.g. general relativity, Brans-Dicke, ...).
    - `generic_species.jl`, `baryons.jl`, `photons.jl`, `neutrinos.jl`, `dark_matter.jl`, `dark_energy.jl`: submodels for different species (e.g. for generic equations of state, baryons, photons, recombination physics, massless and massive neutrinos, dark matter and dark energy).
    - `inflation.jl`: submodel for inflationary power spectrum.
    - `cosmologies.jl`: full models (`ΛCDM`, `RMΛ`, `QCDM`, `w0waCDM`, `BDΛCDM`, ...) created by composing submodels together.
- `src/observables/`: defines ways to compute various observables from `CosmologyProblem`/`CosmologySolution`.
    - `distances.jl`: distance measures, luminosity distance, etc.
    - `fourier.jl`: $k$-space observables, matter power spectrum, arbitrary source functions $S(τ, k)$.
    - `angular.jl`: $l$-space observables, angular CMB power spectra, line-of-sight integration.
- `src/utils.jl`: symbolic-system helper functions, spline utilities, grid builders.
- `test/runtests.jl`: the test suite as a single `@testset`-organized file, to be run as `using Pkg; Pkg.test()` or by including the relevant tests directly in the REPL.
- `docs/src/`: documentation source files for Documenter.jl.

## Creating new physical models

- If possible, create new physical models as incomplete submodels that will be composed with other submodels and place them under `src/models/`.
- See existing models under `src/models/` and the unstructured ΛCDM model in `docs/src/LCDM.md` for examples on how to build models.
- Always create a model as a ModelingToolkit `System(eqs, τ, vars, pars; initial_conditions, initialization_eqs, constraints, guesses)` with the necessary arguments.
    - `eqs` is an array of symbolic equations for both derivative equations and observed equations.
    - `τ` is the symbolic conformal time variable defined inside the `SymBoltz` module.
    - `vars` and `pars` are symbolic variables and parameters created with `@variables` and `@parameters`, respectively.
    - `initial_conditions` is an array of `x => y` mappings that is substituted during initialization.
    - `initialization_eqs` is an array of `x ~ y` equations that is solved during initialization.
    - `constraints` is an array of `x ~ y` equations that the shooting method will solve the background for repeatedly to determine an equal number of shooting variables.
    - `guesses` is an array of `x => y` approximate guesses for variables that are solved for during initialization or in the shooting method.
