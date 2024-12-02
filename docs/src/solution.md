# Solution handling

Once a [cosmological model](@ref "Models") has been constructed, it can be solved for some given parameters:

```@docs
solve(M::CosmologyModel, pars)
solve(M::CosmologyModel, pars, ks::AbstractArray)
```

For example:
```@example sol
using SymBoltz, Unitful, UnitfulAstro
M = SymBoltz.ΛCDM()
pars = SymBoltz.parameters_Planck18(M)
ks = 10 .^ range(-5, +1, length=100) / u"Mpc"
sol = solve(M, pars, ks)
```

## Accessing the solution

The returned solution `sol` can be conveniently accessed to obtain any variable `y` of the model `M`:

- `sol(t, y)` returns the *background* variable(s) $y(t)$ as a function of conformal time(s) $t$. It interpolates between time points using the ODE solver's custom-tailored interpolator.
- `sol(k, t, y)` returns the *perturbation* variable(s) $y(k,t)$ as a function of the wavenumber(s) $k$ and conformal time(s) $t$. It also interpolates linearly between the logarithms of the wavenumbers passed to `solve`.

Note that `y` can be any symbolic variables in the model `M`, and even expressions thereof.
*Unknown* variables are part of the state vector integrated by the ODE solver are returned directly from its solution,
while *observed* variables and expressions are functions of the unknowns are automatically calculated from their definitions in `M`.
For example:

```@example sol
# TODO: document callable solution when this is fixed: https://github.com/JuliaDocs/Documenter.jl/issues/558 # hide
ts = sol[M.t] # get time points used in the background solution
ts = exp.(range(log.(extrema(ts))..., length=1000)) # 1000 uniformly log-spaced times
ks = [1e-3, 1e-2, 1e-1, 1e0] / u"Mpc" # wavenumbers
as = sol(ts, M.g.a) # scale factors
Ωms = sol(ts, (M.b.ρ + M.c.ρ) / M.G.ρ) # matter-to-total density ratios
τs = sol(ts, M.b.rec.τ) # optical depths
Φs = sol(ks, ts, M.g.Φ) # metric potentials
Φs_over_Ψs = sol(ks, ts, M.g.Φ / M.g.Ψ) # ratio between metric potentials
nothing # hide
```

## Choice of solver

In principle, models can be solved with any [DifferentialEquations.jl ODE solver](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).
But most cosmological models have very stiff Einstein-Boltzmann equations that can only be solved by implicit solvers, while explicit solvers usually fail.
For the [standard ΛCDM model](@ref "Standard ΛCDM"), some good solvers are:

- `Rodas5P`: Slow for large systems. Very accurate. Handles severe stiffness. Default background solver.
- `KenCarp4` (and `KenCarp47`): Fast. Handles medium stiffness. Default perturbation solver.
- `Kvaerno5`: Behaves similar to `KenCarp4`. Slightly more accurate. Slightly slower.
- `TRBDF2`: Very fast. Decent accuracy. Handles severe stiffness.

See also the [solver benchmarks](@ref "Benchmarks").
