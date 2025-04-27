# Solving models

Once a symbolic [cosmological model](@ref "Models") `M` has been constructed, it can be turned into a numerical problem that can be solved for some given parameters:

```@docs
CosmologyProblem
parameters(prob::CosmologyProblem)
remake(prob::CosmologyProblem, pars::Dict)
solve(prob::CosmologyProblem, ks::AbstractArray)
```

For example:
```@example sol
using SymBoltz, Unitful, UnitfulAstro
M = SymBoltz.ΛCDM()
pars = SymBoltz.parameters_Planck18(M)
prob = CosmologyProblem(M, pars, Dict(M.Λ.Ω₀ => 0.5), [M.g.ℰ ~ 1])
ks = 10 .^ range(-5, +1, length=100) / u"Mpc"
sol = solve(prob, ks)
```

## Accessing the solution

The returned solution `sol` can be conveniently accessed to obtain any variable `y` of the model `M`:

- `sol(t, y)` returns the *background* variable(s) $y(t)$ as a function of conformal time(s) $t$. It interpolates between time points using the ODE solver's custom-tailored interpolator.
- `sol(k, t, y)` returns the *perturbation* variable(s) $y(k,t)$ as a function of the wavenumber(s) $k$ and conformal time(s) $t$. It also interpolates linearly between the logarithms of the wavenumbers passed to `solve`.

Note that `y` can be any symbolic variables in the model `M`, and even expressions thereof.
*Unknown* variables are part of the state vector integrated by the ODE solver are returned directly from its solution.
*Observed* variables or expressions are functions of the unknowns, and are automatically calculated from the equations that define them in the symbolic model.
For example:

```@example sol
# TODO: document callable solution when this is fixed: https://github.com/JuliaDocs/Documenter.jl/issues/558 # hide
ts = sol[M.t] # get time points used in the background solution
ks = [1e-3, 1e-2, 1e-1, 1e0] / u"Mpc" # wavenumbers
as = sol(ts, M.g.a) # scale factors
Ωms = sol(ts, (M.b.ρ + M.c.ρ) / M.G.ρ) # matter-to-total density ratios
κs = sol(ts, M.b.rec.κ) # optical depths
Φs = sol(ks, ts, M.g.Φ) # metric potentials
Φs_over_Ψs = sol(ks, ts, M.g.Φ / M.g.Ψ) # ratio between metric potentials
nothing # hide
```

## Plotting the solution

SymBoltz.jl includes [plot recipes](https://docs.juliaplots.org/latest/recipes/) for easily visualizing the solution.
It works similarly to the solution accessing: call `plot(sol, [wavenumber(s),] x_expr, y_expr)` to plot `y_expr` as a function of `x_expr`.
For example, to plot some of the same quantities that we obtained above:
```@example sol
using Plots
p1 = plot(sol, log10(M.g.a), (M.b.ρ + M.c.ρ) / M.G.ρ)
p2 = plot(sol, log10(M.g.a), log10(abs(M.b.rec.κ)))
p3 = plot(sol, ks[1:3], log10(M.g.a), M.g.Φ / M.g.Ψ) # exclude last k, where Φ and Ψ cross 0
plot(p1, p2, p3, layout=(3, 1), size=(600, 800))
```

## Shooting method parametrization

It is common to parametrize some models not by initial conditions or constant parameters, but by values of variables at some (non-initial) time, like today.
This is exactly like the boundary conditions of a boundary value problem.
SymBoltz.jl supports such parametrizations with the shooting method:

```@docs
shoot
```

Example usage is shown on the [models page](@ref "Models").

## Choice of solver

In principle, models can be solved with any [DifferentialEquations.jl ODE solver](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).
But most cosmological models have very stiff Einstein-Boltzmann equations that can only be solved by implicit solvers, while explicit solvers usually fail.
For the [standard ΛCDM model](@ref "Standard ΛCDM"), some good solvers are:

- `Rodas4P`: Slow for large systems. Very accurate. Handles extreme stiffness. Default background solver.
- `Rodas5P`: Slow for large systems. Very accurate. Handles severe stiffness.
- `KenCarp4` (and `KenCarp47`): Fast. Handles medium stiffness. Default perturbation solver.
- `Kvaerno5`: Behaves similar to `KenCarp4`. Slightly more accurate. Slightly slower.
- `TRBDF2`: Very fast. Decent accuracy. Handles severe stiffness.

See also the [solver benchmarks](@ref "Benchmarks").
