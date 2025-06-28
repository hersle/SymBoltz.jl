# Solving models

First create a symbolic cosmological model:
```@example sol
using SymBoltz
M = SymBoltz.ΛCDM()
nothing # hide
```

## Creating the problem

Once the symbolic [cosmological model](@ref "Models") `M` has been constructed, it can be turned into a numerical problem:
For example:
```@example sol
pars = SymBoltz.parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
```

```@docs
CosmologyProblem
```

## Updating the parameters

Constructing a `CosmologyProblem` is an **expensive** operation that compiles all the symbolics down to numerics.
It is not necessary to repeat this just to update parameter values.
To do so, use the function `parameter_updater` that returns a function that quickly creates new problems with updated parameter values:

```@example sol
probmaker = SymBoltz.parameter_updater(prob, [M.g.h, M.c.Ω₀]) # fast factory function
prob = probmaker([0.70, 0.27]) # create updated problem
```

```@docs
parameter_updater
```

## Solving the problem

The (updated) problem can now be solved for some wavenumbers:
```@example sol
using Unitful, UnitfulAstro
ks = 10 .^ range(-5, +1, length=100) / u"Mpc"
sol = solve(prob, ks)
```

```@docs
solve(prob::CosmologyProblem, ks::AbstractArray)
```

## Accessing the solution

The returned solution `sol` can be conveniently accessed to obtain any variable `y` of the model `M`:

- `sol(y, τ)` returns the *background* variable(s) $y(τ)$ as a function of conformal time(s) $τ$. It interpolates between time points using the ODE solver's custom-tailored interpolator.
- `sol(y, τ, k)` returns the *perturbation* variable(s) $y(τ,k)$ as a function of the wavenumber(s) $k$ and conformal time(s) $τ$. It also interpolates linearly between the logarithms of the wavenumbers passed to `solve`.

Note that `y` can be any symbolic variables in the model `M`, and even expressions thereof.
*Unknown* variables are part of the state vector integrated by the ODE solver, and are returned directly from its solution.
*Observed* variables or expressions are functions of the unknowns, and are automatically calculated from the equations that define them in the symbolic model.
For example:

```@example sol
# TODO: document callable solution when this is fixed: https://github.com/JuliaDocs/Documenter.jl/issues/558 # hide
τs = sol[M.τ] # get time points used in the background solution
ks = [1e-3, 1e-2, 1e-1, 1e0] / u"Mpc" # wavenumbers
as = sol(M.g.a, τs) # scale factors
Ωms = sol((M.b.ρ + M.c.ρ) / M.G.ρ, τs) # matter-to-total density ratios
κs = sol(M.b.rec.κ, τs) # optical depths
Φs = sol(M.g.Φ, τs, ks) # metric potentials
Φs_over_Ψs = sol(M.g.Φ / M.g.Ψ, τs, ks) # ratio between metric potentials
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
p3 = plot(sol, log10(M.g.a), M.g.Φ / M.g.Ψ, ks[1:3]) # exclude last k, where Φ and Ψ cross 0
plot(p1, p2, p3, layout=(3, 1), size=(600, 800))
```

More examples are shown on the [models page](@ref "Models").

## Solve background and perturbations directly

For lower-level control, you can solve the background and perturbations separately:
```@docs
solvebg
solvept
```

## Choice of ODE solver

In principle, models can be solved with any [DifferentialEquations.jl ODE solver](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).
But most cosmological models have very stiff Einstein-Boltzmann equations that can only be solved by implicit solvers, while explicit solvers usually fail.
For the [standard ΛCDM model](@ref "Standard ΛCDM"), some good solvers are:

- `Rodas4P`: Slow for large systems. Very accurate. Handles extreme stiffness. Default background solver.
- `Rodas5P`: Slow for large systems. Very accurate. Handles severe stiffness.
- `KenCarp4` (and `KenCarp47`): Fast. Handles medium stiffness. Default perturbation solver.
- `Kvaerno5`: Behaves similar to `KenCarp4`. Slightly more accurate. Slightly slower.
- `TRBDF2`: Very fast. Decent accuracy. Handles severe stiffness.

See also the [solver benchmarks](@ref "Performance and benchmarks").
