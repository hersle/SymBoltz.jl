# Solving models

First create a symbolic cosmological model:
```@example sol
using SymBoltz
M = ΛCDM()
nothing # hide
```

## Creating the problem

Once the symbolic [cosmological model](@ref "Cosmologies (full models)") `M` has been constructed, it can be turned into a numerical problem:
For example:
```@example sol
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars; jac = true, sparse = true)
```
The keyword arguments generate a analytical and sparse Jacobian matrix, so solving large perturbation systems is efficient.

```@docs
CosmologyProblem
```

## Updating the parameters

Constructing a `CosmologyProblem` is an **expensive** operation that compiles all the symbolics down to numerics.
It is not necessary to repeat this just to update parameter values.
To do so, use the function `parameter_updater` that returns a function that quickly creates new problems with updated parameter values:

```@example sol
probmaker = parameter_updater(prob, [M.g.h, M.c.Ω₀]) # fast factory function
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
κs = sol(M.b.κ, τs) # optical depths
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
p2 = plot(sol, log10(M.g.a), log10(abs(M.b.κ)))
p3 = plot(sol, log10(M.g.a), M.g.Φ / M.g.Ψ, ks[1:3]) # exclude last k, where Φ and Ψ cross 0
plot(p1, p2, p3, layout=(3, 1), size=(600, 800))
```

More examples are shown on the [models page](@ref "Cosmologies (full models)").

## Shooting method

Some problems require tuning parameters or initial conditions to satisfy constraints at a later time.
This is handled by the shooting method, which uses a rootfinder to repeatedly solve the background for different parameters to find the values that satisfies the constraints.
As a trivial example, we can construct a model where the continuity equation for the cosmological constant is integrated numerically and not analytically:
```@example sol
g = SymBoltz.metric()
Λ = SymBoltz.cosmological_constant(g; analytical = false)
M = ΛCDM(; g, Λ)
equations(background(M.Λ))
```
We specify to shoot for $\rho_\Lambda(\tau_\text{ini})$ and give an initial guess that is used as the starting point in Newton's method.
We also specify that the constraint $H/H₀ = 1$ (in code units) must hold today:
```@example sol
shoot = Dict(M.Λ.ρ => 0.0)
conditions = [M.g.H ~ 1]
prob = CosmologyProblem(M, pars, shoot, conditions)
sol = solve(prob; verbose = true)
@assert sol[M.γ.Ω + M.ν.Ω + M.c.Ω + M.b.Ω + M.h.Ω + M.Λ.Ω][end] ≈ 1 # hide
nothing # hide
```
You can specify any number of shooting variables and conditions, but they must be equal in number to form a well-defined rootfinding problem.
When there is only one shooting variable, we can also use bracketing rootfinders instead of Newton's method.
To do this, replace the scalar guess with an interval:
```@example sol
shoot = Dict(M.Λ.ρ => (0.0, 0.5))
prob = CosmologyProblem(M, pars, shoot, conditions)
sol = solve(prob; verbose = true)
@assert sol[M.γ.Ω + M.ν.Ω + M.c.Ω + M.b.Ω + M.h.Ω + M.Λ.Ω][end] ≈ 1 # hide
nothing # hide
```

## Solve background and perturbations directly

For lower-level control, you can solve the background and perturbations separately:
```@docs
solvebg
solvept
```

## Choice of ODE solver

In principle, models can be solved with any [OrdinaryDiffEq.jl ODE solver](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).
But most cosmological models have very stiff Einstein-Boltzmann equations that can only be solved by implicit solvers, while explicit solvers usually fail.
For the stiff [standard ΛCDM model](@ref "Standard ΛCDM"), we find success with these solvers (from best to worst):

1. **[Rosenbrock methods](https://docs.sciml.ai/DiffEqDocs/latest/api/ordinarydiffeq/semiimplicit/Rosenbrock):** `Rodas5P`, `Rodas4P`, `Rodas5`, `Rodas4`.
2. **[ESDIRK methods](https://docs.sciml.ai/DiffEqDocs/latest/api/ordinarydiffeq/implicit/SDIRK):** `KenCarp4`, `KenCarp47`, `Kvaerno5`, `TRBDF2`.
3. **[BDF methods](https://docs.sciml.ai/DiffEqDocs/latest/api/ordinarydiffeq/implicit/BDF):** `FBDF`, `QNDF`.

See the [solver benchmarks](@ref "Performance and benchmarks") for comparisons between them.
