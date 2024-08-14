# Creating extended models

This tutorial shows how to create extended (beyond-ΛCDM) models.

As a simple example, we will replace the cosmological constant Λ with equation of state $P(t) / \rho(t) = w(t) = -1$ in the pure ΛCDM model
with the w₀wₐ-parametrization from [Chevallier, Polarski (2000)](https://arxiv.org/abs/gr-qc/0009008) and [Linder (2002)](https://arxiv.org/abs/astro-ph/0208512) with equation of state
```math
w(t) = w_0 + w_a (1 - a(t)).
```
With this equation of state, solving the continuity equation
```math
\frac{\mathrm{d} \rho(t)}{\mathrm{d}t} = -3 H(t) \, (ρ(t) + P(t)) = -3 H(t) \, \rho(t) \, (1 + w(t))
```
with the ansatz $\rho(t) = \rho(t_0) \, a(t)^m \exp(n (1 - a(t)))$
gives $m = -3 (1 + w_0 + w_a)$, $n = -3 w_a$ and hence
```math
\rho(t) = \rho(t_0) \, a(t)^{1 + w_0 + w_a} \exp(w_a (1 - a(t))).
```
TODO: use t for cosmic time, τ for conformal time

## Create the symbolic component

SymBoltz.jl represents each component of the Einstein-Boltzmann equations as a partial/incomplete [`ModelingToolkit.jl`](https://docs.sciml.ai/ModelingToolkit/) `ODESystem`.
These are effectively "chunks" of logically related parameters, variables, equations and initial conditions
that are composed together into a full/complete cosmological model.
Let us create a component for the w₀wₐ-parametrization:
```@example 1
using ModelingToolkit # must be loaded to create custom components
using SymBoltz: t, ϵ # load conformal time and perturbation book-keeper

# Create a function that creates a symbolic w0wa dark energy component
function w0wa(g; kwargs...)
    # 1. Create parameters
    pars = @parameters w0 wa ρ0 Ω0

    # 2. Create variables
    vars = @variables ρ(t) P(t) w(t) δ(t) σ(t)

    # 3. Specify background equations (~ means equality in ModelingToolkit)
    eqs0 = [
        w ~ w0 + wa * (1 - g.a) # equation of state
        ρ ~ ρ0 * g.a^(-3 * (1 + w0 + wa)) * exp(-3 * wa * (1 - g.a)) # energy density
        P ~ w * ρ # pressure
    ] .|> SymBoltz.O(ϵ^0) # O(ϵ⁰) multiplies all equations by 1 (no effect), marking them as background equations

    # 4. Specify perturbation (O(ϵ^0)) equations # TODO: include
    eqs1 = [
        δ ~ 0 # energy overdensity
        σ ~ 0 # shear stress
    ] .|> SymBoltz.O(ϵ^1) # O(ϵ^1) multiplies all equations by ϵ, marking them as perturbation equations

    # 5. Specify relationships between parameters
    defaults = [
        ρ0 => 3/8π * Ω0
    ]

    # 6. Pack everything into an ODE system
    return ODESystem([eqs0; eqs1], t, vars, pars; defaults, kwargs...)
end

# TODO: should illustrate initialization_eqs behavior
```

## Create the standard and extended models

We can now create both the standard ΛCDM model,
and then recreate it with the cosmological constant replaced by the w₀wₐ-component to construct the extended w₀wₐCDM model:
```@example 1
# TODO: start with M1 from the very top, then add M2 later
using SymBoltz
M1 = ΛCDM(name = :ΛCDM)

X = w0wa(M1.g; name = :X)
M2 = ΛCDM(Λ = X, name = :w0waCDM)
```

## Solve and compare the models

Now set some parameters and solve both models:
```@example 1
θ1 = [
    M1.γ.Ω0 => 5e-5
    M1.c.Ω0 => 0.27
    M1.b.Ω0 => 0.05
    M1.ν.Neff => 3.0
    M1.g.h => 0.7
    M1.b.rec.Yp => 0.25
]
θ2 = [
    θ1; # extend previous parameter list
    M2.X.w0 => -0.9
    M2.X.wa => 0.2
]
prob1, prob2 = CosmologyProblem(M1), CosmologyProblem(M2)
sol1, sol2 = solve(prob1, θ1), solve(prob2, θ2)
```

Let us compare $H_2 / H_1$ as a function of the scale factor $a$:
```@example 1
# get common extrema for two iterators
function common_extrema(itr1, itr2)
    min1, max1 = extrema(itr1)
    min2, max2 = extrema(itr2)
    return max(min1, min2), min(max1, max2)
end
t1s, t2s = sol1[SymBoltz.t], sol2[SymBoltz.t]
ts = exp.(range(log.(common_extrema(t1s, t2s))..., length=1000))
lga1s, lga2s = sol1(ts, log10(M1.g.a)), sol2(ts, log10(M2.g.a))
H1s, H2s = sol1(ts, M1.g.H), sol2(ts, M2.g.H)

# evaluate at common scale factors
# TODO: define convenience method for doing spline conversion
using DataInterpolations
lgas = range(-3, 0, length=500)
H1s = CubicSpline(H1s, lga1s; extrapolate=true)(lgas)
H2s = CubicSpline(H2s, lga2s; extrapolate=true)(lgas)

using Plots
plot(lgas, H2s ./ H1s; xlabel = "lg(a)", ylabel = "H₂ / H₁", label = nothing)
```
