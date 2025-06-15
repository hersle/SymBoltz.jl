# Extending models

This tutorial shows how to create extended (beyond-ΛCDM) models.
As a simple example, we replace the cosmological constant with equation of state ``w(τ) = -1`` by w₀wₐ dark energy with equation of state
```math
w(τ) = \frac{P(τ)}{ρ(τ)} = w_0 + w_a (1 - a(τ)),
```
turning the ΛCDM model into the w₀wₐCDM model,
as suggested by [Chevallier, Polarski (2000)](https://arxiv.org/abs/gr-qc/0009008) and [Linder (2002)](https://arxiv.org/abs/astro-ph/0208512).

SymBoltz represents each component of the Einstein-Boltzmann equations as a [`ModelingToolkit.ODESystem`](https://docs.sciml.ai/ModelingToolkit/stable/systems/ODESystem/#ModelingToolkit.ODESystem).
These are effectively incomplete "chunks" of logically related parameters, variables, equations and initial conditions
that are composed together into a complete cosmological model.
A cosmological model can be modified by adding or replacing components.

## 1. Create the reference model

First, we create the reference model that is to be extended.
We will modify the standard ΛCDM model:
```@example ext
using SymBoltz
M1 = ΛCDM()
hierarchy(M1; describe = true)
```
Everything belonging to the cosmological constant species is stored in the component `M1.Λ`:
```@example ext
equations(M1.Λ)
```
This is what we will replace.

## 2. Create the extended component

The background continuity equation
``ρ̇(τ) = -3 H(τ) \, ρ(τ) \, (1 + w(τ))``
can be solved analytically with the w₀wₐ equation of state and the ansatz
```math
ρ(τ) = ρ(τ₀) \, a(τ)ᵐ \exp(n (1 - a(τ))),
```
giving ``m = -3 (1 + w₀ + wₐ)`` and ``n = -3 wₐ``.
The perturbation equations are
```math
\begin{aligned}
δ̇ &= -(1 + w) (θ - 3 Φ̇) - 3 \frac{ȧ}{a} (c_s^2 - w) δ, \\
θ̇ &= -\frac{ȧ}{a} (1 - 3w) θ - \frac{ẇ}{1+w} θ + \frac{c_s^2}{1+w} k^2 δ - k^2 σ + k^2 Ψ,
\end{aligned}
```
with (adiabatic) initial conditions ``δ = -\frac{3}{2} (1+w) Ψ`` and ``θ = \frac{1}{2} k^2 τ Ψ``,
following [Bertschinger and Ma (equation 30)](https://arxiv.org/pdf/astro-ph/9506072#%5B%7B%22num%22%3A70%2C%22gen%22%3A0%7D%2C%7B%22name%22%3A%22FitH%22%7D%2C387%5D).

Next, we must simply pack this into a symbolic component that represents the w₀wₐ dark energy species:
```@example ext
using ModelingToolkit # load to create custom components
using SymBoltz: τ, D, k # load conformal time, derivative and perturbation wavenumber
g = M1.g # reuse metric of original model

# 1. Create parameters (will resurface as tunable numbers in the full model)
pars = @parameters w₀ wₐ cₛ² Ω₀

# 2. Create variables
vars = @variables ρ(τ) P(τ) w(τ) δ(τ, k) θ(τ, k) σ(τ, k)

# 3. Specify equations (~ means equality in ModelingToolkit)
eqs = [
    # Background equations
    w ~ w₀ + wₐ * (1 - g.a) # equation of state
    ρ ~ 3/(8*Num(π))*Ω₀ * g.a^(-3*(1+w₀+wₐ)) # D(ρ) ~ -3 * g.ℰ * ρ * (1 + w) # energy density (ρ₀ => 3/(8*Num(π)) * exp(+3*wₐ) * Ω₀)
    P ~ w * ρ # pressure

    # Perturbation equations
    D(δ) ~ -(1 + w) * (θ - 3*g.Φ) - 3 * g.ℰ * (cₛ² - w) * δ # energy overdensity
    D(θ) ~ -g.ℰ * (1 - 3*w) * θ - D(w) / (1 + w) * θ + cₛ² / (1 + w) * k^2 * δ - k^2 * σ + k^2 * g.Ψ # momentum
    σ ~ 0 # shear stress
]

# 4. Specify initial conditions (for perturbations)
initialization_eqs = [
    δ ~ -3/2 * (1+w) * g.Ψ
    θ ~ 1/2 * (k^2/g.ℰ) * g.Ψ
]

# 5. Pack into an ODE system called "X"
description = "w₀wₐ (CPL) dynamical dark energy"
@named X = System(eqs, τ, vars, pars; initialization_eqs, description)
```

Note that the w₀wₐ component only knows about itself (and the metric),
but is completely unaware of the theory of gravity, other species and other components.
Its "job" is only to expose the variables like `ρ`, `P`, `δ` and `σ` that source the Einstein equations.
This connection is made when the component is used to create a full cosmological model, as we will do next.

## 3. Create the extended model

Finally, we create a new `ΛCDM` model, but replace `Λ` by `X`, and call it the `w0waCDM` model:
```@example ext
M2 = ΛCDM(Λ = X, name = :w0waCDM)
hierarchy(M2; describe = true)
```
Now `M2.Λ` no longer exists, but `M2.X` contains our new dark energy species:
```@example ext
equations(M2.X)
```

## Solve and compare the models

To test, let us set some parameters and solve both models with one perturbation mode.
For the ΛCDM model:
```@example ext
θ1 = SymBoltz.parameters_Planck18(M1)
prob1 = CosmologyProblem(M1, θ1)
ks = 1.0
sol1 = solve(prob1, ks)
```
And for the w₀wₐCDM model:
```@example ext
θ2 = merge(θ1, Dict(
    M2.X.w₀ => -0.9,
    M2.X.wₐ => 0.2,
    M2.X.cₛ² => 1.0
))
prob2 = CosmologyProblem(M2, θ2)
sol2 = solve(prob2, ks)
```
Let us compare ``H(τ)`` and ``Ψ(k,τ)`` at equal scale factors ``a(τ)``:
```@example ext
lgas = range(-3, 0, length=500) # log10(a)
H1s = sol1(log10(M1.g.a) => lgas, M1.g.H)
H2s = sol2(log10(M2.g.a) => lgas, M2.g.H)
Ψ1s = sol1(ks, log10(M1.g.a) => lgas, M1.g.Ψ)
Ψ2s = sol2(ks, log10(M2.g.a) => lgas, M2.g.Ψ)

using Plots
plot(lgas, H2s ./ H1s; xlabel = "lg(a)", label = "H₂ / H₁")
plot!(lgas, Ψ1s ./ Ψ2s; xlabel = "lg(a)", label = "Ψ₂ / Ψ₁")
```
