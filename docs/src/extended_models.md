# Creating extended models

This tutorial shows how to create extended (beyond-ΛCDM) models.
As a simple example, we replace the cosmological constant with equation of state ``w(t) = -1`` by w₀wₐ dark energy with equation of state
```math
w(t) = \frac{P(t)}{ρ(t)} = w_0 + w_a (1 - a(t)),
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
M1 = ΛCDM(Λanalytical = false)
```
The cosmological constant species is stored in the component `M1.Λ`.
Here are the equations we will remove:
```@example ext
equations(M1.Λ)
```

## 2. Create the extended component

The background continuity equation
``\dot{ρ}(t) = -3 H(t) \, ρ(t) \, (1 + w(t))``
can be solved analytically with the w₀wₐ equation of state and the ansatz
```math
ρ(t) = ρ(t_0) \, a(t)^m \exp(n (1 - a(t))),
```
giving ``m = -3 (1 + w_0 + w_a)`` and ``n = -3 w_a``.
The perturbation equations are
```math
\begin{aligned}
\dot{δ} &= -(1 + w) (θ - 3 \dot{Φ}) - 3 \frac{\dot{a}}{a} (c_s^2 - w) δ, \\
\dot{θ} &= -\frac{\dot{a}}{a} (1 - 3w) θ - \frac{\dot{w}}{1+w} θ + \frac{c_s^2}{1+w} k^2 δ - k^2 σ + k^2 Ψ,
\end{aligned}
```
with (adiabatic) initial conditions ``δ = -\frac{3}{2} (1+w) Ψ`` and ``θ = \frac{1}{2} k^2 t Ψ``,
following [Bertschinger and Ma (equation 30)](https://arxiv.org/pdf/astro-ph/9506072#%5B%7B%22num%22%3A70%2C%22gen%22%3A0%7D%2C%7B%22name%22%3A%22FitH%22%7D%2C387%5D).

Next, we must simply pack this into a symbolic component that represents the w₀wₐ dark energy species:
```@example ext
using ModelingToolkit # load to create custom components
using SymBoltz: t, D, k, ϵ # load conformal time, derivative, perturbation wavenumber and perturbation book-keeper variables
g = M1.g # reuse metric component from the original model

# 1. Create parameters (will resurface as tunable numbers in the full model)
pars = @parameters w₀ wₐ cₛ² # ρ₀ Ω₀ # TODO: compare to analytical solution in the end

# 2. Create variables
vars = @variables ρ(t) P(t) w(t) δ(t) θ(t) σ(t)

# 3. Specify equations (~ means equality in ModelingToolkit)
eqs = [
    # Background equations (of order O(ϵ⁰)
    w ~ w₀ + wₐ * (1 - g.a) # equation of state
    D(ρ) ~ -3 * g.ℰ * ρ * (1 + w) # energy density
    P ~ w * ρ # pressure

    # Perturbation equations (mulitiplied by ϵ to mark them as order O(ϵ¹))
    D(δ) * ϵ ~ (-(1 + w) * (θ - 3*g.Φ) - 3 * g.ℰ * (cₛ² - w) * δ) * ϵ # energy overdensity
    D(θ) * ϵ ~ (-g.ℰ * (1 - 3*w) * θ - D(w) / (1 + w) * θ + cₛ² / (1 + w) * k^2 * δ - k^2 * σ + k^2 * g.Ψ) * ϵ # momentum
    σ * ϵ ~ 0 # shear stress
]

# 4. Specify initial conditions (for perturbations)
initialization_eqs = [
    δ * ϵ ~ -3/2 * (1+w) * g.Ψ * ϵ
    θ * ϵ ~ 1/2 * (k^2*t) * g.Ψ * ϵ
]

# 5. Pack into an ODE system called "X"
@named X = ODESystem(eqs, t, vars, pars; initialization_eqs)
```

Note that the w₀wₐ component only knows about itself (and the metric),
but is completely unaware of the theory of gravity, other species and other components.
Its "job" is only to expose the variables like `ρ`, `P`, `δ` and `σ` that source the Einstein equations.
This connection is made when the component is used to create a full cosmological model, as we will do next.

## 3. Create the extended model

Finally, we create a new `ΛCDM` model, but replace `Λ` by `X`, and call it the `w0waCDM` model:
```@example ext
M2 = ΛCDM(Λ = X, name = :w0waCDM, Λanalytical = false)
```
Now `M2.Λ` no longer exists, but `M2.X` contains our new equations:
```@example ext
equations(M2.X)
```

## Solve and compare the models

To test, let us set some parameters and solve both models with one perturbation mode.
For the ΛCDM model:
```@example ext
θ1 = Dict(
    M1.γ.T₀ => 2.7,
    M1.c.Ω₀ => 0.27,
    M1.b.Ω₀ => 0.05,
    M1.ν.Neff => 3.0,
    M1.g.h => 0.7,
    M1.b.rec.Yp => 0.25,
    M1.h.m => 0.06 * SymBoltz.eV / SymBoltz.c^2
)
ks = 1.0
sol1 = solve(M1, θ1, ks)
```
And for the w₀wₐCDM model:
```@example ext
θ2 = merge(θ1, Dict(
    M2.X.w₀ => -0.9,
    M2.X.wₐ => 0.2,
    M2.X.cₛ² => 1.0
))
sol2 = solve(M2, θ2, ks)
```
Let us compare ``H(t)`` and ``Ψ(k,t)`` at equal scale factors ``a(t)``:
```@example ext
lgas = range(-3, 0, length=500) # log10(a)
H1s = sol1(log10(M1.g.a), lgas, M1.g.H)
H2s = sol2(log10(M2.g.a), lgas, M2.g.H)
Ψ1s = sol1(log10(M1.g.a), ks, lgas, M1.g.Ψ)
Ψ2s = sol2(log10(M2.g.a), ks, lgas, M2.g.Ψ)

using Plots
# TODO: make convenience plot method that compares two models
plot(lgas, H2s ./ H1s; xlabel = "lg(a)", label = "H₂ / H₁")
plot!(lgas, Ψ1s ./ Ψ2s; xlabel = "lg(a)", label = "Ψ₂ / Ψ₁")
```
