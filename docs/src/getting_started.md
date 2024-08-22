# Getting started

[Install SymBoltz.jl](@ref "Installation") and load it with
```@example 1
using SymBoltz
``` 

## 1. Create the model

The first step is to define our cosmological model.
This is a symbolic representation of the variables and equations that describe the various components in the universe, such as the theory of gravity (like general relativity) and particle species (like the cosmological constant, cold dark matter, photons and baryons).
It will be used both to create the numerical problem to solve, and to access variables in its solution.

To get started, we will simply load the standard ΛCDM model:
```@example 1
M = ΛCDM()
```

As shown, the symbolic model is structured with a hierarchy of components, each of which contains chunks of the Einstein-Boltzmann system (it can also be displayed graphically with `using Plots; plot(M)`).
.
In this case, the components are
the spacetime metric (`g`, for $g_{\mu\nu}$),
the theory of gravity (`G`, here general relativity),
photons (`γ`),
massless neutrinos (`ν`),
cold dark matter (`c`),
baryons (`b`),
massive neutrinos (`h`),
dark energy as a cosmological constant (`Λ`)
and recombination (`rec`).

Internally, these components are combined into a full Einstein-Boltzmann system.
The full system is also separated into sequential computational stages, like the background and perturbation systems.

The hierarchical structure can be inspected interactively as a [ModelingToolkit](https://docs.sciml.ai/ModelingToolkit) system by evaluating `M`, `M.<TAB>`, `M.G` and so on (with TAB-completion).
For example, to see all equations for the theory of gravity:
```@example 1
equations(M.G)
```

## 2. Solve the model

Next, we set our desired cosmological parameters and wavenumbers, and solve the model:
```@example 1
pars = [
    M.γ.Ω0 => 5e-5
    M.c.Ω0 => 0.27
    M.b.Ω0 => 0.05
    M.ν.Neff => 3.0
    M.g.h => 0.7
    M.b.rec.Yp => 0.25
]
ks = 10 .^ range(-3, 0, length=100) / SymBoltz.k0 # TODO: improve on units
sol = solve(M, pars, ks)
```

To solve only the background, you can simply omit the `ks` argument: `solve(prob, pars)`.

## 3. Access the solution

You are now free to do whatever you want with the solution object.
It can be called like
- `sol(t, y)` to get the variable(s) $y(t)$ as a function of conformal time(s) $t$, or
- `sol(k, t, y)` to get the variable(s) $y(k,t)$ as a function of the wavenumber(s) $k$ and conformal time(s) $t$.
`y` can be *any* symbolic variable(s) in the model `M`.
The solution will automatically interpolate linearly in $\log(k)$, and with the ODE solver's custom interpolator in $t$.

For example, to get the reduced Hubble function $E(t) = H(t) / H_0$ for 1000 log-spaced conformal times:
```@example 1
ts = exp.(range(log.(extrema(sol[SymBoltz.t]))..., length=1000)) # TODO: replace with M.t
Es = sol(ts, M.g.E)
```
Similarly, to get $\Phi(k,t)$ for the 100 wavenumbers we solved for and the same 1000 log-spaced conformal times:
```@example 1
Φs = sol(ks, ts, M.g.Φ)
```

This can be plotted with `using Plots; plot(log10.(ts), transpose(Φs))`, but this is even easier with the included plot recipe:
```@example 1
using Plots
plot(sol, ks[begin:5:end], log10(M.g.a), M.g.Φ) # lg(a) vs. Φ for every 5th wavenumber
```

We can also calculate the power spectrum for a desired species (here: cold dark matter with `M.c`):
```@example 1
Ps = power_spectrum(sol, M.c, ks)
plot(log10.(ks), log10.(Ps); xlabel="lg(k)", ylabel="lg(P)", label=nothing)
```