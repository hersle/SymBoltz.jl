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
nothing # hide
```

The symbolic model has a hierarchical structure with several components of the Einstein-Boltzmann system:
```@example 1
using Plots
plot(M.sys)
```
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
Internally, the full model is separated into sequential computational stages, like the background and perturbation systems.

The hierarchical structure can be inspected interactively as a [ModelingToolkit](https://docs.sciml.ai/ModelingToolkit) system by evaluating `M.sys`, `M.sys.<TAB>`, `M.sys.G` and so on (with TAB-completion).
For example, to see all equations for the theory of gravity:
```@example 1
equations(M.sys.G)
```

See the [tutorial on Creating extended models](@ref "Creating extended models") to learn more about the model representation and how to define extended cosmological models.

## 2. Solve the model

Next, we set our desired cosmological parameters and wavenumbers, and solve the model:
```@example 1
pars = [
    M.sys.γ.Ω0 => 5e-5
    M.sys.c.Ω0 => 0.27
    M.sys.b.Ω0 => 0.05
    M.sys.ν.Neff => 3.0
    M.sys.g.h => 0.7
    M.sys.b.rec.Yp => 0.25
]
ks = 10 .^ range(-3, 0, length=100) / SymBoltz.k0 # TODO: improve on units
sol = solve(M, pars, ks)
```

To solve only the background, you can simply omit the `ks` argument: `solve(prob, pars, ks)`.

## 3. Access the solution

You are now free to do whatever you want with the solution object.
It can be called like
- `sol(t, var)` to access the background variable(s) `var` as a function of the conformal time(s) `t`, or
- `sol(k, t, var)` to access the perturbed variable(s) `var` as a function of the wavenumber(s) `k` and conformal time(s) `t`.
This will interpolate linearly in $\log(k)$, and with the ODE solver's interpolator in $t$.

For example, to get the reduced Hubble function $E(t) = H(t) / H_0$ for 1000 log-spaced conformal times:
```@example 1
ts = exp.(range(log.(extrema(sol[SymBoltz.t]))..., length=1000)) # TODO: replace with M.t
Es = sol(ts, M.sys.g.E)
```
Similarly, to get $\Phi(k,t)$ for the 100 wavenumbers we solved for and the same 1000 log-spaced conformal times:
```@example 1
Φs = sol(ks, ts, M.sys.g.Φ)
```

This can be plotted with `using Plots; plot(log10.(ts), transpose(Φs))`, but this is even easier with the included plot recipe:
```@example 1
using Plots
plot(sol, ks[begin:5:end], log10(M.sys.g.a), M.sys.g.Φ) # lg(a) vs. Φ for every 5th wavenumber
```

We can also calculate the power spectrum for a desired species (here: cold dark matter with `M.c`):
```@example 1
Ps = power_spectrum(sol, M.sys.c, ks)
plot(log10.(ks), log10.(Ps); xlabel="lg(k)", ylabel="lg(P)", label=nothing)
```