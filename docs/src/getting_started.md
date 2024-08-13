# Getting started

[Install SymBoltz.jl](@ref "Installation") and load it with
```@example 1
using ModelingToolkit
import SymBoltz
``` 

## 1. Define the model

The first step is to define our cosmological model.
This is a symbolic representation of the variables and equations that describe the various components in the universe, such as the theory of gravity (like general relativity) and particle species (like the cosmological constant, cold dark matter, photons and baryons).
It will be used both to create the numerical problem to solve, and to access variables in its solution.

To get started, we will simply load the standard ΛCDM model:
```@example 1
@named M = SymBoltz.ΛCDM()
nothing # hide
```

The symbolic model has a hierarchical structure, with multiple subsystems representing different components of the Einstein-Boltzmann system:
```@example 1
using Plots
plot(M)
```
In this case, the components are:
| Symbol | Component                        |
| ------ | -------------------------------- |
| `g`    | spacetime metric ($g_{\mu\nu}$)  |
| `G`    | gravity                          |
| `γ`    | photons                          |
| `ν`    | massless neutrinos               |
| `c`    | cold dark matter                 |
| `b`    | baryons                          |
| `h`    | massive neutrinos                |
| `Λ`    | cosmological constant            |
| `rec`  | recombination                    |

The hierarchical structure of the whole symbolic model can be inspected interactively as a [ModelingToolkit](https://docs.sciml.ai/ModelingToolkit) system by evaluating `M`, `M.<TAB>`, `M.G` and so on (with TAB-completion).
For example, to see all equations for the theory of gravity:
```@example 1
equations(M.G)
```

See TODO to learn more about the model representation and how to define extended cosmological models.

## 2. Build the problem

The next step is to create a cosmological problem from the symbolic model.
This separates the full model into sequential computational stages, like the background and perturbation systems.

```@example 1
prob = SymBoltz.CosmologyProblem(M)
nothing # hide
```

## 3. Solve the problem

Finally, we set our desired cosmological parameters and wavenumbers, and solve the problem:
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
sol = SymBoltz.solve(prob, pars, ks)
```

To solve only the background, you can simply omit the `ks` argument: `SymBoltz.solve(prob, pars, ks)`.

## 4. Access the solution

You are now free to do whatever you want with the solution object.
It can be called like
- `sol(t, var)` to access the background variable(s) `var` as a function of the conformal time(s) `t`, or
- `sol(k, t, var)` to access the perturbed variable(s) `var` as a function of the wavenumber(s) `k` and conformal time(s) `t`.
This will interpolate linearly in $\log(k)$, and with the ODE solver's interpolator in $t$.

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
Ps = SymBoltz.P(sol, M.c, ks)
plot(log10.(ks), log10.(Ps); xlabel="lg(k)", ylabel="lg(P)", label=nothing)
```