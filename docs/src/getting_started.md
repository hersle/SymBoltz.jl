# Getting started

[Install SymBoltz.jl](@ref "Installation") and load it with
```@example getting_started
using SymBoltz
``` 

## 1. Create the model

The first step is to define our cosmological model.
This is a symbolic representation of the variables and equations that describe the various components in the universe, such as the theory of gravity (like general relativity) and particle species (like the cosmological constant, cold dark matter, photons and baryons).
It will be used both to create the numerical problem to solve, and to access variables in its solution.

To get started, we will simply load the standard ΛCDM model:
```@example getting_started
M = ΛCDM()
```

As shown, the symbolic model is structured with a hierarchy of components, each of which contains chunks of the Einstein-Boltzmann system (it can also be displayed graphically with `using Plots; plot(M)`).
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

The hierarchical structure can be inspected interactively (with TAB-completion) as a [ModelingToolkit](https://docs.sciml.ai/ModelingToolkit) system by evaluating `M`, `M.G`, `M.g.ρ` and so on.
For example, to see all equations for the theory of gravity:
```@example getting_started
equations(M.G)
```

## 2. Solve the model

Next, we set our desired cosmological parameters and wavenumbers, and solve the model:
```@example getting_started
using Unitful, UnitfulAstro # for interfacing without internal code units
pars = Dict(
    M.γ.T0 => 2.7,
    M.c.Ω0 => 0.27,
    M.b.Ω0 => 0.05,
    M.ν.Neff => 3.0,
    M.g.h => 0.7,
    M.b.rec.Yp => 0.25
)
ks = 10 .^ range(-3, 0, length=100) / u"Mpc"
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

For example, to get the reduced Hubble function $E(t) = H(t) / H_0$ for 300 log-spaced conformal times:
```@example getting_started
ts = exp.(range(log.(extrema(sol[M.t]))..., length=300))
Es = sol(ts, M.g.E)
```
Similarly, to get $\Phi(k,t)$ for the 100 wavenumbers we solved for and the same 300 log-spaced conformal times:
```@example getting_started
Φs = sol(ks, ts, M.g.Φ)
```

This can be plotted with `using Plots; plot(log10.(ts), transpose(Φs))`, but this is even easier with the included plot recipe:
```@example getting_started
using Plots
plot(sol, [1e-3, 1e-2, 1e-1, 1e-0] / u"Mpc", log10(M.g.a), M.g.Φ) # lg(a) vs. Φ for 4 wavenumbers
```

We can also calculate the power spectrum for a desired species (here: cold dark matter with `M.c`):
```@example getting_started
Ps = power_spectrum(sol, M.c, ks)
plot(log10.(ks/u"1/Mpc"), log10.(Ps/u"Mpc^3"); xlabel = "lg(k/Mpc⁻¹)", ylabel = "lg(P/Mpc³)", label = nothing)
```
Similarly, we can calculate the angular CMB power spectrum:
```@example getting_started
ls = 10:10:1000
Cls = Cl(M, pars, ls)
plot(ls, @. Cls * ls * (ls + 1) / 2π / 1e-6^2; xlabel = "l", ylabel = "l (l+1) Cₗ / 2π / (μK)²", label = nothing) # TODO: fix
```
