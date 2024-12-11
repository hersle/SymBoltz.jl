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
    M.γ.T₀ => 2.7,
    M.c.Ω₀ => 0.27,
    M.b.Ω₀ => 0.05,
    M.ν.Neff => 3.0,
    M.g.h => 0.7,
    M.b.rec.Yp => 0.25,
    M.h.m => 0.06 * SymBoltz.eV / SymBoltz.c^2
)
ks = 10 .^ range(-4, 1, length=100) / u"Mpc"
sol = solve(M, pars, ks)
```

To solve only the background, you can simply omit the `ks` argument: `solve(prob, pars)`.

## 3. Access the solution

You are now free to [do whatever you want with the solution object](@ref "Solution handling").
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
ks_plot = [1e-3, 1e-2, 1e-1, 1e-0] / u"Mpc"
plot(sol, ks_plot, log10(M.g.a), M.g.Φ) # lg(a) vs. Φ for 4 wavenumbers
```

We can also calculate the power spectrum:
```@example getting_started
Ps = spectrum_matter(sol, ks)
plot(log10.(ks/u"1/Mpc"), log10.(Ps/u"Mpc^3"); xlabel = "lg(k/Mpc⁻¹)", ylabel = "lg(P/Mpc³)", label = nothing)
```
Similarly, we can calculate the angular CMB power spectrum:
```@example getting_started
ls = 10:10:1000
Cls = spectrum_cmb(:TT, M, pars, ls)
plot(ls, @. Cls * ls * (ls + 1) / 2π / 1e-6^2; xlabel = "l", ylabel = "l (l+1) Cₗ / 2π / (μK)²", label = nothing)
```

And here is a condensed plot with several quantities:
```@example getting_started
p = plot(layout=(3, 3), size=(900, 700), tickfontsize=6, labelfontsize=6, legendfontsize=5)
plot!(p[1], sol, log10(M.g.a), [M.b.ρ, M.c.ρ, M.γ.ρ, M.ν.ρ, M.h.ρ, M.Λ.ρ] ./ M.G.ρ)
plot!(p[2], sol, log10(M.g.a), [M.b.w, M.c.w, M.γ.w, M.ν.w, M.h.w, M.Λ.w])
plot!(p[3], sol, log10(M.g.a), log10(M.g.H / M.g.H₀))
plot!(p[4], sol, log10(M.g.a), [M.b.rec.XHe⁺⁺, M.b.rec.XHe⁺, M.b.rec.XH⁺, M.b.rec.Xe])
plot!(p[5], sol, log10(M.g.a), log10.([M.b.rec.Tγ, M.b.rec.Tb] ./ M.γ.T₀))
plot!(p[6], sol, log10(M.g.a), log10(abs(M.b.rec.τ)))
plot!(p[7], sol, ks_plot, log10(M.g.a), [M.g.Φ, M.g.Ψ])
plot!(p[8], sol, ks_plot, log10(M.g.a), log10.(abs.([M.b.δ, M.c.δ, M.γ.δ, M.ν.δ, M.h.δ])), klabel = false)
plot!(p[9], sol, ks_plot, log10(M.g.a), log10.(abs.([M.b.θ, M.c.θ, M.γ.θ, M.ν.θ])), klabel = false)
```