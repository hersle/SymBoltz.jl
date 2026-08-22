# Plotting and visualization

!!! note
    This page is work-in-progress.

First, set up a cosmological problem to solve:
```@example plot
using SymBoltz, Unitful, UnitfulAstro
M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
```

## Static plot recipes

Use SymBoltz' included plot recipes to plot the evolution of background and perturbation quantities over time:
```@example plot
import Plots
ks = [1e-3, 1e-2, 1e-1, 1e0] / u"Mpc"
sol = solve(prob, ks)
p1 = Plots.plot(sol, M.χ, M.g.a)
p2 = Plots.plot(sol, log10(M.g.a), [M.g.Φ, M.g.Ψ], ks)
Plots.plot(p1, p2; layout = (2, 1), size = (600, 600))
```

Plot the evolution of $Φ(k,τ)$ over conformal time $τ$:
```@example plot
ks = 10 .^ range(-3, 0, length=100) / u"Mpc"
sol = solve(prob, ks)
τs = [0.0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0]
pal = Plots.palette([:black, :orange], length(τs))
color = permutedims([pal[i] for i in eachindex(τs)])
labels = permutedims(map(τ -> "τ/H₀⁻¹ = $τ", τs))
Plots.plot(log10.(ks*u"Mpc"), permutedims(sol(M.g.Φ, τs, ks)); xlabel = "k/Mpc⁻¹", ylabel = "Φ(τ,k)", color, labels)
```

Visualize the CMB source function $S₀(k,τ)$ in a 3D plot:
```@example plot
ks = range(0.0, 0.3, length=100)[2:end] / u"Mpc"
sol = solve(prob, ks)
Plots.surface(sol, M.τ, M.k, M.ST; xlims = (0.03, 0.09))
```

## Interactive visualization

The excellent [Makie plotting library](https://docs.makie.org/stable/) can be used to interactively visualize results:
```julia
using GLMakie
```
The function `plot_interactive` builds a figure with one slider per parameter.
Whenever a slider is moved, it updates the problem with the new parameter values and replots the data returned by a user-supplied function of the new problem:
```@example plot
using CairoMakie # for showing static doc image # hide

obspars = [
    M.g.h => 0.60:0.01:0.70,
    M.c.Ω₀ => 0.20:0.01:0.30,
    M.b.Ω₀ => 0.02:0.01:0.10,
    M.γ.T₀ => 2.50:0.01:3.00,
    M.h.m_eV => 0.01:0.01:0.15,
    M.b.YHe => 0.20:0.01:0.30,
    M.ν.Neff => 2.90:0.01:3.10
]
fig = plot_interactive(prob, obspars; xlabel = "lg(a)", ylabel = "Xₑ") do prob′
    sol = solve(prob′)
    τ = SymBoltz.timeseries(sol)
    collect(zip(log10.(sol(M.g.a, τ)), sol(M.b.Xe, τ))) # [(x1, y1), (x2, y2), ...]
end
```

!!! info
    This plot is static due to limitations in the documentation building system.
    It is interactive when you execute the code locally.

We can make a similar plot for the matter power spectrum $P(k; θ)$ as a function of cosmological parameters:

```@example plot
using DataInterpolations # for smoothing
obspars = [
    obspars; # extend vector from above
    M.I.ln_As1e10 => 2.0:0.1:4.0
    M.I.ns => 0.90:0.01:1.10
]
fig = plot_interactive(prob, obspars; xlabel = "lg(k / Mpc⁻¹)", ylabel = "lg(P / Mpc³)") do prob′
    lgks = unique([-4:0.5:-3; -3:0.2:-2; -2:0.05:0]) # as few points as possible
    ks = 10 .^ lgks / u"Mpc"
    Ps = spectrum_matter(prob′, ks; ptopts = (alg = SymBoltz.TRBDF2(), reltol = 1e-4, abstol = 1e-4))
    lgPs = log10.(Ps/u"Mpc^3")

    # smoothen with spline and sample more densely
    lgPspline = CubicSpline(lgPs, lgks)
    lgks = range(lgks[begin], lgks[end]; step = 0.01)
    lgPs = lgPspline(lgks)

    collect(zip(lgks, lgPs)) # [(x1, y1), (x2, y2), ...]
end
```

!!! info
    This plot is static due to limitations in the documentation building system.
    It is interactive when you execute the code locally.
