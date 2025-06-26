# Plotting and visualization

!!! note
    This page is work-in-progress.

First, set up a cosmological problem to solve:
```@example plot
using SymBoltz, Unitful, UnitfulAstro
M = SymBoltz.ΛCDM()
pars = SymBoltz.parameters_Planck18(M)
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
using CairoMakie
τs = range(0.05, 0.08; length=50)
ks = range(0.0, 0.3, length=100) / u"Mpc"
sol = solve(prob, ks)

xs = τs
ys = ks*u"Mpc"
zs = sol(M.ST0, τs, ks)

fig = Figure()
ax = Axis3(fig[1,1], azimuth = π/4, xlabel = "k/Mpc⁻¹", ylabel = "τ/H₀⁻¹", zlabel = "S₀(τ,k)")
cmax = min(-minimum(filter(!isnan, zs)), maximum(filter(!isnan, zs))) # saturate both ends of color scale
surface!(ax, xs, ys, zs; alpha = 0.9, colormap = :seismic, colorrange = (-cmax, +cmax))

fig
```

## Interactive visualization

The excellent [Makie plotting library](https://docs.makie.org/stable/) can be used to interactively visualize results:
```julia
using GLMakie
```
We can now build a simple interactive interface that updates a plotted function when the user drags parameter value sliders:
```@example plot
using CairoMakie # for showing static doc image # hide

# TODO: move into SymBoltz plotting recipes once MakieCore supports @lift etc.
# TODO: take in (xyfunc) as a function that can plot anything # hide
# TODO: handle wavenumber perturbations etc. # hide
# TODO: plot multiple y vars # hide
function plot_interactive(prob::CosmologyProblem, xyfunc::Function, obspars...; xlabel = "", ylabel = "", kwargs...)
    fig = Figure(size = (800, 800), fontsize = 18)
    ax = Axis(fig[1, 1]; xlabel, ylabel)
    sg = SliderGrid(fig[2, 1], ((label = SymBoltz.displayname(par), range = range, startvalue = prob.bg.ps[par]) for (par, range) in obspars)...)

    obs = Observable([pars[par] for (par, _) in obspars]) # observable array for all parameters
    obss = [Observable(pars[par]) for (par, _) in obspars] # array of observables for each parameter

    for i in 1:length(obspars)
        # update full parameter array when individual parameters change
        on(obss[i]) do val
            vals = obs[]
            vals[i] = val
            obs[] = vals
        end
        connect!(obss[i], sg.sliders[i].value)
    end

    lines!(ax, @lift(xyfunc($obs)); kwargs...)

    return fig
end
function plot_interactive(prob::CosmologyProblem, xvar::SymBoltz.Num, yvar::SymBoltz.Num, obspars...; N = 2000, kwargs...)
    # select points according to cumulative distribution of ODE solver's adaptive time steps
    function τs(sol)
        τ = sol.bg.t
        τcum(n) = τ[Int(floor(n))] + (τ[min(Int(floor(n))+1,length(τ))] - τ[Int(floor(n))]) * (n - floor(n)) # <3
        return τcum.(range(1, length(τ), length = N))
    end
    probgen = SymBoltz.parameter_updater(prob, [par for (par, _) in obspars])
    function xyfunc(θ)
        prob = probgen(θ)
        sol = solve(prob)
        τ = τs(sol)
        xs = sol(xvar, τ)
        ys = sol(yvar, τ)
        return collect(zip(xs, ys)) # [(x1, y1), (x2, y2), ...]
    end
    return plot_interactive(prob, xyfunc, obspars...; xlabel = SymBoltz.displayname(xvar), ylabel = SymBoltz.displayname(yvar), kwargs...)
end

obspars = [
    M.g.h => 0.60:0.01:0.70,
    M.c.Ω₀ => 0.20:0.01:0.30,
    M.b.Ω₀ => 0.02:0.01:0.10,
    M.γ.T₀ => 2.50:0.01:3.00,
    M.h.m_eV => 0.01:0.01:0.15,
    M.b.rec.Yp => 0.20:0.01:0.30,
    M.ν.Neff => 2.90:0.01:3.10
]
fig = plot_interactive(prob, log10(M.g.a), M.b.rec.Xe, obspars...)
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
probgen = SymBoltz.parameter_updater(prob, [par for (par, _) in obspars])
function xyfunc(θ)
    prob = probgen(θ)
    lgks = unique([-4:0.5:-3; -3:0.2:-2; -2:0.05:0]) # as few points as possible
    ks = 10 .^ lgks / u"Mpc"
    Ps = spectrum_matter(prob, ks; ptopts = (alg = SymBoltz.TRBDF2(), reltol = 1e-4, abstol = 1e-4))
    lgPs = log10.(Ps/u"Mpc^3")

    # smoothen with spline and sample more densely
    lgPspline = CubicSpline(lgPs, lgks)
    lgks = range(lgks[begin], lgks[end]; step = 0.01)
    lgPs = lgPspline(lgks)

    return collect(zip(lgks, lgPs)) # [(x1, y1), (x2, y2), ...]
end
fig = plot_interactive(prob, xyfunc, obspars...; xlabel = "lg(k / Mpc⁻¹)", ylabel = "lg(P / Mpc³)")
```

!!! info
    This plot is static due to limitations in the documentation building system.
    It is interactive when you execute the code locally.
