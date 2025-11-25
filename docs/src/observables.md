# Observable quantities

This page shows observable quantities that can be derived from solutions of the Einstein-Boltzmann system, such as power spectra and distances.

## Primordial power spectra

```@docs
SymBoltz.spectrum_primordial
```

#### Example

```@example
using SymBoltz, Unitful, UnitfulAstro, Plots
M = ΛCDM()
pars = Dict(M.g.h => 0.7, M.I.As => 2e-9, M.I.ns => 0.95)
ks = 10 .^ range(-5, +1, length=100) / u"Mpc"
Ps = spectrum_primordial(ks, M, pars)
plot(log10.(ks*u"Mpc"), log10.(Ps/u"Mpc^3"); xlabel = "log10(k/Mpc⁻¹)", ylabel = "log10(P/Mpc³)")
```

## Matter power spectra

```@docs
SymBoltz.spectrum_matter
SymBoltz.spectrum_matter_nonlinear
```

#### Example

With explicitly chosen wavenumbers:

```@example matter
using SymBoltz, Unitful, UnitfulAstro, Plots
M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
ks = 10 .^ range(-5, +2, length=100) / u"Mpc"
sol = solve(prob, ks)

# Linear power spectrum
Ps = spectrum_matter(sol, ks)
Plots.plot(log10.(ks*u"Mpc"), log10.(Ps/u"Mpc^3"); xlabel = "log10(k/Mpc⁻¹)", ylabel = "log10(P/Mpc³)", label = "linear (SymBoltz)", marker = :circle, markersize = 2)

# Nonlinear power spectrum (from halofit)
Ps = spectrum_matter_nonlinear(sol, ks)
Plots.plot!(log10.(ks*u"Mpc"), log10.(Ps/u"Mpc^3"); label = "non-linear (halofit)", marker = :circle, markersize = 2, legend_position = :bottomleft)
```

With adaptively chosen wavenumbers on an interval:

```@example matter
ks, Ps = spectrum_matter(prob, (1e0, 1e3))
plot(log10.(ks), log10.(Ps); xlabel = "k / (H₀/c)", ylabel = "P / (c/H₀)³", label = "$(length(ks)) × k (adaptive)", marker = :circle, markersize = 2)
```

## CMB power spectra

```@docs
SymBoltz.spectrum_cmb
```

#### Example

```@example
# TODO: more generic, source functions, ... # hide
using SymBoltz, Plots
M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)

ls = 25:25:3000 # 25, 50, ..., 3000
jl = SphericalBesselCache(ls)
modes = [:TT, :EE, :TE, :ψψ, :ψT, :ψE]
Dls = spectrum_cmb(modes, prob, jl; normalization = :Dl)

Plots.plot(ls, log10.(abs.(Dls)); xlabel = "l", ylabel = "lg(Dₗ)", label = permutedims(String.(modes)))
```

## Two-point correlation function

```@docs
SymBoltz.correlation_function
```

#### Example

```@example
using SymBoltz, Unitful, UnitfulAstro, Plots
M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
ks = 10 .^ range(-5, +3, length=300) / u"Mpc"
sol = solve(prob, ks)
rs, ξs = correlation_function(sol)
rs = rs / (SymBoltz.k0*sol.bg.ps[:h]) * u"Mpc" # TODO: auto units
plot(rs, @. ξs * rs^2; xlims = (0, 200), xlabel = "r", ylabel = "r² ξ")
```

## Matter density fluctuations

```@docs
SymBoltz.variance_matter
SymBoltz.stddev_matter
```

```@example
using SymBoltz, Unitful, UnitfulAstro, Plots
M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
ks = 10 .^ range(-5, +3, length=300) / u"Mpc"
sol = solve(prob, ks)

h = sol[M.g.h]
Rs = 10 .^ range(0, 2, length=100) * u"Mpc"
σs = stddev_matter.(sol, Rs)
plot(log10.(Rs/(u"Mpc"/h)), log10.(σs); xlabel = "lg(R / (Mpc/h))", ylabel = "lg(σ)", label = nothing)

R8 = 8 * u"Mpc"/h
σ8 = stddev_matter(sol, R8)
scatter!((log10(R8/(u"Mpc"/h)), log10(σ8)), series_annotation = text("  σ₈ = $(round(σ8; digits=3))", :left), label = nothing)
```

## Luminosity distance

```@docs
SymBoltz.distance_luminosity
```

```@example
using SymBoltz, Plots
M = RMΛ(K = SymBoltz.curvature(SymBoltz.metric()))
pars = Dict(
    M.r.Ω₀ => 5e-5,
    M.m.Ω₀ => 0.3,
    M.K.Ω₀ => 0.1,
    M.r.T₀ => NaN,
    M.g.h => 0.7
)
prob = CosmologyProblem(M, pars)
sol = solve(prob)

zs = 0.0:1.0:10.0
τs = SymBoltz.timeseries(sol, M.g.z, zs) # times at given redshifts
dLs = distance_luminosity(sol, τs) / SymBoltz.Gpc
@assert isapprox(dLs[begin], 0.0; atol = 1e-14) || zs[begin] != 0.0 # ensure bug does not reappear # hide
plot(zs, dLs; marker=:dot, xlabel="z", ylabel="dL / Gpc", label=nothing)
```

## Sound horizon (BAO scale)

```@docs
sound_horizon
```

```@example
using SymBoltz, Plots
M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
sol = solve(prob)
τs = sol[M.τ]
rs = sound_horizon(sol)
plot(τs, rs; xlabel = "τ / H₀⁻¹", ylabel = "rₛ / (c/H₀)")
```

## Source functions

```@docs
source_grid
source_grid_adaptive
```

```@example
using SymBoltz, Plots, DataInterpolations
M = ΛCDM(h = nothing, ν = nothing)
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
sol = solve(prob)

τs = sol[M.τ] # conformal times in background solution
ks = [1.0, 2000.0] # initial coarse grid
ks, Ss = source_grid_adaptive(prob, [M.ST0], τs, ks; atol = 5.0)
iτ = argmax(sol[M.b.v]) # index of decoupling time
iτs = iτ-75:iτ+75 # indices around decoupling
p1 = surface(ks, τs[iτs], Ss[1, iτs, :]; camera = (45, 25), xlabel = "k", ylabel = "τ", zlabel = "S", colorbar = false)

lgas = -6.0:0.2:0.0
τs = LinearInterpolation(sol[M.τ], sol[log10(M.g.a)])(lgas) # τ at given lg(a)
ks = 5.0:5.0:100.0
Ss = source_grid(prob, [M.g.Ψ], τs, ks)
p2 = wireframe(ks, lgas, Ss[1, :, :]; camera = (75, 20), xlabel = "k", ylabel = "lg(a)", zlabel = "Φ")

plot(p1, p2)
```

## Line-of-sight integration

```@docs
SymBoltz.los_integrate
```
