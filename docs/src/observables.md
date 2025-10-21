# Observable quantities

This page shows observable quantities that can be derived from solutions of the Einstein-Boltzmann system, such as power spectra and distances.

## Primordial power spectra

```@docs
SymBoltz.spectrum_primordial
```

#### Example

```@example
using SymBoltz, Unitful, UnitfulAstro, Plots
M = SymBoltz.ΛCDM()
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

```@example
using SymBoltz, Unitful, UnitfulAstro, Plots
M = SymBoltz.ΛCDM()
pars = SymBoltz.parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
ks = 10 .^ range(-5, +2, length=200) / u"Mpc"
sol = solve(prob, ks)

# Linear power spectrum
Ps = spectrum_matter(sol, ks)
plot(log10.(ks*u"Mpc"), log10.(Ps/u"Mpc^3"); xlabel = "log10(k/Mpc⁻¹)", ylabel = "log10(P/Mpc³)", label = "linear (SymBoltz)")

# Nonlinear power spectrum (from halofit)
Ps = spectrum_matter_nonlinear(sol, ks)
plot!(log10.(ks*u"Mpc"), log10.(Ps/u"Mpc^3"); label = "non-linear (halofit)", legend_position = :bottomleft)
```

## CMB power spectra

```@docs
SymBoltz.spectrum_cmb
```

#### Example

```@example
# TODO: more generic, source functions, ... # hide
using SymBoltz, Unitful, UnitfulAstro, Plots
M = SymBoltz.ΛCDM()
pars = SymBoltz.parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
ls = 10:5:1500

jl = SphericalBesselCache(ls)
Dls = spectrum_cmb([:TT, :EE, :TE], prob, jl; normalization = :Dl, unit = u"μK")
pTT = plot(ls, Dls[:, 1]; ylabel = "Dₗᵀᵀ")
pEE = plot(ls, Dls[:, 2]; ylabel = "Dₗᴱᴱ")
pTE = plot(ls, Dls[:, 3]; ylabel = "Dₗᵀᴱ", xlabel = "l")
plot(pTT, pEE, pTE, layout = (3, 1), size = (600, 700), legend = nothing)
```

## Two-point correlation function

```@docs
SymBoltz.correlation_function
```

#### Example

```@example
using SymBoltz, Unitful, UnitfulAstro, Plots
M = SymBoltz.ΛCDM()
pars = SymBoltz.parameters_Planck18(M)
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
M = SymBoltz.ΛCDM()
pars = SymBoltz.parameters_Planck18(M)
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
M = SymBoltz.RMΛ(K = SymBoltz.curvature(SymBoltz.metric()))
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
dLs = SymBoltz.distance_luminosity(sol, τs) / SymBoltz.Gpc
@assert isapprox(dLs[begin], 0.0; atol = 1e-14) || zs[begin] != 0.0 # ensure bug does not reappear # hide
plot(zs, dLs; marker=:dot, xlabel="z", ylabel="dL / Gpc", label=nothing)
```

## Sound horizon (BAO scale)

```@docs
sound_horizon
```

```@example
using SymBoltz, Plots
M = SymBoltz.ΛCDM()
pars = SymBoltz.parameters_Planck18(M)
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
