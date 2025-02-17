# Observable quantities

This page shows observable quantities that can be derived from solutions of the Einstein-Boltzmann system, such as power spectra and distances.

## Primordial power spectra

```@docs
SymBoltz.spectrum_primordial
```

#### Example

```@example
using SymBoltz, Unitful, UnitfulAstro, Plots
M = SymBoltz.ΛCDM(K = nothing)
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
M = SymBoltz.ΛCDM(K = nothing)
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
M = SymBoltz.ΛCDM(K = nothing)
pars = SymBoltz.parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
ls = 10:5:1500

DlTTs, DlEEs, DlTEs = spectrum_cmb([:TT, :EE, :TE], prob, ls; normalization = :Dl, unit = u"μK")
# TODO: compare EE and TE with CLASS # hide
pTT = plot(ls, DlTTs; ylabel = "Dₗᵀᵀ")
pEE = plot(ls, DlEEs; ylabel = "Dₗᴱᴱ")
pTE = plot(ls, DlTEs; ylabel = "Dₗᵀᴱ", xlabel = "l")
plot(pTT, pEE, pTE, layout = (3, 1), size = (600, 700), legend = nothing)
```

## Two-point correlation function

```@docs
SymBoltz.correlation_function
```

#### Example

```@example
using SymBoltz, Unitful, UnitfulAstro, Plots
M = SymBoltz.ΛCDM(K = nothing)
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
M = SymBoltz.ΛCDM(K = nothing)
pars = SymBoltz.parameters_Planck18(M)
prob = CosmologyProblem(M, pars, Dict(M.Λ.Ω₀ => 0.5), [M.g.ℰ ~ 1])
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
M = SymBoltz.RMΛ()
pars = Dict(
    M.r.Ω₀ => 5e-5,
    M.m.Ω₀ => 0.3,
    M.K.Ω₀ => 0.1,
    M.r.T₀ => NaN,
    M.g.h => 0.7
)
pars[M.Λ.Ω₀] = 1 - pars[M.r.Ω₀] - pars[M.m.Ω₀] - pars[M.K.Ω₀]
prob = CosmologyProblem(M, pars)
sol = solve(prob)

zs = 0.0:1.0:10.0
ts = SymBoltz.timeseries(sol, M.g.z, zs) # times at given redshifts
dLs = SymBoltz.distance_luminosity(sol, ts) / SymBoltz.Gpc
plot(zs, dLs; mark=:dot, xlabel="z", ylabel="dL / Gpc", label=nothing)
```
