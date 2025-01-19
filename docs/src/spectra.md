# Sources and spectra

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
prob = CosmologyProblem(M, pars, Dict(M.Λ.Ω₀ => 0.5), [M.g.ℰ ~ 1])
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
prob = CosmologyProblem(M, pars, Dict(M.Λ.Ω₀ => 0.5), [M.g.ℰ ~ 1])
ls = 10:5:1500

ClTTs, ClEEs, ClTEs = spectrum_cmb([:TT, :EE, :TE], prob, ls)
DlTTs, DlEEs, DlTEs = [SymBoltz.Dl(Cls, ls) for Cls in [ClTTs, ClEEs, ClTEs]]

pTT = plot(ls, DlTTs / 1e-12; ylabel = "Dₗᵀᵀ / 10⁻¹²")
pEE = plot(ls, DlEEs / 1e-12; ylabel = "Dₗᴱᴱ / 10⁻¹²")
pTE = plot(ls, DlTEs / 1e-12; ylabel = "Dₗᵀᴱ / 10⁻¹²", xlabel = "l")
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
prob = CosmologyProblem(M, pars, Dict(M.Λ.Ω₀ => 0.5), [M.g.ℰ ~ 1])
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
