# Observable quantities

This page shows observable quantities that can be derived from solutions of the Einstein-Boltzmann system, such as power spectra and distances.

## Primordial power spectra

```@docs
SymBoltz.spectrum_primordial
```

#### Example

```@example
using SymBoltz, Plots
M = ΛCDM()
pars = Dict(M.g.h => 0.7, M.I.As => 2e-9, M.I.ns => 0.95)
ks = 10 .^ range(-2, +4, length=100)
Ps = spectrum_primordial(ks, M, pars)
plot(log10.(ks), log10.(Ps); xlabel = "log10(k / (H₀/c))", ylabel = "log10(P / (c/H₀)³)")
```

## Matter power spectra

```@docs
SymBoltz.spectrum_matter
SymBoltz.spectrum_matter_nonlinear
```

#### Example

With explicitly chosen wavenumbers:

```@example matter
using SymBoltz, Plots
M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
ks = 10 .^ range(-2, +5, length=100)
sol = solve(prob, ks)

# Linear power spectrum
modes = [:m, :cb, :h]
Ps = spectrum_matter(modes, sol, ks)
plot(
    log10.(ks), transpose(log10.(Ps));
    xlabel = "log10(k / (H₀/c))", ylabel = "log10(P / (c/H₀)³)",
    label = permutedims("linear (SymBoltz), " .* string.(modes)), ylims = (-15, -6),
    linestyle = [:solid :dash :dot :dashdot :dashdotdot], legend_position = :bottomleft
)

# Nonlinear power spectrum (from halofit)
Ps = spectrum_matter_nonlinear(sol, ks)
plot!(
    log10.(ks), log10.(Ps);
    label = "non-linear (halofit), matter", legend_position = :bottomleft
)
```

As a function of conformal time and redshift:

```@example matter
τs = range(sol[M.τ][end], 0.5, length=10)
Ps = spectrum_matter(sol, ks, τs)
zs = sol(M.g.z, τs) # corresponding redshifts
plot(
    log10.(ks), transpose(log10.(Ps));
    xlabel = "log10(k / (H₀/c))", ylabel = "log10(P / (c/H₀)³)",
    label = permutedims("z=".*string.(round.(zs;digits=1))),
    legend_position = :bottomleft
)
```

## CMB power spectra

```@docs
SymBoltz.spectrum_cmb
```

#### Example

```@example
using SymBoltz, Plots
M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)

ls = 25:25:3000 # 25, 50, ..., 3000
jl = SphericalBesselCache(ls)
modes = [:TT, :EE, :TE, :ψψ, :ψT, :ψE]
Dls = spectrum_cmb(modes, prob, jl; normalization = :Dl)

plot(ls, log10.(abs.(Dls)); xlabel = "l", ylabel = "lg(Dₗ)", label = permutedims(String.(modes)))
```

## Two-point correlation function

```@docs
SymBoltz.correlation_function
```

#### Example

```@example
using SymBoltz, Plots
M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
ks = 10 .^ range(-2, +6, length=300)
sol = solve(prob, ks)
rs, ξs = correlation_function(sol)
rs = rs * L100 # convert from c/H₀ to Mpc/h
plot(rs, @. ξs * rs^2; xlims = (0, 200), xlabel = "r / (Mpc/h)", ylabel = "r² ξ / (Mpc/h)²")
```

## Matter density fluctuations

```@docs
SymBoltz.variance_matter
SymBoltz.stddev_matter
```

```@example
using SymBoltz, Plots
M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
ks = 10 .^ range(-2, +6, length=300)
sol = solve(prob, ks)

Rs = 10 .^ range(0.5, 2.5, length=100) # Mpc/h
σs = stddev_matter.(sol, Rs / L100) # Mpc/h to c/H₀
plot(log10.(Rs), log10.(σs); xlabel = "lg(R / (Mpc/h))", ylabel = "lg(σ)", label = nothing)

σ8 = stddev_matter(sol, 8 / L100) # 8 Mpc/h to c/H₀
scatter!((log10(8), log10(σ8)), series_annotation = text("  σ₈ = $(round(σ8; digits=3))", :left), label = nothing)
```

## Distance measures

```@docs
distance
```

```@example
using SymBoltz, Plots
M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
sol = solve(prob)

τ0 = today(sol)
τs = range(0.5*τ0, τ0, length = 100) # conformal times back in time
zs = sol(M.g.z, τs) # corresponding redshifts
modes = [:χ, :M, :A, :L, :V] # distances from today
Ds = distance(modes, sol, τs)
labels = ["Dχ (lookback distance)" "DM (transverse comoving distance)" "DA (angular diameter distance)" "DL (luminosity distance)" "DV (volume-averaged distance)"]
plot(zs, transpose(Ds); xlabel = "z", ylabel = "D / (c/H₀)", label = labels, xlims = extrema(zs), ylims = (0, 3))
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
```

```@example
using SymBoltz, Plots, DataInterpolations
M = ΛCDM(h = nothing, ν = nothing)
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
sol = solve(prob)

τs = sol[M.τ] # conformal times in background solution
ks = exp.(range(log(1.0), log(2000.0), length = 50)) # logarithmic k-grid
Ss = source_grid(prob, M.ST, τs, ks)
iτ = argmax(sol[M.b.v]) # index of decoupling time
iτs = iτ-75:iτ+75 # indices around decoupling
p1 = surface(ks, τs[iτs], Ss[iτs, :]; camera = (45, 25), xlabel = "k", ylabel = "τ", zlabel = "S", colorbar = false)

lgas = -6.0:0.2:0.0
τs = LinearInterpolation(sol[M.τ], sol[log10(M.g.a)])(lgas) # τ at given lg(a)
ks = 5.0:5.0:100.0
Ss = source_grid(prob, M.g.Ψ, τs, ks)
p2 = wireframe(ks, lgas, Ss; camera = (75, 20), xlabel = "k", ylabel = "lg(a)", zlabel = "Φ")

plot(p1, p2)
```

## Line-of-sight integration

```@docs
SymBoltz.los_integrate
```
