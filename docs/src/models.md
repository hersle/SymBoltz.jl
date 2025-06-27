# Cosmologies (full models)

## Free radiation, matter and cosmological constant (RMΛ)

```@docs
SymBoltz.RMΛ
```

```@example RMΛ
using SymBoltz, Unitful, UnitfulAstro, Plots
M = RMΛ()
pars = Dict(M.r.Ω₀ => 5e-5, M.m.Ω₀ => 0.3, M.g.h => 1.0, M.r.T₀ => NaN) # TODO: don't pass h and T₀ to avoid infinite loop
prob = CosmologyProblem(M, pars)
ks = [1e-3, 1e-2, 1e-1, 1e-0] / u"Mpc"
sol = solve(prob, ks)
p1 = plot(sol, log10(M.g.a), [M.r.ρ, M.m.ρ, M.Λ.ρ, M.G.ρ] ./ M.G.ρ)
p2 = plot(sol, log10(M.g.a), M.g.Φ, ks)
plot(p1, p2, layout = (2, 1), size = (600, 600))
```

## Standard ΛCDM

```@docs
SymBoltz.ΛCDM
```

```@example ΛCDM
using SymBoltz, Plots, Unitful, UnitfulAstro
M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
ks = [1e-3, 1e-2, 1e-1, 1e-0] / u"Mpc"
sol = solve(prob, ks)
p1 = plot(sol, log10(M.g.a), [M.γ.ρ, M.ν.ρ, M.h.ρ, M.b.ρ, M.c.ρ, M.Λ.ρ, M.G.ρ] ./ M.G.ρ)
p2 = plot(sol, log10(M.g.a), M.g.Φ, ks)
plot(p1, p2, layout = (2, 1), size = (600, 600))
```

## w₀wₐCDM (CPL parametrization)

```@docs
SymBoltz.w0waCDM
```

```@example w0waCDM
using SymBoltz, Plots, Unitful, UnitfulAstro
M = w0waCDM()
pars = merge(parameters_Planck18(M), Dict(
    M.X.w0 => -0.9,
    M.X.wa => 0.2,
    M.X.cₛ² => 1.0
))
ks = [1e-3, 1e-2, 1e-1, 1e-0] / u"Mpc"
prob = CosmologyProblem(M, pars)
sol = solve(prob, ks)
p1 = plot(sol, log10(M.g.a), M.X.w)
p2 = plot(sol, log10(M.g.a), M.X.δ, ks)
plot(p1, p2, layout = (2, 1), size = (600, 600))
```

## Brans-Dicke ΛCDM

```@docs
SymBoltz.BDΛCDM
```

Solve background such that `E = G = 1` today, and plot scalar field and Hubble function:
```@example BDΛCDM
using SymBoltz, Unitful, UnitfulAstro, Plots
M = BDΛCDM()
D = Differential(M.τ)
ks = [1e-3, 1e-2, 1e-1, 1e-0] / u"Mpc"
pars = merge(parameters_Planck18(M), Dict(M.G.ω => 100.0, D(M.G.ϕ) => 0.0)) # unspecified: M.Λ.Ω₀, M.G.ϕ
prob = CosmologyProblem(M, pars, Dict(M.G.ϕ => 0.95, M.Λ.Ω₀ => 0.5), [M.g.ℰ ~ 1, M.G.G ~ 1])
sol = solve(prob, ks; verbose = true)
p1 = plot(sol, log10(M.g.a), [M.g.ℰ, M.G.G], ylims = (0.8, 1.2))
p2 = plot(sol, log10(M.g.a), M.G.δϕ, ks)
plot(p1, p2, layout = (2, 1), size = (600, 600))
```

## Brans-Dicke RMΛ

```@example BDRMΛ
using SymBoltz, Unitful, UnitfulAstro, Plots
M = SymBoltz.BDRMΛ()
D = Differential(M.τ)
pars = Dict(M.r.Ω₀ => 5e-5, M.m.Ω₀ => 0.3, M.g.h => 1.0, M.r.T₀ => 0.0, M.G.ω => 10.0, D(M.G.ϕ) => 0.0) # unspecified: M.Λ.Ω₀, M.G.ϕ
prob = CosmologyProblem(M, pars, Dict(M.G.ϕ => 0.8, M.Λ.Ω₀ => 0.8), [M.g.ℰ ~ 1, M.G.G ~ 1])
k = 1e-0 / u"Mpc"
sol = solve(prob, k; verbose = true)
p1 = plot(sol, log10(M.g.a), M.G.G)
p2 = plot(sol, log10(M.g.a), M.G.δϕ, k)
plot(p1, p2, layout = (2, 1))
```

## Quintessence-CDM

```@docs
SymBoltz.QCDM
```

```@example QCDM
using SymBoltz, Plots
@parameters V0 N
V = ϕ -> V0 * ϕ^N
M = QCDM(V, I = nothing)
D = Differential(M.τ)
pars = merge(parameters_Planck18(M), Dict(M.Q.ϕ => 1, D(M.Q.ϕ) => 1.0, M.Q.V0 => 1e-2, M.Q.N => 2))
prob = CosmologyProblem(M, pars)
sol = solve(prob) # TODO: shoot so M.g.ℰ ~ 1 today
plot(sol, M.Q.ϕ, M.Q.V, line_z = log10(M.g.a)) # plot V(ϕ(τ))
```
