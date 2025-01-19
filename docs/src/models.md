# Models

## Free radiation, matter and cosmological constant (RMΛ)

```@docs
SymBoltz.RMΛ
```

```@example RMΛ
using SymBoltz, Unitful, UnitfulAstro, Plots
M = RMΛ()
pars = Dict(M.r.Ω₀ => 5e-5, M.m.Ω₀ => 0.3, M.g.h => 1.0, M.r.T₀ => NaN) # TODO: don't pass h and T₀ to avoid infinite loop
pars[M.Λ.Ω₀] = 1 - pars[M.r.Ω₀] - pars[M.m.Ω₀]
prob = CosmologyProblem(M, pars)
ks = [1e-3, 1e-2, 1e-1, 1e-0] / u"Mpc"
sol = solve(prob, ks)
p1 = plot(sol, log10(M.g.a), [M.r.ρ, M.m.ρ, M.Λ.ρ, M.G.ρ] ./ M.G.ρ; Nextra = 39) # plot extra points to resolve large time steps
p2 = plot(sol, ks, log10(M.g.a), M.g.Φ; Nextra = 9)
plot(p1, p2, layout = (2, 1), size = (600, 600))
```

## Standard ΛCDM

```@docs
SymBoltz.ΛCDM
```

```@example ΛCDM
using SymBoltz, Plots, Unitful, UnitfulAstro
M = ΛCDM()
prob = CosmologyProblem(M, parameters_Planck18(M), Dict(M.Λ.Ω₀ => 0.5), [M.g.ℰ ~ 1])
pars = merge(pars, shoot(prob, Dict(M.Λ.Ω₀ => 0.5), [M.g.ℰ ~ 1]))
ks = [1e-3, 1e-2, 1e-1, 1e-0] / u"Mpc"
sol = solve(prob, ks)
p1 = plot(sol, log10(M.g.a), [M.γ.ρ, M.ν.ρ, M.h.ρ, M.b.ρ, M.c.ρ, M.Λ.ρ, M.G.ρ] ./ M.G.ρ)
p2 = plot(sol, ks, log10(M.g.a), M.g.Φ)
plot(p1, p2, layout = (2, 1), size = (600, 600))
```

## w₀wₐCDM (CPL parametrization)

```@docs
SymBoltz.w0waCDM
```

```@example w0waCDM
using SymBoltz, Plots, Unitful, UnitfulAstro
M = w0waCDM()
pars = merge(parameters_Planck18(M), Dict( # TODO: make parameters_planck18 take args... of additional parameters
    M.X.w0 => -0.9,
    M.X.wa => 0.2,
    M.X.cₛ² => 1.0
))
ks = [1e-3, 1e-2, 1e-1, 1e-0] / u"Mpc"
prob = CosmologyProblem(M, pars, Dict(M.X.Ω₀ => 0.5), [M.g.ℰ ~ 1])
sol = solve(prob, ks)
p1 = plot(sol, log10(M.g.a), M.X.w)
p2 = plot(sol, ks, log10(M.g.a), M.X.δ)
plot(p1, p2, layout = (2, 1), size = (600, 600))
```

## Brans-Dicke ΛCDM

```@docs
SymBoltz.BDΛCDM
```

Solve background such that `E = G = 1` today, and plot scalar field and Hubble function:
```@example BDΛCDM
using SymBoltz, ModelingToolkit, Unitful, UnitfulAstro, Plots
M = BDΛCDM()
D = Differential(M.t)
ks = [1e-3, 1e-2, 1e-1, 1e-0] / u"Mpc"
pars = merge(parameters_Planck18(M), Dict(M.G.ω => 100.0, D(M.G.ϕ) => 0.0)) # unspecified: M.Λ.Ω₀, M.G.ϕ
prob = CosmologyProblem(M, pars, Dict(M.G.ϕ => 0.95, M.Λ.Ω₀ => 0.5), [M.g.ℰ ~ 1, M.G.G ~ 1])
sol = solve(prob, ks)
p1 = plot(sol, log10(M.g.a), [M.g.ℰ, M.G.G], ylims=(0.8, 1.2))
p2 = plot(sol, ks, log10(M.g.a), M.G.δϕ)
plot(p1, p2, layout = (2, 1), size = (600, 600))
```

## Brans-Dicke RMΛ

```@example BDRMΛ
using SymBoltz, ModelingToolkit, Unitful, UnitfulAstro, Plots
M = SymBoltz.BDRMΛ()
D = Differential(M.t)
pars = Dict(M.r.Ω₀ => 5e-5, M.m.Ω₀ => 0.3, M.g.h => 1.0, M.r.T₀ => 0.0, M.G.ω => 10.0, D(M.G.ϕ) => 0.0) # unspecified: M.Λ.Ω₀, M.G.ϕ
prob = CosmologyProblem(M, pars, Dict(M.G.ϕ => 0.8, M.Λ.Ω₀ => 0.8), [M.g.ℰ ~ 1, M.G.G ~ 1])
k = 1e-0 / u"Mpc"
sol = solve(prob, k; shootopts = (alg = SymBoltz.TrustRegion(), reltol = 1e-2), verbose = true)
p1 = plot(sol, log10(M.g.a), M.G.G)
p2 = plot(sol, k, log10(M.g.a), M.G.δϕ)
plot(p1, p2, layout = (2, 1))
```

## Quintessence-CDM

```@docs
SymBoltz.QCDM
```

```@example QCDM
# TODO: fix!
#=
using SymBoltz, ModelingToolkit, Plots
@parameters V0 N
V = ϕ -> V0 * ϕ^N
M = QCDM(V, h = nothing, I = nothing) # TODO: remove h = nothing
D = Differential(M.t)
pars = merge(parameters_Planck18(M), Dict(M.Q.ϕ => 1, M.Q.V0 => 1e-2, M.Q.N => 2))
prob = CosmologyProblem(M, pars, Dict(D(M.Q.ϕ) => +1.0), [M.g.ℰ ~ 1])
sol = solve(prob)
plot(sol, M.Q.ϕ, M.Q.V, line_z = log10(M.g.a)) # plot V(ϕ(t))
=#
```
