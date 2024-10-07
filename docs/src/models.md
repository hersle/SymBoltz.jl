# Models

## RMΛ

```@docs
SymBoltz.RMΛ
```

```@example RMΛ
using SymBoltz, Unitful, UnitfulAstro, Plots
M = RMΛ()
pars = [M.r.Ω0 => 5e-5, M.m.Ω0 => 0.3, M.g.h => 1.0, M.r.T0 => 0.0] # TODO: don't pass h and T0 to avoid infinite loop
ks = [1e-3, 1e-2, 1e-1, 1e-0] / u"Mpc"
sol = solve(M, pars, ks)
p1 = plot(sol, log10(M.g.a), [M.r.ρ, M.m.ρ, M.Λ.ρ, M.G.ρ] ./ M.G.ρ; N = 10000)
p2 = plot(sol, ks, log10(M.g.a), M.g.Φ; N = 10000)
plot(p1, p2, layout = (2, 1), size = (600, 600))
```

## ΛCDM

```@docs
SymBoltz.ΛCDM
```

```@example ΛCDM
using SymBoltz, Plots, Unitful, UnitfulAstro
M = ΛCDM()
pars = parameters_Planck18(M)
ks = [1e-3, 1e-2, 1e-1, 1e-0] / u"Mpc"
sol = solve(M, pars, ks)
p1 = plot(sol, log10(M.g.a), [M.γ.ρ, M.ν.ρ, M.h.ρ, M.b.ρ, M.c.ρ, M.Λ.ρ, M.G.ρ] ./ M.G.ρ)
p2 = plot(sol, ks, log10(M.g.a), M.g.Φ)
plot(p1, p2, layout = (2, 1), size = (600, 600))
```

## Brans-Dicke ΛCDM

```@docs
SymBoltz.BDΛCDM
```

Shoot for parameters that give `E = G = 1` today:
```@example BDΛCDM
using SymBoltz, ModelingToolkit
M = BDΛCDM()
D = Differential(M.t)

pars_fixed = [parameters_Planck18(M); M.G.ω => 100.0; D(M.G.ϕ) => 0.0] # unspecified: M.Λ.Ω0, M.G.ϕ
pars_guess = [M.G.ϕ => 0.95, M.Λ.Ω0 => 0.7] # initial guesses for shooting method
pars_shoot = shoot(M, pars_fixed, pars_guess, [M.g.ℰ ~ 1, M.G.G ~ 1]; thermo = false, backwards = false) # exact solutions
pars = [pars_fixed; pars_shoot] # merge fixed and shooting parameters
```
Solve background and plot scalar field and Hubble function:
```@example BDΛCDM
# TODO: solve for perturbed scalar field
# TODO: (must add proper pressure perturbation to all species first)
using Unitful, UnitfulAstro, Plots
ks = [1e-0] / u"Mpc"
sol = solve(M, pars, ks, backwards = false)
p1 = plot(sol, log10(M.g.a), [M.g.ℰ, M.G.G], ylims=(0.85, 1.15))
p2 = plot(sol, ks, log10(M.g.a), M.G.δϕ)
plot(p1, p2, layout = (2, 1), size = (600, 600))
```

## Brans-Dicke RMΛ

```@example BDRMΛ
using SymBoltz, ModelingToolkit, Unitful, UnitfulAstro, Plots
M = SymBoltz.BDRMΛ()
D = Differential(M.t)

pars_fixed = [M.r.Ω0 => 5e-5, M.m.Ω0 => 0.3, M.g.h => 1.0, M.r.T0 => 0.0, M.G.ω => 10.0, D(M.G.ϕ) => 0.0] # unspecified: M.Λ.Ω0, M.G.ϕ
pars_guess = [M.G.ϕ => 0.95, M.Λ.Ω0 => 0.7] # initial guesses for shooting method
pars_shoot = shoot(M, pars_fixed, pars_guess, [M.g.ℰ ~ 1, M.G.G ~ 1]; thermo = false, backwards = false) # exact solutions
pars = [pars_fixed; pars_shoot] # merge fixed and shooting parameters

# TODO: add perturbations to ΛCDM example, too (need to define h.cs²)
ks = [1e-0] / u"Mpc"
sol = solve(M, pars, ks, backwards = false) # TODO: set background integration direction during model creation
p1 = plot(sol, log10(M.g.a), M.G.G; N = 10000)
p2 = plot(sol, ks[1], log10(M.g.a), M.G.δϕ; N = 50000)
plot(p1, p2, layout = (2, 1))
```

## Quintessence-CDM

```@docs
SymBoltz.QCDM
```

```@example QCDM
using SymBoltz, ModelingToolkit, Plots
@parameters V0 N
V(ϕ) = V0 * ϕ^N
M = QCDM(V)
D = Differential(M.t)
pars = [parameters_Planck18(M); M.Q.ϕ => 1; M.Q.V0 => 1e-2; M.Q.N => 2]
sol = solve(M, pars, thermo = false; guesses = [D(M.Q.ϕ) => +1.0])
plot(sol, M.Q.ϕ, M.Q.V, line_z = log10(M.g.a)) # plot V(ϕ(t))
```
