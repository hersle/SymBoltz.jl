# Components

```@autodocs
Modules = [SymBoltz]
Pages   = ["Components.jl"]
```

## Brans-Dicke

```@docs
brans_dicke
BDΛCDM
```

Shoot for parameters that give `E = G = 1` today:
```@example 1
using SymBoltz, ModelingToolkit
M = BDΛCDM()
D = Differential(M.t)

pars_fixed = [parameters_Planck18(M); M.G.ω => 100.0; D(M.G.ϕ) => 0.0] # unspecified: M.Λ.Ω0, M.G.ϕ
pars_guess = [M.G.ϕ => 0.95, M.Λ.Ω0 => 0.7] # initial guesses for shooting method
pars_shoot = shoot(M, pars_fixed, pars_guess, [M.g.ℰ ~ 1, M.G.G ~ 1]; thermo = false, backwards = false) # exact solutions
pars = [pars_fixed; pars_shoot] # merge fixed and shooting parameters
```
Solve background and plot scalar field and Hubble function:
```@example 1
using Plots
sol = solve(M, pars, thermo = false, backwards = false)
plot(sol, log10(M.g.a), [M.g.ℰ, M.G.G], ylims=(0.8, 1.2))
```
