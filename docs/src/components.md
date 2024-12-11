# Components

```@setup components
# TODO: collapsible equations? https://discourse.julialang.org/t/documenter-jl-and-collapsible-markdown/50898/2
using SymBoltz, ModelingToolkit
```

## Spacetimes

```@docs
SymBoltz.metric
```
```@example components
g = SymBoltz.metric()
```

## Gravitational theories

```@docs
SymBoltz.general_relativity
```
```@example components
G = SymBoltz.general_relativity(g)
```

```@docs
SymBoltz.brans_dicke
```
```@example components
G = SymBoltz.brans_dicke(g)
```

## Species

### Generic

```@docs
SymBoltz.species_constant_eos
```
```@example components
using SymBoltz: t
@variables w(t)
s = SymBoltz.species_constant_eos(g, w)
```

### Relativistic

```@docs
SymBoltz.radiation
```
```@example components
r = SymBoltz.radiation(g)
```

```@docs
SymBoltz.photons
```
```@example components
γ = SymBoltz.photons(g; lmax = 6)
```

### Non-relativistic

```@docs
SymBoltz.matter
```
```@example components
m = SymBoltz.matter(g)
```

```@docs
SymBoltz.cold_dark_matter
```
```@example components
c = SymBoltz.cold_dark_matter(g)
```

```@docs
SymBoltz.baryons
```
```@example components
b = SymBoltz.baryons(g)
```

### Neutrinos

```@docs
SymBoltz.massless_neutrinos
```
```@example components
ν = SymBoltz.massless_neutrinos(g; lmax = 6)
```

```@docs
SymBoltz.massive_neutrinos
```
```@example components
h = SymBoltz.massive_neutrinos(g; nx=4, lmax = 5)
```

### Dark energy

```@docs
SymBoltz.cosmological_constant
```
```@example components
Λ = SymBoltz.cosmological_constant(g)
```

```@docs
SymBoltz.w0wa
```
```@example components
X = SymBoltz.w0wa(g)
```

```@docs
SymBoltz.quintessence
```
```@example components
@variables ϕ(t) V(ϕ) V′(ϕ) V′′(ϕ)
Q = SymBoltz.quintessence(g, V, V′, V′′)
```
