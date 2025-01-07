# Components

```@setup components
# For every component, display collapsible equations as in https://discourse.julialang.org/t/documenter-jl-and-collapsible-markdown/50898/2
using SymBoltz, ModelingToolkit
```

## Spacetimes

```@docs
SymBoltz.metric
```
```@example components
g = SymBoltz.metric()
nothing # hide
```
```@raw html
<details>
<summary><h4 style="display: inline-block">Equations</h4></summary>
```
```@example components
equations(g)
```
```@raw html
</details>
```

## Gravitational theories

```@docs
SymBoltz.general_relativity
```
```@example components
G = SymBoltz.general_relativity(g)
nothing # hide
```
```@raw html
<details>
<summary><h4 style="display: inline-block">Equations</h4></summary>
```
```@example components
equations(G)
```
```@raw html
</details>
<details>
<summary><h4 style="display: inline-block">Initialization equations</h4></summary>
```
```@example components
initialization_equations(G)
```
```@raw html
</details>
```

```@docs
SymBoltz.brans_dicke
```
```@example components
G = SymBoltz.brans_dicke(g)
nothing # hide
```
```@raw html
<details>
<summary><h4 style="display: inline-block">Equations</h4></summary>
```
```@example components
equations(G)
```
```@raw html
</details>
<details>
<summary><h4 style="display: inline-block">Initialization equations</h4></summary>
```
```@example components
initialization_equations(G)
```
```@raw html
</details>
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
nothing # hide
```
```@raw html
<details>
<summary><h4 style="display: inline-block">Equations</h4></summary>
```
```@example components
equations(s)
```
```@raw html
</details>
<details>
<summary><h4 style="display: inline-block">Initialization equations</h4></summary>
```
```@example components
initialization_equations(s)
```
```@raw html
</details>
```

### Relativistic

```@docs
SymBoltz.radiation
```
```@example components
r = SymBoltz.radiation(g)
nothing # hide
```
```@raw html
<details>
<summary><h4 style="display: inline-block">Equations</h4></summary>
```
```@example components
equations(r)
```
```@raw html
</details>
<details>
<summary><h4 style="display: inline-block">Initialization equations</h4></summary>
```
```@example components
initialization_equations(r)
```
```@raw html
</details>
```

```@docs
SymBoltz.photons
```
```@example components
γ = SymBoltz.photons(g; lmax = 6)
nothing # hide
```
```@raw html
<details>
<summary><h4 style="display: inline-block">Equations</h4></summary>
```
```@example components
equations(γ)
```
```@raw html
</details>
<details>
<summary><h4 style="display: inline-block">Initialization equations</h4></summary>
```
```@example components
initialization_equations(γ)
```
```@raw html
</details>
```

### Non-relativistic

```@docs
SymBoltz.matter
```
```@example components
m = SymBoltz.matter(g)
nothing # hide
```
```@raw html
<details>
<summary><h4 style="display: inline-block">Equations</h4></summary>
```
```@example components
equations(m)
```
```@raw html
</details>
<details>
<summary><h4 style="display: inline-block">Initialization equations</h4></summary>
```
```@example components
initialization_equations(m)
```
```@raw html
</details>
```

```@docs
SymBoltz.cold_dark_matter
```
```@example components
c = SymBoltz.cold_dark_matter(g)
nothing # hide
```
```@raw html
<details>
<summary><h4 style="display: inline-block">Equations</h4></summary>
```
```@example components
equations(c)
```
```@raw html
</details>
<details>
<summary><h4 style="display: inline-block">Initialization equations</h4></summary>
```
```@example components
initialization_equations(c)
```
```@raw html
</details>
```

```@docs
SymBoltz.baryons
```
```@example components
b = SymBoltz.baryons(g)
nothing # hide
```
```@raw html
<details>
<summary><h4 style="display: inline-block">Equations</h4></summary>
```
```@example components
equations(b)
```
```@raw html
</details>
<details>
<summary><h4 style="display: inline-block">Initialization equations</h4></summary>
```
```@example components
initialization_equations(b)
```
```@raw html
</details>
```

### Neutrinos

```@docs
SymBoltz.massless_neutrinos
```
```@example components
ν = SymBoltz.massless_neutrinos(g; lmax = 6)
nothing # hide
```
```@raw html
<details>
<summary><h4 style="display: inline-block">Equations</h4></summary>
```
```@example components
equations(ν)
```
```@raw html
</details>
<details>
<summary><h4 style="display: inline-block">Initialization equations</h4></summary>
```
```@example components
initialization_equations(ν)
```
```@raw html
</details>
```

```@docs
SymBoltz.massive_neutrinos
```
```@example components
h = SymBoltz.massive_neutrinos(g; nx=4, lmax = 5)
nothing # hide
```
```@raw html
<details>
<summary><h4 style="display: inline-block">Equations</h4></summary>
```
```@example components
equations(h)
```
```@raw html
</details>
<details>
<summary><h4 style="display: inline-block">Initialization equations</h4></summary>
```
```@example components
initialization_equations(h)
```
```@raw html
</details>
```

### Dark energy

```@docs
SymBoltz.cosmological_constant
```
```@example components
Λ = SymBoltz.cosmological_constant(g)
nothing # hide
```
```@raw html
<details>
<summary><h4 style="display: inline-block">Equations</h4></summary>
```
```@example components
equations(Λ)
```
```@raw html
</details>
<details>
<summary><h4 style="display: inline-block">Initialization equations</h4></summary>
```
```@example components
initialization_equations(Λ)
```
```@raw html
</details>
```

```@docs
SymBoltz.w0wa
```
```@example components
X = SymBoltz.w0wa(g)
nothing # hide
```
```@raw html
<details>
<summary><h4 style="display: inline-block">Equations</h4></summary>
```
```@example components
equations(X)
```
```@raw html
</details>
<details>
<summary><h4 style="display: inline-block">Initialization equations</h4></summary>
```
```@example components
initialization_equations(X)
```
```@raw html
</details>
```

```@docs
SymBoltz.quintessence
```
```@example components
@variables ϕ(t) V(ϕ) V′(ϕ) V′′(ϕ)
Q = SymBoltz.quintessence(g, V, V′, V′′)
nothing # hide
```
```@raw html
<details>
<summary><h4 style="display: inline-block">Equations</h4></summary>
```
```@example components
equations(Q)
```
```@raw html
</details>
<details>
<summary><h4 style="display: inline-block">Initialization equations</h4></summary>
```
```@example components
initialization_equations(Q)
```
```@raw html
</details>
```
