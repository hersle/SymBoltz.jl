---
title: |
    | ![](media/einstein_boltzmann.png){height=55%}
    | SymBoltz.jl: A symbolically modular, approximation-free and differentiable Einstein-Boltzmann solver
    | 
    | \scriptsize ![](media/github.png){height=0.75em} \textcolor{black}{Code \& documentation:} [github.com/hersle/SymBoltz.jl](https://github.com/hersle/SymBoltz.jl)
#titlegraphic: media/components_baumann.pdf
colorlinks: true
header-includes: |
    \usepackage{emoji}
    \usepackage{tikz}
    \usetikzlibrary{positioning}
    \usetikzlibrary{arrows}
    \usetikzlibrary{trees}
    \usepackage{mathtools}
    \usepackage{caption}
    \captionsetup[figure]{labelformat=empty}
pdf-engine: lualatex
monofont: JuliaMono
---

---

# Outline

1. Current state of cosmology

2. Boltzmann solvers

3. New features

   1. Symbolic modularity

   2. Approximation-free

   3. Differentiable

4. Why start from scratch in Julia?

5. Future

# Cosmology is at a crossroads

![ ](media/crossroads_meme_edited.jpg){height=75%}

# ΛCDM model: components

\begin{center}
\begin{tikzpicture}[
   font=\scriptsize,
   comp/.style = {draw, circle, fill=gray, minimum size = 2.0cm, align=center},
   grow cyclic,
   level 1/.append style = {level distance = 3.3cm, sibling angle = 60},
   interaction/.style = {latex-latex, thick},
]
\node[comp] (grav) {Gravity} [interaction]
child { node[comp] (bar) {Baryons}}
child { node[comp] (pho) {Photons}}
child { node[comp] (neu) {Massless\\neutrinos}}
child { node[comp] (mneu) {Massive\\neutrinos}}
child { node[comp] (cdm) {Cold dark\\matter}}
child { node[comp] (cc) {Cosmological\\constant}};
\draw[interaction] (bar) -- node[above=, sloped, align=center] {Thomson\\scattering} (pho);
\end{tikzpicture}
\end{center}

# ΛCDM model: challenges

- [$H₀$ and $S₈$ tensions](https://arxiv.org/abs/2105.05208)

- [$Λ$ problem](https://arxiv.org/abs/2105.05208)

- [Anisotropically distributed quasars](https://arxiv.org/abs/2009.14826)

- [Surprisingly massive early galaxies (JWST)](https://arxiv.org/abs/2309.13100)

- [Excess ISW signal from supervoids](https://arxiv.org/abs/1811.07812) 

- [DESI observations](https://arxiv.org/pdf/2404.08056)

Cosmologists are exploring alternative models.

# Alternative models

![Beyond ΛCDM panel, Oslo 2015 ([1512.05356](https://arxiv.org/abs/1512.05356)).](media/vote_edited.png){width=95%}

# Alternative models: e.g. modified gravity

[![](https://www.tessabaker.space/images/map_slide_v2.pdf){width=95%}](https://www.tessabaker.space/images/map_slide_v2.pdf)

# Boltzmann solvers are fundamental cosmology tools

[![Used in most cosmological analyses](media/hammer_edited.jpeg)](https://www.reddit.com/r/photoshopbattles/comments/cfxzzg/battle_366_bigger_splash_via_previous_winner/)

# An Einstein-Boltzmann solver typically ...

1. reads input parameters

2. solves background ODEs

   - $H^2 = \frac{8 \pi G}{3} \rho$

3. solves thermodynamics ODEs

   - $\frac{\mathrm{d} x_e}{\mathrm{d} t} = \ldots$

4. solves perturbation ODEs

   - $\dot{Φ} = -\frac{4\pi}{3} a^2 ℋ \delta\rho - \frac{k^2}{3 ℋ} Φ - ℋ Ψ$

5. solves line of sight ODEs

6. computes spectra

7. writes output

# An Einstein-Boltzmann solver ...
... solves the (gravitational) Einstein equations
$$
%\underset{\underset{\mathllap{\text{spacetime geometry}}}{\uparrow}}{G_{\mu \nu}} = \frac{8 \pi G}{c^4} \underset{\underset{\mathrlap{\text{energy-momentum content}}}{\uparrow}}{T_{\mu \nu}}
\scriptsize{\textcolor{gray}{\text{geometry} \rightarrow \quad}} G_{\mu\nu} = \frac{8 \pi G}{c^4} T_{\mu\nu} \scriptsize{\textcolor{gray}{\quad \leftarrow \text{content}}}
$$
coupled to particle species, each described by Boltzmann equations
$$\left( p^\mu \frac{\partial}{\partial x^\mu} - \Gamma^\mu_{\alpha \beta} p^\alpha p^\beta \frac{\partial}{\partial p^\mu} \right) \underset{\textcolor{gray}{\underset{\mathclap{\text{distribution function}}}{\uparrow}}}{f_s} = C[f_s] \textcolor{gray}{\quad \leftarrow \text{interactions}}$$
to 1st perturbative order around a homo. and iso. FLRW spacetime
$$\mathrm{d}s^2 = a^2 (-1-2 \Psi) \mathrm{d}t^2 + a^2 (1 - 2 \Phi) \mathrm{d}\mathbf{x}^2$$
(with lots of pre- and post-processing).

# An Einstein-Boltzmann solver solves the ODEs ...

$$
\begin{aligned}
%\dot{\phantom{x}} &= \mathrm{d}/\mathrm{d} \tau, \\
\dot{a} &= \sqrt{\frac{8 \pi}{3} \rho} \, a^2, \\
\ddot{a} &= \frac{\dot{a}^2}{a} - \frac{4\pi}{3} (\rho + 3P) a^3, \\
\dot{Φ} &= -\frac{4\pi}{3} a^2 ℋ \delta\rho - \frac{k^2}{3 ℋ} Φ - ℋ Ψ, \\
k^2 (Φ - Ψ) &= 12\pi a^2 Π, \\
\dot{ρ}_s &= -3 ℋ (ρ_s + P_s), \\
P_s &= w_s \rho_s, \\
\dot{δ}_s &= -(1+w_s) (θ_s-3\dot{Φ}) - 3 ℰ (cₛ²-w_s)δ_s, \\
\dot{θ}_s &= -ℋ(1-3w_s)θ_s - \frac{\dot{w}_s θ_s}{1+w_s} + \frac{cₛ² k^2 δ_s}{1+w_s}  - k^2 σ_s + k^2 Ψ, \\
&\vdots
\end{aligned}
$$

---

# An Einstein-Boltzmann solver outputs e.g. ...

:::::::::::::: {.columns}
::: {.column width="50%"}
![Matter power spectrum](media/matter_power_spectrum.png){height=3.5cm}
:::
::: {.column width="50%"}
![CMB power spectrum](media/cmb_power_spectrum.png){height=3.5cm}
:::
::::::::::::::


# History of Boltzmann solvers

\scriptsize

Year                  Code                                                                                                 Lang.   New features                                     
--------------------- ---------------------------------------------------------------------------------------------------- ------- ------------------------------------------------- ------
1995                   [COSMICS](https://arxiv.org/abs/astro-ph/9506072) [ ](https://arxiv.org/pdf/astro-ph/9506070)       Fortran First proper treatment; seminal paper             \emoji{headstone}
1996                   [CMBFAST](https://arxiv.org/abs/astro-ph/9603033)                                                   Fortran Line-of-sight integration (lower $l_\text{max}$)  \emoji{headstone}
2000                   [CAMB](https://arxiv.org/abs/astro-ph/9911177)                                                      Fortran Further development, closed models                \emoji{sports-medal}
2003                   [CMBEASY](https://arxiv.org/abs/astro-ph/0302138)                                                   C++     Code structure, object-oriented                   \emoji{headstone}
2011                   [CLASS](https://arxiv.org/abs/1104.2932)                                                            C/Py    User-friendliness, flexibility, accuracy control  \emoji{sports-medal}
2017                   [PyCosmo](https://arxiv.org/abs/1708.05177)                                                         Py/C++  Symbolics, C++ code gen., sparsity optim.         \emoji{hatching-chick}
2021                   [Bolt](https://github.com/xzackli/Bolt.jl/)                                                         Julia   Differentiable, approx.-free?                     \emoji{hatching-chick}
2024                   [DISCO-EB](https://arxiv.org/abs/2311.03291)                                                        Py/Jax  Differentiable, approx.-free?                     \emoji{hatching-chick}
2025                   [SymBoltz](https://github.com/hersle/SymBoltz.jl)?                                                  Julia   Symbolic modularity, approx.-free, differentiable \emoji{hatching-chick}

... and all forks thereof; e.g. [EFTCAMB](https://arxiv.org/abs/1312.5742), [HiCLASS](https://arxiv.org/abs/1909.01828), ...

# Feature 1: symbolic modularity

\scriptsize

```
julia> M = ΛCDM(); hierarchy(M)
ΛCDM: Standard cosmological constant and cold dark matter cosmological model
├─ g: Spacetime FLRW metric in Newtonian gauge
├─ G: General relativity gravity
├─ γ: Photon radiation
├─ ν: Massless neutrinos
├─ c: Cold dark matter
├─ b: Baryonic matter
│  └─ rec: Baryon-photon recombination thermodynamics (RECFAST)
├─ h: Massive neutrino
└─ Λ: Cosmological constant
```


# Add a new species: CLASS

Official advice: CTRL+F for `"species"`

TODO: add picture of `grep _fld`?

- Add input parameter hooks in `input.c`
- Define indices `background.c` and `perturbations.c` (maybe `thermodynamics.c`)
- Write down governing equations in `background.c` and `perturbations.c` (maybe `thermodynamics.c`)
- Define columns and save output in `background.c` and `perturbations.c` (maybe `thermodynamics.c`)
- Save output in `background.c`
- Handle all combos of independent/dependent input parameters (e.g. shooting) in `input.c`
- Always wrap in `if (has_species == _TRUE_) {...}`
- Add input parameter hooks to Python wrapper, too.

Observation: most of these should be common to a modeling library

# Add a new species: SymBoltz

Everything related to *one species* should be in *one place*:

\tiny

```julia
M = ΛCDM()

g = M.g
pars = @parameters w₀ wₐ ρ₀ Ω₀ cₛ²
vars = @variables ρ(t) P(t) w(t) δ(t) θ(t) σ(t)
eqs = [
    O(ϵ^0)(w ~ w₀ + wₐ * (1 - g.a))
    O(ϵ^0)(ρ ~ ρ₀ * g.a^(-3 * (1 + w₀ + wₐ)) * exp(-3 * wₐ * (1 - g.a)))
    O(ϵ^0)(P ~ w * ρ)
    O(ϵ^1)(D(δ) ~ -(1 + w) * (θ - 3*g.Φ) - 3 * g.ℰ * (cₛ² - w) * δ)
    O(ϵ^1)(D(θ) ~ -g.ℰ * (1 - 3*w) - D(w) / (1 + w) * θ + cₛ² / (1 + w) * k^2 * δ - k^2 * σ + k^2 * g.Ψ)
    O(ϵ^1)(σ ~ 0)
]
initialization_eqs = [
    O(ϵ^1)(δ ~ -3/2 * (1+w) * g.Ψ)
    O(ϵ^1)(θ ~ 1/2 * (k^2*t) * g.Ψ)
]
defaults = [
    ρ₀ => 3/8π * Ω₀
]
@named X = ODESystem(eqs, t, vars, pars; initialization_eqs, defaults)

M′ = ΛCDM(Λ = X; name = :w₀wₐCDM)
```

# Architecture: monolithic vs. modular

:::::::::::::: {.columns}
::: {.column width="50%"}
![monolithic](media/monolithic.png)
:::
::: {.column width="50%"}
![modular](media/modular.png)
:::
::::::::::::::

# Architecure: SymBoltz / modular

:::::::::::::: {.columns}
::: {.column width="50%"}

![](media/organization_class.png){height=50%}

:::
::: {.column width="50%"}

\begin{tikzpicture}[
   font=\scriptsize,
   comp/.style = {draw, circle, fill=gray, minimum size = 2.0cm, align=center},
   grow cyclic,
   level 1/.append style = {level distance = 3.3cm, sibling angle = 60},
   interaction/.style = {latex-latex, thick},
]
\node[comp] (grav) {Gravity} [interaction]
child { node[comp] (bar) {Baryons}}
child { node[comp] (pho) {Photons}}
child { node[comp] (neu) {Massless\\neutrinos}}
child { node[comp] (mneu) {Massive\\neutrinos}}
child { node[comp] (cdm) {Cold dark\\matter}}
child { node[comp] (cc) {Cosmological\\constant}};
\draw[interaction] (bar) -- node[above=, sloped, align=center] {Thomson\\scattering} (pho);
\end{tikzpicture}

:::
::::::::::::::


# In my *opinion*:

- CLASS' monolithic structure doesn't scale well in model space

  - massive number of forks leads to fragmentation (TODO: fork graph?)


- User should be able to just write down the full set of (background + perturbation) equations
- Code should automatically read input and write output
- Code should automatically transform full equations into necessary computational stages (background, thermodynamics, perturbations, LOS, ...)


# Feature 2: approximation-free

![Problem: Einstein-Boltzmann ODEs are *very* stiff! [ ](https://peakptandwellness.com/blog/16707/Helpful-Stretches-to-Combat-Stiff-Legs-from-Your-Friendly-Physical-Therapist)](media/stiffness_edited.png){width=60%}

![](media/stiff_ode.png){width=60%}

# Solution 1: explicit integrators with approximations

:::::::::::::: {.columns}
::: {.column width="61%"}
Remove stiffness with approximations for

- baryon-photon tight coupling (TCA),

- radiation streaming (RSA),

- ultra-relativistic fluid (UFA),

- analog for massive neutrinos,

- analog for interacting dark radiation,

- analog for other extended models.
:::
::: {.column width="39%"}
![[CLASS approximations](https://arxiv.org/pdf/1104.2933)](media/approximations.png){height=50%}
:::
::::::::::::::

\emoji{green-square} Fast

\emoji{red-square} Scales poorly in model space, complicates physics and numerics

# Solution 2: implicit integrators without approximations

General Runge-Kutta method for ODE $\dot{y} = f(t, y)$ with $s$ stages:
$$\textstyle y_{n+1} = y_n + h \sum_{i=1}^s b_i k_i, \quad k_i = f(t_n + c_i h, y_n + \sum_{j=1}^s a_{ij} k_j)$$

- $A = [a_{ij}]$ lower triangular $\implies$ $\mathbf{k}$ given **explicitly**.

- $A = [a_{ij}]$ general: $\implies$ $\mathbf{k}$ given **implicitly** from $\mathbf{z}$ in
  $$\underbrace{\left[ I - (A \otimes I) J(\mathbf{z}_n) \right]}_\text{$sd \times sd$ matrix} \mathbf{z}_{n+1} = -\mathbf{z}_n + (A \otimes I) F(\mathbf{z}_n)$$
  - Need accurate $J_{ij} = \partial f_i / \partial y_j$: automatic differentiation > finite differences

  - Reduce cost by LU factorizing and reusing $I - (A \otimes I) J(\mathbf{z}_n)$

  - Sparse matrix methods become important for large systems!

\emoji{green-square} Scales well in model space

\emoji{red-square} Slower (but fast enough?)

# Feature 3: differentiability

![Derivatives of matter power spectrum wrt. parameters](media/power_spectrum_diff.png){height=80%}

# Method 1: manual differentiation

```
L(w₁,w₂) = w₂*log(w₁) + √(w₂*log(w₁))

∂L₁(w₁,w₂) = w₂/w₁ + (w₂/w₁) / (2*√(w₂*log(w₁)))
∂L₂(w₁,w₂) = log(w₁) + log(w₁) / (2*√(w₂))
∂L(w₁,w₂) = [∂L₁(w₁,w₂), ∂L₂(w₁,w₂)]
```

\emoji{green-square} Exact

\emoji{red-square} Tedious, scales horribly, only explicit expressions

# Method 2: symbolic differentiation

```
L(w₁,w₂) = w₂*log(w₁) + √(w₂*log(w₁))

using Symbolics # CAS package
function ∂L(w₁,w₂)
   @variables W₁ W₂
   ∂L₁ = expand_derivatives(Differential(W₁)(L(W₁,W₂)))
   ∂L₂ = expand_derivatives(Differential(W₂)(L(W₁,W₂)))
   return substitute([∂L₁, ∂L₂], Dict(W₁ => w₁, W₂ => w₂))
end
```

\emoji{green-square} Exact

\emoji{red-square} Scales horribly, only explicit expressions


# Method 3: numerical differentiation (finite differences)

```
L(w₁,w₂) = w₂*log(w₁) + √(w₂*log(w₁))

∂L₁(w₁,w₂; ϵ=1e-5) = (L(w₁+ϵ/2,w₂) - L(w₁-ϵ/2,w₂)) / ϵ
∂L₂(w₁,w₂; ϵ=1e-5) = (L(w₁,w₂+ϵ/2) - L(w₁,w₂-ϵ/2)) / ϵ
∂L(w₁,w₂; ϵ=1e-5) = [∂L₁(w₁,w₂;ϵ), ∂L₂(w₁,w₂;ϵ)]
```

\emoji{green-square} Simple, reuses original code

\emoji{red-square} Depends on $\epsilon$, unstable, time \emph{or} space complexity $\mathcal{O(n)}$ for $\nabla f(x_1,\ldots,x_n)$


# Method 4: automatic differentiation ("autodiff")

```
L(w₁,w₂) = w₂*log(w₁) + √(w₂*log(w₁))

using ForwardDiff
∂L(w₁,w₂) = ForwardDiff.gradient(L, [w₁,w₂]) # ???
```

\emoji{green-square} Simple, reuses code, stable, no $\epsilon$, handles any function

\emoji{red-square} Runtime $\mathcal{O}(n)$, needs full source code

# Automatic differentiation: motivation

Gradients are important in

- Optimization (e.g. parameter inference; Metropolis-Hastings $\rightarrow$ Hamiltonian Monte Carlo)

- Machine learning tools (e.g. emulators)

- Implicit ODE methods (e.g. stiff Boltzmann equations)

- Any calculation with derivatives as input (e.g. Fisher forecasts)

- Pedagogical aspect: understand how output changes with input

# Automatic differentiation: intuition

Any computer program is a (big) composite function $f(x) = f_n(f_{n-1}(\cdots f_2(f_1(x))))$ of elementary operations $f_i(x)$.

Let the compiler transform code for $f(x)$ into code for $f\prime(x)$.

Chain rule on $f(x) = f_3(f_2(f_1(x)))$ can be traversed in two ways:

- $\displaystyle \frac{\partial f}{\partial x} = \left( \frac{\partial f_3}{\partial f_2} \left( \frac{\partial f_2}{\partial f_1} \left( \frac{\partial f_1}{\partial x} \right) \right) \right) \quad\rightarrow\quad \text{forward-mode}$

- $\displaystyle \frac{\partial f}{\partial x} = \left( \left( \left( \frac{\partial f_3}{\partial f_2} \right) \frac{\partial f_2}{\partial f_1} \right) \frac{\partial f_1}{\partial x} \right) \quad\rightarrow\quad \text{reverse-mode}$

# Automatic differentiation: forward-mode

# Automatic differentiation: reverse-mode


# Differentiable programming

AD needs source access: spawned new "programming paradigm"

Reverse-mode AD largely turned into a compiler and language design problem

# Differentiation methods

[![](media/differentiation.png)](https://fmin.xyz/docs/methods/Autograd.html)

# Cheap gradient principle

$\mathcal{O}(\text{simulation}) \sim \mathcal{O}(\text{optimization})$
*independent of $n$!*

# Why start from scratch in Julia?

- Differentiable programs need entire source code
- Modular (SymBoltz) vs. monolithic (CLASS)
- Simplify construction of alternative models
- Want to contribute to cosmology environment in Julia
- Julia has excellent ODE, modeling and AD libraries: best tool for the job?
- Perturbation theory is computer food
- "Bottom-up" alternative to "top-down" generalized EFT approaches
- Focus on high-level equations instead of low-level perturbations


# DifferentialEquations.jl: "#1 differential equations library"

![](media/de_comparison.pdf)

# ModelingToolkit.jl: symbolic-numeric modeling language

- Component-based model building
- Automatic symbolic index management
- Generate system $\symbfit{y}^\prime(t) = \symbfit{f}(\symbfit{y}(t))$ or $\symbfit{F}(\symbfit{y}^\prime(t), \symbfit{y}(t), t) = \symbf{0}$
- Generate analytical Jacobian matrix

# ModelingToolkit.jl: my contributions

- Still maturing

![](media/issues.png){width=49%}
![](media/prs.png){width=49%}

# ForwardDiff.jl: forward mode automatic differentiation

# New features

\begin{center}
% https://stackoverflow.com/a/56695563, https://tex.stackexchange.com/a/56591
\begin{tikzpicture}[
   pillar/.style={circle, draw, minimum size=2.5cm, align=center, font=\scriptsize},
   joiner/.style={midway, sloped, black, font=\scriptsize},
   arrow/.style={-latex, thick, red},
]
\node[pillar, fill=gray] (sym) at (0:0cm) {Easy symbolic\\extensibility};
\node[pillar, fill=gray] (app) at (240:6cm) {Approximation-\\free};
\node[pillar, fill=gray] (dif) at (300:6cm) {Automatic\\differentiation};
\end{tikzpicture}
\end{center}

# New features are self-reinforcing

\begin{center}
% https://stackoverflow.com/a/56695563, https://tex.stackexchange.com/a/56591
\begin{tikzpicture}[
   pillar/.style={circle, draw, minimum size=2.5cm, align=center, font=\scriptsize},
   joiner/.style={midway, sloped, black, font=\scriptsize},
   arrow/.style={-latex, thick, red},
]
\node[pillar, fill=gray] (sym) at (0:0cm) {Easy symbolic\\extensibility};
\node[pillar, fill=gray] (app) at (240:6cm) {Approximation-\\free};
\node[pillar, fill=gray] (dif) at (300:6cm) {Automatic\\differentiation};
\draw[arrow] (sym.305) -- node[joiner, above] {$J$ sparsity detection} (dif.115);
\draw[arrow] (sym.245) -- node[joiner, below] {Generate equations?} (app.55);
\draw[arrow] (dif.185) -- node[joiner, below] {Generate accurate $J$} (app.355);
\draw[arrow] (dif.125) -- node[joiner, below] {Model-agnostic shooting} (sym.295);
\draw[arrow] (app.65) -- node[joiner, above] {Easy to specify model} (sym.235);
\draw[arrow] (app.5) -- node[joiner, above] {No discontinuities} (dif.175);
\end{tikzpicture}
\end{center}

# Future plans:

- Symbolic tensor library for model building
   - Compute $\Gamma^\alpha_{\beta\gamma}$, $R^\alpha_{\beta\gamma\delta}$, $R_{\alpha\beta}$, $G_{\alpha\beta}$, ... easily
   - Specify full equations of motion; then expand perturbatively
- Handle initial conditions symbolically
- Unit checking
- CMB power spectrum

# Miscellaneous

- [Standard reference for cosmological perturbations](https://arxiv.org/abs/astro-ph/9506072)

- [Why use ModelingToolkit.jl](https://discourse.julialang.org/t/what-are-the-advantages-of-using-modelingtoolkit-to-make-and-solve-odeproblems/103563)

- [Causal vs Acausal Modeling By Example: Why Julia ModelingToolkit.jl Scales](https://www.youtube.com/watch?v=ZYkojUozeC4)

- How to do equation-specific tweaks?

- *Why* re-implement everything from scratch, instead of building upon existing works?

- 2nd order code: [`song`](https://github.com/coccoinomane/song)

- CAMB has a "[symbolic toolkit notebook](https://camb.readthedocs.io/en/latest/ScalEqs.html)"

- [Lesgourges' CLASS presentations](https://lesgourg.github.io/courses.html): [content](https://lesgourg.github.io/class-tour/Padova/CLASS_Padova_Content.pdf), [coding](https://lesgourg.github.io/class-tour/Padova/CLASS_Padova_Coding.pdf), [exercises](https://lesgourg.github.io/class-tour/Padova/CLASS_Padova_Exercises.pdf)

- [Implicit ODE solver source](https://www.epfl.ch/labs/anchp/wp-content/uploads/2018/05/part2-1.pdf)

---

- Chris Rackauckas' StackExchange answers [BDF vs implicit](https://scicomp.stackexchange.com/questions/27178/bdf-vs-implicit-runge-kutta-time-stepping) and [stiff solver time complexity](https://scicomp.stackexchange.com/questions/28617/comparing-algorithmic-complexity-ode-solvers-big-o?noredirect=1&lq=1)

- Search for [Creative Commons media](https://openverse.org/)

- [Hans' Boltzmann solver output](https://cmb.wintherscoming.no/calc.php)

- [Øyvind's presentation](https://www.mn.uio.no/fysikk/english/research/groups/theory/theory-seminars/2024_christiansen.pdf)

- [CLASS: "TCA could even be removed](https://lesgourg.github.io/class-tour/Padova/CLASS_Padova_Content.pdf#page=47)

# Miscellaneous

Definitions
$$
\begin{aligned}
u^\mu &= \ldots, \\
T^{\mu\nu} &= (\rho + P) u^\mu u^\nu + P g^{\mu\nu}, \\
G_{\mu\nu} &= \ldots,
\end{aligned}
$$

Equations of motion
$$
\begin{aligned}
G_{\mu\nu} &= 8 \pi T_{\mu\nu}, \\
\nabla_\mu (n u^\mu) &= 0, \\
\nabla_\mu T^{\mu\nu} &= 0, \\
\ldots
\end{aligned}
$$

Relativistic fluid mechanics on a manifold?

Handle interactions on symbolic level;
then generate efficient numerical functions.
