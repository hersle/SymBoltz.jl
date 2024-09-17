---
title: |
    [![](components_baumann.pdf)](http://physics.bu.edu/~schmaltz/PY555/baumann_notes.pdf)  
    \   
    A modular and symbolic Einstein-Boltzmann solver
#titlegraphic: components_baumann.pdf
colorlinks: true
header-includes:
   - \usepackage{emoji}
   - \usepackage{tikz}
---

---



# Why?

- Perturbation theory is computer food

- Simplify construction of alternative models

- "Bottom-up" alternative to "top-down" generalized EFT approaches

- Contribute to modern cosmology environment in Julia

- Focus on high-level equations instead of low-level perturbations

# Main idea

Write a new Einstein-Boltzmann solver in Julia:

                   [CAMB](https://arxiv.org/abs/astro-ph/9911177) [CLASS](https://arxiv.org/abs/1104.2932) [Bolt.jl](https://github.com/xzackli/Bolt.jl) [DISCO-DJ](https://arxiv.org/abs/2311.03291) This
------------------ ---------------------------------------------- ---------------------------------------- --------------------------------------------- -------------------------------------------- --------------------
**Release**        2000                                           2011                                     ~2023                                         2023                                         Future
**Author**         A. Lewis                                       J. Lesg.                                 Z. Li                                         O. Hahn                                      Me
**Language**       Fortran                                        C                                        Julia                                         Python+Jax                                   Julia
**Readable**       \emoji{red-square}                             \emoji{green-square} (?)                 \emoji{green-square}                          \emoji{green-square}                         \emoji{green-square}?
**Flexible**       \emoji{red-square}                             \emoji{green-square}                     \emoji{green-square}                          \emoji{green-square}                         \emoji{green-square}?
**Differentiable** \emoji{red-square}                             \emoji{red-square}                       \emoji{green-square}                          \emoji{green-square}                         \emoji{green-square}?
**Symbolic**       \emoji{red-square}                             \emoji{red-square}                       \emoji{red-square}                            \emoji{red-square}                           \emoji{green-square}?
**Gauge**          Synchr.                                        Synchr.                                  Newt.                                         Synchr.                                      Any?
**Abstraction**    Low                                            Low                                      Low                                           Low                                          High?
**Extensibility**  \emoji{hot-face}                               \emoji{hot-face}                         \emoji{hot-face}                              \emoji{hot-face}                             \emoji{smiling-face-with-sunglasses}?
**Raw speed**      \emoji{dash}\emoji{dash}                       \emoji{dash}\emoji{dash}\emoji{dash}     \emoji{dash}\emoji{dash}                      \emoji{dash}\emoji{dash}                     \emoji{dash}\emoji{dash}?

# Why differentiable?

Compute change of observable $O$ with cosmological parameters $\theta$:
$$ \partial O / \partial \theta_i \qquad \text{(e.g. $\partial P(k, z \,|\, \Omega_{m0}, A_s, \ldots) / \partial \Omega_{m0}$)} $$

- Optimize $\mathcal{L}$-sampling in MCMC methods \newline
  (Metropolis-Hastings $\rightarrow$ Hamiltonian Monte Carlo)
- Understand how observables change with parameters


# Why symbolic?

1. Specify space-time gauge $g_{\mu\nu}$ and cosmological species $T^{(s)}_{\mu\nu}$
2. Specify high-level EOMs, e.g. Einstein & Boltzmann equations
   $$G_{\mu\nu} = 8 \pi T_{\mu \nu}, \quad \left( p^\mu \frac{\partial}{\partial x^\mu} - \Gamma^\mu_{\alpha \beta} p^\alpha p^\beta \frac{\partial}{\partial p^\mu} \right) f = C[f].$$
3. Generate and solve ODE systems order-by-order
   1. $\mathcal{O}(\epsilon^0)$: background ($k$-independent)
   2. $\mathcal{O}(\epsilon^1)$: linear perturbations (on a $k$-range)
   3. $\mathcal{O}(\epsilon^2)$: 2nd order perturbations (on a $k$-grid)

$\hrulefill$

- Conceptually clearer to add new species, modified gravity, ...
  - More efficient exploration of new models
- Systematic perturbations and togglable approximations
  - 2nd order, quasi-static limit, sub-horizon limit, ...
- Work in any gauge (by regenerating system with new metric)

# Components

::: incremental

1. [`Symbolics.jl`](https://docs.sciml.ai/Symbolics/): symbolic equation representation

   - Designed for "symbolics-numerics"
   
   - Compute $\Gamma^\alpha_{\beta\gamma}$, $R^\alpha_{\beta\gamma\delta}$, $R_{\alpha\beta}$, $G_{\alpha\beta}$, ... easily

   - Specify full equations of motion; then expand perturbatively

2. [`ModelingToolkit.jl`](https://docs.sciml.ai/ModelingToolkit/): composable model building

   - [*Compose*](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/rc_circuit/) models from components (`Species`, `Spacetime`, ...)

   - Generate system $\symbfit{y}^\prime(t) = \symbfit{f}(\symbfit{y}(t))$ or $\symbfit{F}(\symbfit{y}^\prime(t), \symbfit{y}(t), t) = \symbf{0}$

   - Generate analytical Jacobian matrix
   
3. [`DifferentialEquations.jl`](https://docs.sciml.ai/DiffEqDocs/): numerical (O)DE integration

   - ["#1 differential equations library"](https://www.stochasticlifestyle.com/wp-content/uploads/2019/11/de_solver_software_comparsion.pdf)

   - Fantastic selection of explicit and implicit methods

   - Produces auto-differentiable output $\partial y / \partial \symbf{\theta}$

:::

# Why is Julia suited to the task?

- `DifferentialEquations.jl` is the leading ODE library:
  ![](https://www.stochasticlifestyle.com/wp-content/uploads/2019/11/de_solver_software_comparsion.pdf){height=60%}
  - Implicitly handle stiff systems (fewer approximations)
  - Use DAE solvers on 2nd order systems
- `Symbolics.jl` is [designed for symbolic-numerics](https://symbolics.juliasymbolics.org/stable/comparison/)
- Marriage between autodiff+implicit+symbolics
- Expand to systematic LPT expansion for non-linear corrections, ...?

# How

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

# Unanswered questions and challenges

- How to handle initial conditions?

- How to handle thermodynamics / reactions?

  - Can this be done with [MomentClosure.jl](https://docs.sciml.ai/MomentClosure/dev/)?

- How easy to add the things I don't include?

- How to limit its scope?

- How to do equation-specific tweaks?

- *Why* re-implement everything from scratch, instead of building upon existing works?

- 2nd order code: [`song`](https://github.com/coccoinomane/song)

- CAMB has a "[symbolic toolkit notebook](https://camb.readthedocs.io/en/latest/ScalEqs.html)"

# References


- [Standard reference for cosmological perturbations](https://arxiv.org/abs/astro-ph/9506072)

- [Why use ModelingToolkit.jl](https://discourse.julialang.org/t/what-are-the-advantages-of-using-modelingtoolkit-to-make-and-solve-odeproblems/103563)

- [Causal vs Acausal Modeling By Example: Why Julia ModelingToolkit.jl Scales](https://www.youtube.com/watch?v=ZYkojUozeC4)

# Interconnected features

Triangle with differentiable, approximation-free and symbolical modularity
and how they enable each other.

# Finite differences

$f(\mathbf{x}) = f(x_1, \ldots, x_n)$
calculate $f$ and $\nabla f$ in $\mathcal{O}(2n+1)$ time

# Automatic differentiation

Cheap Gradient Principle:
$\mathcal{O}(\text{simulation}) \sim \mathcal{O}(\text{optimization})$
*independent of $n$!*

1. `f(x) = x^2`'

