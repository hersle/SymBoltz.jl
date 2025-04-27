# Code conventions

The conventions in the code are inspired by the [seminal paper by Ma and Bertschinger (1995)](https://arxiv.org/abs/astro-ph/9506072).

## Sign conventions

- Metric tensor: $g = a^2 \mathrm{diag}(-1-2ϵΨ, +1-2ϵΦ, +1-2ϵΦ, +1-2ϵΦ)$ (i.e. $Ψ = Φ$ under no stress and $\mathrm{sign}(g) = \mathrm{diag}(-,+,+,+)$).

## Variable names

- `τ` ($t$) is *conformal time*.
- `κ` ($κ$) is *optical depth*.

## Dimensions

- In background and perturbation equations, all numbers are dimensionless, defined as their dimensionful counterparts in units of $H_0$, $G$ and $c$.
  For example, `ρ` is mass density in units of $H_0^2/G$, i.e. the dimensionless number $\rho_\text{code} = \rho_\text{phys} / (H_0^2/G)$.
- In recombination equations, all quantities are in SI units. (TODO: make everything consistent in Mpc units?)
