# Solution interface

!!! note
    This page is a work in progress.

```@docs
solve(M::CosmologyModel, pars)
solve(M::CosmologyModel, pars, ks::AbstractArray)
```

For the `ΛCDM()` model, some good solvers from `OrdinaryDiffEq` are

- `TRBDF2`: TODO note
- `KenCarp4`, `KenCarp47`: TODO note
- `Kvaerno5`: TODO note
- `Rodas5P`: TODO note
