# Getting started

[Install SymBoltz.jl](@ref "Installation") and load it with
```@example getting_started
using SymBoltz
```

## 1. Create the model

The first step is to define our cosmological model.
This is a *symbolic* representation of the variables and equations that describe the various components in the universe, such as the theory of gravity (like general relativity) and particle species (like the cosmological constant, cold dark matter, photons and baryons).

To get started, we will simply load the standard ΛCDM model:
```@example getting_started
M = ΛCDM()
hierarchy(M; describe = true)
```
As shown, the model is structured as a hierarchy of the physical components.
Each of these components contains a block of self-contained variables and equations that are independent from the other components, making up interchangeable modules of the entire Einstein-Boltzmann system.
A full model joins several such incomplete blocks into a complete set of equations for the entire Einstein-Boltzmann system.

The hierarchical structure can be inspected interactively (with TAB-completion) in the Julia REPL by evaluating `M`, `M.G`, `M.G.ρ` and so on.
For example, to see all equations for the theory of gravity:
```@example getting_started
equations(M.G)
```

## 2. Solve the problem

Next, we create a *numerical* representation of the cosmological problem we want to solve.
This splits the full symbolic model into computational stages (like the background, thermodynamics and perturbations), assigns input values to parameters and defines any parameters that are solved for with the shooting method by matching conditions at the final time.
```@example getting_started
using Unitful, UnitfulAstro # for interfacing without internal code units
pars = Dict(
    M.γ.T₀ => 2.7,
    M.c.Ω₀ => 0.27,
    M.b.Ω₀ => 0.05,
    M.ν.Neff => 3.0,
    M.g.h => 0.7,
    M.b.rec.YHe => 0.25,
    M.h.m_eV => 0.06,
    M.I.ln_As1e10 => 3.0,
    M.I.ns => 0.95
)
ks = 10 .^ range(-5, 1, length=500) / u"Mpc"
prob = CosmologyProblem(M, pars)
```
Finally, we can simply solve the problem:

```@example getting_started
sol = solve(prob, ks) # or just solve(prob) to solve only the background
```

!!! warning
    Constructing the `CosmologyProblem` is an expensive operation!
    It compiles symbolics to numerics, and should not be used to update parameter values.
    To just update parameter values, use the [parameter updater](@ref "Updating the parameters") function.

!!! tip
    For maximum performance, use the optimal BLAS backend for your platform, such as [MKL](https://github.com/JuliaLinearAlgebra/MKL.jl) for Intel processors.
    See the [performance page](@ref "Performance and benchmarks") for more information.
    For example:
    ```julia
    using MKL
    ```

## 3. Use the solution

You are now free to [do whatever you want with the solution object](@ref "Solving models").
For example, to get the time points used by the solver and corresponding values of the scale factor $a(τ)$:
```@example getting_started
τs = sol[M.τ]
as = sol(M.g.a, τs)
```
Similarly, to get $\Phi(k,τ)$ for the 500 wavenumbers we solved for at the same times:
```@example getting_started
Φs = sol(M.g.Φ, τs, ks)
```

You could plot this with `using Plots; plot(log10.(as), transpose(Φs))`, but this is more convenient with the included plot recipe:
```@example getting_started
using Plots
ks_plot = [1e-3, 1e-2, 1e-1, 1e-0] / u"Mpc"
plot(sol, log10(M.g.a), M.g.Φ, ks_plot) # lg(a) vs. Φ for 4 wavenumbers
```

We can also calculate the matter power spectrum:
```@example getting_started
Ps = spectrum_matter(sol, ks)
plot(log10.(ks/u"1/Mpc"), log10.(Ps/u"Mpc^3"); xlabel = "lg(k/Mpc⁻¹)", ylabel = "lg(P/Mpc³)", label = nothing)
```
Similarly, we can calculate the angular CMB (TT) power spectrum:
```@example getting_started
ls = 10:10:1000
Dls = spectrum_cmb(:TT, prob, ls; normalization = :Dl, unit = u"μK")
plot(ls, Dls; xlabel = "l", ylabel = "l (l+1) Cₗ / 2π", label = nothing)
```

And here is a condensed plot with several quantities:
```@example getting_started
p = plot(layout=(3, 3), size=(900, 700), tickfontsize=6, labelfontsize=6, legendfontsize=5)
plot!(p[1], sol, log10(M.g.a), [M.b.ρ, M.c.ρ, M.γ.ρ, M.ν.ρ, M.h.ρ, M.Λ.ρ] ./ M.G.ρ)
plot!(p[2], sol, log10(M.g.a), [M.b.w, M.c.w, M.γ.w, M.ν.w, M.h.w, M.Λ.w])
plot!(p[3], sol, log10(M.g.a), log10(M.g.E))
plot!(p[4], sol, log10(M.g.a), [M.b.rec.XHe⁺⁺, M.b.rec.XHe⁺, M.b.rec.XH⁺, M.b.rec.Xe])
plot!(p[5], sol, log10(M.g.a), log10.([M.b.rec.Tγ, M.b.rec.Tb] ./ M.γ.T₀))
plot!(p[6], sol, log10(M.g.a), log10(abs(M.b.rec.κ)))
plot!(p[7], sol, log10(M.g.a), [M.g.Φ, M.g.Ψ], ks_plot)
plot!(p[8], sol, log10(M.g.a), log10.(abs.([M.b.δ, M.c.δ, M.γ.δ, M.ν.δ, M.h.δ])), ks_plot; klabel = false)
plot!(p[9], sol, log10(M.g.a), log10.(abs.([M.b.θ, M.c.θ, M.γ.θ, M.ν.θ, M.h.θ])), ks_plot; klabel = false)
```
