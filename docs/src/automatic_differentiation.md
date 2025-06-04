# Using automatic differentiation

This tutorial shows how to compute the power spectrum $P(k; \theta)$
and its (logarithmic) derivatives
```math
\frac{\partial \lg P}{\partial \lg \theta_i}
```
using automatic differentiation with [ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl).
This technique [can also differentiate any other quantity](@ref "General approach").

## 1. Wrap the evaluation

First, we must decide which parameters $\theta$ the power spectrum $P(k; \theta)$ should be considered a function of.
To do so, let us write a small wrapper function that calculates the power spectrum as a function of the parameters $(T_{\gamma 0}, \Omega_{c0}, \Omega_{b0}, N_\textrm{eff}, h, Y_p)$, following the [Getting started tutorial](@ref "Getting started"):
```@example ad
using SymBoltz
M = ΛCDM(K = nothing)
pars = [M.γ.T₀, M.c.Ω₀, M.b.Ω₀, M.ν.Neff, M.g.h, M.b.rec.Yp, M.h.m_eV, M.I.As, M.I.ns]
prob0 = CosmologyProblem(M, Dict(pars .=> NaN))

probgen = SymBoltz.parameter_updater(prob0, pars)
P(k, θ) = spectrum_matter(probgen(θ), k; verbose = true, ptopts = (reltol = 1e-3,))
```
It is now easy to evaluate the power spectrum:
```@example ad
using Unitful, UnitfulAstro
θ = [2.7, 0.27, 0.05, 3.0, 0.7, 0.25, 0.06, 2e-9, 0.95]
ks = 10 .^ range(-3, 0, length=100) / u"Mpc"
Ps = P(ks, θ)
```
This can be plotted with
```@example ad
using Plots
plot(log10.(ks/u"1/Mpc"), log10.(Ps/u"Mpc^3"); xlabel = "lg(k/Mpc⁻¹)", ylabel = "lg(P/Mpc³)", label = nothing)
```

## 2. Calculate the derivatives

To get $\partial \lg P / \partial \lg \theta$, we can simply pass the wrapper function `P(k, θ)` through [`ForwardDiff.jacobian`](https://juliadiff.org/ForwardDiff.jl/stable/user/api/#ForwardDiff.jacobian):
```@example ad
using ForwardDiff
lgP(lgθ) = log10.(P(ks, 10 .^ lgθ) / u"Mpc^3") # in log-space
dlgP_dlgθs = ForwardDiff.jacobian(lgP, log10.(θ))
```
The matrix element `dlgP_dlgθs[i, j]` now contains $\partial \lg P(k_i) / \partial \lg \theta_j$.
We can plot them all at once:
```@example ad
plot(
    log10.(ks/u"1/Mpc"), dlgP_dlgθs;
    xlabel = "lg(k/Mpc⁻¹)", ylabel = "∂ lg(P) / ∂ lg(θᵢ)",
    labels = "θᵢ=" .* ["Tγ0" "Ωc0" "Ωb0" "Neff" "h" "Yp" "mh" "As" "ns"]
)
```

## Get values and derivatives together

The above example showed how to calculate the power spectrum *values* and their *derivatives* through *two separate calls*.
If you need both, it is faster to calculate them simultaneously with the package [`DiffResults.jl`](https://juliadiff.org/DiffResults.jl/stable/):
```@example ad
using DiffResults

# Following DiffResults documentation:
Pres = DiffResults.JacobianResult(ks/u"1/Mpc", θ) # allocate buffer for values+derivatives for a function with θ-sized input and ks-sized output
Pres = ForwardDiff.jacobian!(Pres, lgP, log10.(θ)) # evaluate values+derivatives of lgP(log10.(θ)) and store the results in Pres
lgPs = DiffResults.value(Pres) # extract values
dlgP_dlgθs = DiffResults.jacobian(Pres) # extract derivatives

p1 = plot(
   log10.(ks/u"1/Mpc"), lgPs;
   ylabel = "lg(P/Mpc³)", label = nothing
)
p2 = plot(
   log10.(ks/u"1/Mpc"), dlgP_dlgθs;
   xlabel = "lg(k/Mpc⁻¹)", ylabel = "∂ lg(P) / ∂ lg(θᵢ)",
   labels = "θᵢ=" .* ["Tγ0" "Ωc0" "Ωb0" "Neff" "h" "Yp" "mh" "As" "ns" "ΩΛ0"]
)
plot(p1, p2, layout=(2, 1), size = (600, 600))
```

## General approach

The technique shown here can be used to calculate the derivative of any SymBoltz.jl output quantity:

1. Write a wrapper function `output(input)` that calculates the desired output quantities from the desired input quantities.
2. Use [`ForwardDiff.derivative(output, input)`](https://juliadiff.org/ForwardDiff.jl/stable/user/api/#ForwardDiff.derivative) (scalar-to-scalar), [`ForwardDiff.gradient(output, input)`](https://juliadiff.org/ForwardDiff.jl/stable/user/api/#ForwardDiff.gradient) (vector-to-scalar) or [`ForwardDiff.jacobian(output, input)`](https://juliadiff.org/ForwardDiff.jl/stable/user/api/#ForwardDiff.jacobian) (vector-to-vector) to evaluate the derivative of `output` at the values `input`.
   Or use the similar functions in [`DiffResults`](https://juliadiff.org/DiffResults.jl/stable/) to calculate the value *and* derivatives simultaneously.
