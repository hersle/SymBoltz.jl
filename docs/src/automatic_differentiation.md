# Automatic differentiation

This tutorial shows how to compute the (cold dark matter) power spectrum $P(k; \theta)$
and its (logarithmic) derivatives
```math
\frac{\partial \lg P}{\partial \lg \theta_i}
```
using automatic differentiation with [ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl).
The technique [can also be used to differentiate any other quantity](@ref "General approach").

## 1. Wrapping the evaluation

*We* must first decide which parameters $\theta$ the power spectrum $P(k; \theta)$ will be considered a function of.
To do so, let us write a small wrapper function that calculates the power spectrum as a function of the parameters $(\Omega_{\gamma 0}, \Omega_{c0}, \Omega_{b0}, N_\textrm{eff}, h, Y_p$, following the [Getting started tutorial](@ref "Getting started"):
```@example 1
using ModelingToolkit
import SymBoltz

@named M = SymBoltz.ΛCDM()
prob = SymBoltz.CosmologyProblem(M)

# define ordering and values of parameters
θ_syms = [M.γ.Ω0, M.c.Ω0, M.b.Ω0, M.ν.Neff, M.g.h, M.b.rec.Yp] # symbolic indices
θ_strs = ["Ωγ0" "Ωc0" "Ωb0" "Neff" "h" "Yp"] # plot labels
θ0 = [5e-5, 0.27, 0.05, 3.0, 0.7, 0.25] # numerical values

P(k, θ) = SymBoltz.P(SymBoltz.solve(prob, θ_syms .=> θ, k), M.c, k)
```
It is now easy to evaluate the power spectrum:
```@example 1
ks = 10 .^ range(-3, 0, length=100) / SymBoltz.k0 # TODO: improve on units
Ps = P(ks, θ0)
```
This can be plotted with
```@example 1
using Plots
plot(log10.(ks), log10.(Ps); xlabel = "lg(k)", ylabel = "lg(P)", label = nothing)
```

## 2. Calculate the derivatives

To get $\partial \lg P / \partial \lg \theta$, we can simply pass the wrapper function `P(k, θ)` through [`ForwardDiff.jacobian`](https://juliadiff.org/ForwardDiff.jl/stable/user/api/#ForwardDiff.jacobian):
```@example 1
using ForwardDiff
lgP(lgθ) = log10.(P(ks, 10 .^ lgθ)) # ForwardDiff.jacobian needs an array -> array function
dlgP_dlgθs = ForwardDiff.jacobian(lgP, log10.(θ0))
```
The matrix element `dlgP_dlgθs[i, j]` now contains $\partial \lg P(k_i) / \partial \lg \theta_j$.
We can plot them all at once:
```@example 1
plot(log10.(ks), dlgP_dlgθs; xlabel = "lg(k)", ylabel = "∂ lg(P) / ∂ lg(θᵢ)", labels = "θᵢ=" .* θ_strs)
```

## General approach

The technique shown here can be used to calculate the derivative of any SymBoltz.jl output quantity:

1. Write a wrapper function `output(input)` that calculates the desired output quantities from the desired input quantities.
2. Use [`ForwardDiff.derivative(output, input)`](https://juliadiff.org/ForwardDiff.jl/stable/user/api/#ForwardDiff.derivative) (scalar-to-scalar), [`ForwardDiff.gradient(output, input)`](https://juliadiff.org/ForwardDiff.jl/stable/user/api/#ForwardDiff.gradient) (vector-to-scalar) or [`ForwardDiff.jacobian(output, input)`](https://juliadiff.org/ForwardDiff.jl/stable/user/api/#ForwardDiff.jacobian) (vector-to-vector) to evaluate the derivative of `output` at the values `input`.