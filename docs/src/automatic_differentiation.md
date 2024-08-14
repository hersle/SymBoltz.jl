# Automatic differentiation

This tutorial shows how to compute the (cold dark matter) power spectrum $P(k; \theta)$
and its (logarithmic) derivatives
```math
\frac{\partial \lg P}{\partial \lg \theta_i}
```
using automatic differentiation with [ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl).
This technique [can also differentiate any other quantity](@ref "General approach").

## 1. Wrapping the evaluation

*We* must first decide which parameters $\theta$ the power spectrum $P(k; \theta)$ will be considered a function of.
To do so, let us write a small wrapper function that calculates the power spectrum as a function of the parameters $(\Omega_{\gamma 0}, \Omega_{c0}, \Omega_{b0}, N_\textrm{eff}, h, Y_p)$, following the [Getting started tutorial](@ref "Getting started"):
```@example 1
using SymBoltz

M = ΛCDM()
prob = CosmologyProblem(M)

# define ordering and values of parameters
θ_syms = [M.γ.Ω0, M.c.Ω0, M.b.Ω0, M.ν.Neff, M.g.h, M.b.rec.Yp] # symbolic indices
θ_strs = ["Ωγ0" "Ωc0" "Ωb0" "Neff" "h" "Yp"] # plot labels
θ0 = [5e-5, 0.27, 0.05, 3.0, 0.7, 0.25] # numerical values

P(k, θ) = power_spectrum(solve(prob, θ_syms .=> θ, k), M.c, k)
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

## Get values and derivatives together

The above example showed how to calculate the power spectrum *values* and their *derivatives* through *two separate calls*.
If you need both, it is faster to calculate them simultaneously with the package [`DiffResults.jl`](https://juliadiff.org/DiffResults.jl/stable/):
```@example 1
using DiffResults

Pres = DiffResults.JacobianResult(ks, θ0) # allocate buffer for value+derivatives for a function with θ0-sized input and ks-sized output
Pres = ForwardDiff.jacobian!(Pres, lgP, log10.(θ0)) # evaluate value+derivatives of lgP(log10.(θ0)) and store the results in Pres
lgPs = DiffResults.value(Pres) # extract value
dlgP_dlgθs = DiffResults.jacobian(Pres) # extract derivatives

p1 = plot(log10.(ks), lgPs; ylabel = "lg(P)", label = nothing)
p2 = plot(log10.(ks), dlgP_dlgθs; xlabel = "lg(k)", ylabel = "∂ lg(P) / ∂ lg(θᵢ)", labels = "θᵢ=" .* θ_strs)
plot(p1, p2, layout=(2, 1), size = (600, 600))
```

## General approach

The technique shown here can be used to calculate the derivative of any SymBoltz.jl output quantity:

1. Write a wrapper function `output(input)` that calculates the desired output quantities from the desired input quantities.
2. Use [`ForwardDiff.derivative(output, input)`](https://juliadiff.org/ForwardDiff.jl/stable/user/api/#ForwardDiff.derivative) (scalar-to-scalar), [`ForwardDiff.gradient(output, input)`](https://juliadiff.org/ForwardDiff.jl/stable/user/api/#ForwardDiff.gradient) (vector-to-scalar) or [`ForwardDiff.jacobian(output, input)`](https://juliadiff.org/ForwardDiff.jl/stable/user/api/#ForwardDiff.jacobian) (vector-to-vector) to evaluate the derivative of `output` at the values `input`.
   Or use the similar functions in [`DiffResults`](https://juliadiff.org/DiffResults.jl/stable/) to calculate the value *and* derivatives simultaneously.
