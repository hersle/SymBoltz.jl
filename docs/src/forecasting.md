# Fisher forecasting

Fisher forecasting is an important tool for estimating what constraints on parameters that observations with some given uncertainties can place.
This example shows how SymBoltz can be used to perform a Fisher forecast on a CMB (TT) survey limited only by cosmic variance.

First, create a base ΛCDM cosmological model and problem:
```@example forecast
# inspiration: e.g. https://github.com/xzackli/fishchips-public/blob/master/notebooks/Introduction%20to%20Fisher%20Forecasting.ipynb # hide
# TODO: start by getting equal ad/fd results with these, then include more parameters # hide
using SymBoltz, Plots
M = ΛCDM(K = nothing) # flat
pars = Dict(
    M.g.h => 0.70,
    M.c.Ω₀ => 0.27,
    M.b.Ω₀ => 0.05,
    M.γ.T₀ => 2.7,
    M.ν.Neff => 3.0,
    M.b.YHe => 0.25,
    M.h.m_eV => 0.02,
    M.I.ln_As1e10 => 3.0,
    M.I.ns => 0.96
)
# TODO: more parameters, try one-by-one: Neff is a bit iffy # hide
pars_varying = [M.g.h, M.c.Ω₀, M.b.Ω₀, M.b.YHe, M.I.ln_As1e10, M.I.ns] # parameters to be varied; others are fixed
prob0 = CosmologyProblem(M, merge(pars, Dict(pars_varying .=> NaN))) # set varying to NaN
```

Next, create a function for computing $Cₗ$ of the CMB TT power spectrum.
Since $Cₗ$ is an expensive but smooth function of $l$, we make one function for it exactly on a coarse grid of $l$ and another for interpolating it to a finer grid:
```@example forecast
# TODO: ω0 better than Ω0? # hide
probgen = parameter_updater(prob0, pars_varying)
ls, ls′ = 40:1:1000, 40:20:1000
function Cl(θ; bgopts = (alg = SymBoltz.Rodas4P(), reltol = 1e-9, abstol = 1e-9), ptopts = (alg = SymBoltz.KenCarp4(), reltol = 1e-8, abstol = 1e-8))
    return spectrum_cmb(:TT, probgen(θ), ls, ls′; bgopts, ptopts)
end
```
We can now compute $Cₗ$ and the cosmic variance uncertainties
```math
σₗ = \sqrt{\frac{2}{2l+1}} Cₗ
```
and plot them with error bars:
```@example forecast
θ0 = [pars[par] for par in pars_varying]
Cls = Cl(θ0)
σs = @. √(2/(2ls+1)) * Cls # cosmic variance
plot(ls, Cls.*ls.*(ls.+1)/2π; ribbon = σs.*ls.*(ls.+1)/2π, xlabel = "l", ylabel = "l(l+1)Cₗ/2π", label = "Dₗ ± ΔDₗ")
```
The likelihood (logarithm) function (given model parameters $θ$ and measured $̄\bar{C}ₗ$) is
```math
\log L(θ) = -\frac12 \sum_i \left( \frac{Cₗ(θ)-C̄ₗ}{σᵢ} \right)².
```
To calculate parameter covariances, we first need the [Fisher information matrix](https://en.wikipedia.org/wiki/Fisher_information):
```math
Fᵢⱼ = -\frac{1}{2} ∑ₗ \left⟨ \frac{∂² \log L}{∂θᵢ ∂θⱼ} \right⟩ = ∑ₗ \frac{∂Cₗ}{∂θᵢ} \frac{1}{σₗ²} \frac{∂Cₗ}{∂θⱼ}.
```
Notice that the Fisher matrix is independent of the measured $C̄ₗ$, and depends only on the uncertainties and the *derivatives* of the $Cₗ$ with respect to the cosmological parameters.
We can compute the derivatives using automatic differentiation, and compare them to those found with finite differences:
```@example forecast
using ForwardDiff, FiniteDiff
dCl_dθ_ad = ForwardDiff.jacobian(Cl, θ0)
dCl_dθ_fd = FiniteDiff.finite_difference_jacobian(Cl, θ0, Val{:central}; relstep = 5e-3) # TODO: 4e-2 is good for m_eV

θnames = replace.(string.(pars_varying), "₊" => ".")
color = permutedims(eachindex(θnames))
hline([NaN NaN], color = :black, linestyle = [:solid :dash], xlabel = "l", label = ["AD" "FD"])
plot!(ls, dCl_dθ_ad ./ Cls; color, linestyle = :solid, label = "∂(Cₗ)/∂(" .* permutedims(θnames) .* ")")
plot!(ls, dCl_dθ_fd ./ Cls; color, linestyle = :dash, label = nothing)
```
We can now compute $F$ and the parameter covariance matrix
```math
C = F⁻¹,
```
which is simply the inverse of $F$:
```@example forecast
fisher_matrix(dCl_dθ) = [sum(dCl_dθ[il,i]*dCl_dθ[il,j]/σs[il]^2 for il in eachindex(ls)) for i in eachindex(θ0), j in eachindex(θ0)]
F_fd = fisher_matrix(dCl_dθ_fd)
F_ad = fisher_matrix(dCl_dθ_ad)
C_fd = inv(F_fd)
C_ad = inv(F_ad)
```
Finally, we can plot ellipses for the forecasted parameter constraints in the multi-dimensional parameter space:
```@example forecast
using Plots

function ellipse(C::Matrix, i, j, c = (0.0, 0.0); nstd = 1, N = 33)
    σᵢ², σⱼ², σᵢⱼ = C[i,i], C[j,j], C[i,j]
    θ = (atan(2σᵢⱼ, σᵢ²-σⱼ²)) / 2
    a = √((σᵢ²+σⱼ²)/2 + √((σᵢ²-σⱼ²)^2/4+σᵢⱼ^2))
    b = √(max(0.0, (σᵢ²+σⱼ²)/2 - √((σᵢ²-σⱼ²)^2/4+σᵢⱼ^2)))

    a *= nstd # TODO: correct?
    b *= nstd # TODO: correct?

    cx, cy = c
    ts = range(0, 2π, length=N)
    xs = cx .+ a*cos(θ)*cos.(ts) - b*sin(θ)*sin.(ts)
    ys = cy .+ a*sin(θ)*cos.(ts) + b*cos(θ)*sin.(ts)
    return xs, ys
end

function plot_ellipses!(p, C; label = nothing, kwargs...)
    for i in eachindex(IndexCartesian(), C)
        ix, iy = i[1], i[2]
        if iy == 1 || iy > size(p)[1] + 1 || ix > size(p)[2]
            continue # out of bounds; skip
        end
        subplot = p[iy-1, ix]
        if ix >= iy
            # upper triangular part
            _label = (iy-1, ix) == (1, size(p)[2]) ? label : nothing
            hline!(subplot, [NaN]; framestyle = :none, label = _label, legendfontsize = 10, kwargs...)
        else
            # lower triangular part
            μx = θ0[ix]
            μy = θ0[iy]
            xlabel = iy == length(θ0) ? θnames[ix] : ""
            ylabel = ix == 1 ? θnames[iy] : ""
            for nstd in 1:2
                xs, ys = ellipse(C, ix, iy, (μx, μy); nstd)
                plot!(subplot, xs, ys; xlabel, ylabel, label = nothing, kwargs...)
            end
        end
    end
    return p
end

p = plot(layout = (length(θ0)-1, length(θ0)-1), size = (1000, 1000), aspect = 1)
plot_ellipses!(p, C_ad; color = :blue, linestyle = :solid, linewidth = 2, label = "automatic differentiation")
plot_ellipses!(p, C_fd; color = :red, linestyle = :dash, linewidth = 2, label = "finite differences")
```

