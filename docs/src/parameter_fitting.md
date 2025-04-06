# Fitting parameters and forecasting

This tutorial shows how to perform Bayesian parameter inference on a cosmological model by fitting it to data, and how to forecast parameter constraints from a covariance matrix that describes uncertainties and correlations of (unknown) observed data.

## Observed supernova data

We use the data from [Betoule+ (2014)](https://arxiv.org/abs/1401.4064) of redshifts and luminosity distances of supernovae.
Specifically, we use the values and errors in tables F.1 and F.2
(available from the [data files `data/jla_mub.txt` and `data/jla_mub_covmatrix.dat`](http://supernovae.in2p3.fr/sdss_snls_jla/jla_likelihood_v6.tgz)):
```@example fit
# TODO: do things properly
# TODO: use full covariance
# TODO: w₀wₐ
using LinearAlgebra # TODO: remove after using full covar?
data = [ # https://cmb.wintherscoming.no/data/supernovadata.txt
#   z        dL/Gpc    ΔdL/Gpc
    0.010    0.0390    0.0026
    0.012    0.0597    0.0046 
    0.014    0.0587    0.0021 
    0.016    0.0666    0.0022 
    0.019    0.0829    0.0033 
    0.023    0.0972    0.0025 
    0.026    0.1123    0.0032 
    0.031    0.1412    0.0037 
    0.037    0.1637    0.0043 
    0.043    0.1936    0.0067 
    0.051    0.2139    0.0092 
    0.060    0.2701    0.0077 
    0.070    0.3062    0.0093 
    0.082    0.3902    0.0098 
    0.097    0.4473    0.0123 
    0.114    0.5279    0.0091 
    0.134    0.6510    0.0116 
    0.158    0.7384    0.0118 
    0.186    0.9087    0.0134 
    0.218    1.0747    0.0163 
    0.257    1.2971    0.0189 
    0.302    1.5172    0.0274 
    0.355    1.9243    0.0297 
    0.418    2.2813    0.0436 
    0.491    2.7943    0.0507 
    0.578    3.3373    0.0552 
    0.679    4.0789    0.1179 
    0.799    5.0213    0.1262 
    0.940    6.2302    0.1917 
    1.105    7.9947    0.5692 
    1.300    9.2121    0.5874 
]
reverse!(data; dims=1) # make redshift increasing
data = (zs = data[:,1], dLs = data[:,2], ΔdLs = data[:,3], C = Diagonal(data[:,3] .^ 2))

using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "lg(1+z)", ylabel = "dₗ/z/Gpc")
scatter!(ax, log10.(data.zs.+1), data.dLs ./ data.zs; label = "observations")
errorbars!(ax, log10.(data.zs.+1), data.dLs ./ data.zs, data.ΔdLs ./ data.zs) #, xlabel = "lg(1+z)", ylabel = "dₗ/z/Gpc", label = "observations")
fig
```

## Predicting luminosity distances

To predict luminosity distances
```math
d_L = \frac{r}{a} = \chi \, \mathrm{sinc} (\sqrt{k} \chi),
\qquad \text{where} \qquad
\chi = c \, (t_0 - t)
```
theoretically, we solve the standard ΛCDM model:
```@example fit
using SymBoltz, ModelingToolkit

g = SymBoltz.metric()
K = SymBoltz.curvature(g)
M = RMΛ(K = K)
M = change_independent_variable(M, M.g.a; add_old_diff = true)
pars = Dict(M.r.Ω₀ => 9.3e-5, M.m.Ω₀ => 0.26, M.Λ.Ω₀ => 0.66, M.g.h => 0.7, M.t => 0.0, M.r.T₀ => NaN)
pars[M.K.Ω₀] = 1 - pars[M.r.Ω₀] - pars[M.m.Ω₀] - pars[M.Λ.Ω₀]
prob = CosmologyProblem(M, pars; pt = false, ivspan = (1e-8, 1e0))

# TODO: use luminosity_distance
using SymbolicIndexingInterface # TODO: move this stuff into package
geta = getu(prob.bg, M.g.a)
gett = getu(prob.bg, M.t)
getΩk0 = getp(prob.bg, M.K.Ω₀)
geth = getp(prob.bg, M.g.h)
function dL(sol::CosmologySolution, geta, gett, getΩk0, geth)
    a, t, Ωk0, h = geta(sol.bg), gett(sol.bg), getΩk0(sol.bg), geth(sol.bg)
    t0 = t[end]
    χ = t0 .- t
    r = @. real(sinc(√(-Ωk0+0im)*χ/π) * χ) # Julia's sinc(x) = sin(π*x) / (π*x)
    H0 = SymBoltz.H100 * h
    dLs = @. r / a * SymBoltz.c / H0 # to meters
    return dLs / SymBoltz.Gpc # to Gpc
end

# Show example predictions
as = 1 ./ (data.zs .+ 1)
bgopts = (alg = SymBoltz.Tsit5(), reltol = 1e-6, maxiters = 1e3)
sol = solve(prob; bgopts, saveat = as, save_end = true, term = nothing)
dLs = dL(sol, geta, gett, getΩk0, geth)
scatter!(ax, sol[log10(M.g.z+1)], dLs ./ sol[M.g.z]; label = "predictions")
axislegend(ax, position = :rb)
fig
```

## Bayesian inference

To perform bayesian inference, we define a probabilistic model in [Turing.jl](https://turinglang.org/):
```@example fit
using Turing

@model function supernova(zs, dLs, ΔdLs, probgen, geta, gett, getΩk0, geth; bgopts = bgopts, Ωr0 = 9.3e-5, term = nothing, verbose = false)
    # Parameter priors
    h ~ Uniform(0.1, 1.0)
    Ωm0 ~ Uniform(0.0, 1.0)
    ΩΛ0 ~ Uniform(0.0, 1.0)
    Ωk0 = 1 - Ωr0 - Ωm0 - ΩΛ0

    p = [h, Ωm0, ΩΛ0, Ωk0]
    prob = probgen(p)
    as = 1 ./ (zs .+ 1)
    sol = solve(prob; bgopts, term, saveat = as, save_end = true)

    if !issuccess(sol)
        verbose && println("Discarding solution with illegal parameters Ωm0=$Ωm0, Ωk0=$Ωk0")
        Turing.@addlogprob! -Inf
        return nothing # illegal parameters
    end

    # Theoretical prediction
    dLs_pred = dL(sol, geta, gett, getΩk0, geth)[begin:end-1]

    # Compare predictions to data
    return dLs ~ MvNormal(dLs_pred, ΔdLs) # multivariate Gaussian # TODO: full covariance
end

p = [M.g.h, M.m.Ω₀, M.Λ.Ω₀, M.K.Ω₀]
probgen = SymBoltz.parameter_updater(prob, p; build_initializeprob = false)
sn = supernova(data.zs, data.dLs, data.ΔdLs, probgen, geta, gett, getΩk0, geth); # condition model on data
nothing # hide
```
We can now sample from the model to obtain a MCMC chain:
```@example fit
chain = sample(sn, NUTS(), 5000; initial_params = (h = 0.5, Ωm0 = 0.5, ΩΛ0 = 0.5)) # TODO: speed up: https://discourse.julialang.org/t/modelingtoolkit-odesystem-in-turing/115700/
import Plots, StatsPlots # don't collide with Makie
Plots.plot(chain)
```
Finally, we can visualize the high-dimensional chain with a corner plot using PairPlots.jl:
```@example fit
using PairPlots
layout = (
    PairPlots.Scatter(),
    PairPlots.Contourf(sigmas = 1:2),
    PairPlots.MarginHist(),
    PairPlots.MarginDensity(color = :black),
    PairPlots.MarginQuantileText(color = :black, font = :regular),
    PairPlots.MarginQuantileLines(),
)
pp = pairplot(chain => layout)
```

## Forecasting

Planning future cosmological surveys involves forecasting the precision of the constraints they are believed to place on parameters.
Here we show how one can perform forecasting by combining SymBoltz.jl with Turing.jl.

To start, create another probabilistic supernova model, but instead of observed luminosity distances, we now use simulated luminosity distances in a fiducial cosmology:
```@example fit
# TODO: motivate errors # hide
dLs_fid = dLs[begin:end-1]
sn_fc = supernova(data.zs, dLs_fid, data.ΔdLs, probgen, geta, gett, getΩk0, geth);
nothing # hide
```

### MCMC-driven forecasting

A general assumption-free, but expensive method to perform forecasting is to explore the likelihood using MCMC (as before, only against simulated data):
```@example fit
chain_fc = sample(sn_fc, NUTS(), 5000)
Plots.plot(chain_fc)
```
```@example fit
truth = PairPlots.Truth((h = pars[M.g.h], Ωm0 = pars[M.m.Ω₀], ΩΛ0 = pars[M.Λ.Ω₀]))
pp_fc = pairplot(chain_fc => layout, truth)
```

### Fisher forecasting

A less general, but cheaper method to perform forecasting is use a second-order approximation of the likelihood around the mean of the probability distribution (i.e. maximum of the likelihood).
This technique is called Fisher forecasting, and requires calculation of the **Fisher (information) matrix**
```math
Fᵢⱼ = -\left⟨\frac{∂^2 \log P(θ)}{∂θᵢ \, ∂θⱼ}\right⟩.
```
This effectively measures how quickly the likelihood falls from the maximum in different directions in parameter space.
Under certain assumptions, the Fisher matrix is the inverse of the covariance matrix ``C = F⁻¹``.

First, we ask Turing to [estimate the maximum likelihood mode](https://turinglang.org/docs/usage/mode-estimation/) of the probabilistic model:
```@example fit
maxl_fc = maximum_likelihood(sn_fc; initial_params = [0.5, 0.5, 0.5]) # TODO: or MAP?
#loglikelihood(sn_fc, (h = 0.70, Ωm0 = 0.26, Ωk0 = 0.10, w0 = -1.01, wa = -0.07)) # hide
@assert all(isapprox.(maxl_fc.values.array, [pars[M.g.h], pars[M.m.Ω₀], pars[M.Λ.Ω₀]]; atol = 1e-4)) # hide
maxl_fc # hide
```
As expected, the maximum likelihood corresponds to our chosen fiducial parameters.
Nevertheless, the returned mode estimate object offers convenience methods that greatly simplifies our following calculations.
Next, we calculate the Fisher matrix from the maximum likelihood, and invert it to get the covariance matrix:
```@example fit
using StatsBase
F = informationmatrix(maxl_fc)
C = inv(F)
```
Finally, we [derive 68% (1σ) and 95% (2σ) confidence ellipses from the covariance matrix](https://cookierobotics.com/007/) and draw them onto our corner plot:
```@example fit
# https://www.astronomy.ohio-state.edu/weinberg.21/A8824/stats4.pdf # hide
# TODO: "you can use MCMC instead, with your anticipated measurement errors and setting the data equal to the values expected for your fiducial model" # hide
# https://arxiv.org/pdf/1205.3984 # hide
# https://discourse.julialang.org/t/plot-ellipse-in-makie/82814/4 # hide
# https://docs.juliaplots.org/latest/generated/statsplots/#Covariance-ellipses # hide
# https://github.com/marcobonici/FisherPlot.jl/blob/main/src/FisherPlot.jl#L19 # hide
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

for i in eachindex(IndexCartesian(), C)
    ix, iy = i[1], i[2]
    ix >= iy && continue
    μx = maxl_fc.values[ix]
    μy = maxl_fc.values[iy]
    for nstd in 1:2
        xs, ys = ellipse(C.array, ix, iy, (μx, μy); nstd)
        lines!(pp_fc[iy,ix], xs, ys; color = :red)
    end
end

pp_fc
```

## Extension with w₀wₐ dark energy

```@example fit
g = SymBoltz.metric()
K = SymBoltz.curvature(g)
X = SymBoltz.w0wa(g; analytical = true)
M = RMΛ(K = K, Λ = X)
M = change_independent_variable(M, M.g.a; add_old_diff = true)
pars = Dict(M.r.Ω₀ => 9.3e-5, M.m.Ω₀ => 0.26, M.X.Ω₀ => 0.66, M.g.h => 0.7, M.X.wa => 0.0, M.X.w0 => -1.0, M.t => 0.0, M.r.T₀ => NaN, M.X.cₛ² => NaN)
pars[M.K.Ω₀] = 1 - pars[M.r.Ω₀] - pars[M.m.Ω₀] - pars[M.X.Ω₀]
prob = CosmologyProblem(M, pars; pt = false, ivspan = (1e-8, 1e0))

p = [M.g.h, M.m.Ω₀, M.X.Ω₀, M.X.w0, M.X.wa, M.K.Ω₀]
probgen = SymBoltz.parameter_updater(prob, p; build_initializeprob = false)
geta = getu(prob.bg, M.g.a) # must redefine because of new problem # TODO: move into package
gett = getu(prob.bg, M.t)
getΩk0 = getp(prob.bg, M.K.Ω₀)
geth = getp(prob.bg, M.g.h)

@model function supernova(zs, dLs, ΔdLs, probgen, geta, gett, getΩk0, geth; bgopts = bgopts, Ωr0 = 9.3e-5, term = nothing, verbose = false)
    # Parameter priors
    h ~ Uniform(0.1, 1.0)
    Ωm0 ~ Uniform(0.0, 1.0)
    ΩΛ0 ~ Uniform(0.0, 1.0)
    w0 ~ Uniform(-2.0, 0.0)
    wa ~ Uniform(-1.0, +1.0)
    Ωk0 = 1 - Ωr0 - Ωm0 - ΩΛ0

    p = [h, Ωm0, ΩΛ0, w0, wa, Ωk0]
    prob = probgen(p)
    as = 1 ./ (zs .+ 1)
    sol = solve(prob; bgopts, term, saveat = as, save_end = true)

    if !issuccess(sol)
        verbose && println("Discarding solution with illegal parameters Ωm0=$Ωm0, Ωk0=$Ωk0")
        Turing.@addlogprob! -Inf
        return nothing # illegal parameters
    end

    # Theoretical prediction
    dLs_pred = dL(sol, geta, gett, getΩk0, geth)[begin:end-1]

    # TODO: manual chi squared ala https://cmb.wintherscoming.no/milestone1.php
    #Turing.@addlogprob! -sum((dLs_pred[i] - dLs[i])^2 / ΔdLs[i]^2 for i in eachindex(dLs_pred))
    #return nothing

    # Compare predictions to data
    return dLs ~ MvNormal(dLs_pred, ΔdLs) # multivariate Gaussian # TODO: full covariance
end

sn = supernova(data.zs, data.dLs, data.ΔdLs, probgen, geta, gett, getΩk0, geth);
chain = sample(sn, NUTS(), 5000; initial_params = (h = 0.5, Ωm0 = 0.5, ΩΛ0 = 0.5, w0 = -1.0, wa = 0.0)) # TODO: speed up: https://discourse.julialang.org/t/modelingtoolkit-odesystem-in-turing/115700/
Plots.plot(chain)
```
```@example fit
truth = PairPlots.Truth((h = pars[M.g.h], Ωm0 = pars[M.m.Ω₀], ΩΛ0 = pars[M.X.Ω₀], w0 = pars[M.X.w0], wa = pars[M.X.wa]))
pp = pairplot(chain => layout, truth)
```

```@setup fit
using SymBoltz, OrdinaryDiffEq, Turing

function dL_fast(z, Ωm0, Ωk0, h; Ωr0 = 9.3e-5, aini = 1e-8, reltol = 1e-8, alg = Tsit5(), maxiters = 1e3)
    ΩΛ0 = 1 - Ωr0 - Ωm0 - Ωk0
    H0 = SymBoltz.H100 * h
    aH(a) = a * H0 * √(Ωr0/a^4 + Ωm0/a^3 + Ωk0/a^2 + ΩΛ0)
    function f(_, _, b)
        a = exp(b) # b = ln(a)
        return 1 / (aH(a)) # dt/db; t is conformal time
    end
    tini = 1 / (aH(aini))
    bini = log(aini)
    prob = ODEProblem(f, tini, (bini, 0))
    try # seems faster than going the NaNMath and retcode route
        sol = solve(prob; alg, reltol, maxiters)
        a = 1 ./ (z .+ 1)
        b = log.(a)
        t = sol(b)
        t0 = sol.u[end]
        χ = t0 .- t
        r = sinc.(√(-Ωk0+0im)*χ*H0/π) .* χ * SymBoltz.c .|> real # Julia's sinc(x) = sin(π*x) / (π*x)
        return r ./ a / SymBoltz.Gpc # from meters to Gpc
    catch
        return nothing
    end
end

dLs = dL_fast(data.zs, 0.26, 0.08, 0.7)

@model function supernova_fast(data)
    # Parameter priors
    Ωm0 ~ Uniform(0.0, 1.0)
    Ωk0 ~ Uniform(-1.0, +1.0)
    h ~ Uniform(0.5, 1.5)
    dLs = dL_fast(data.zs, Ωm0, Ωk0, h)
    if isnothing(dLs)
        Turing.@addlogprob! -Inf
        return nothing # illegal parameters
    end
    data.dLs ~ MvNormal(dLs, data.C) # multivariate Gaussian # TODO: full covariance
end

sn = supernova_fast(data)
chain = sample(sn, NUTS(), 500)
```
