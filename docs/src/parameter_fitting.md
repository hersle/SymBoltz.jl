# Fitting parameters and forecasting

This tutorial shows how to perform Bayesian parameter inference on a cosmological model by fitting it to data, and how to forecast parameter constraints from a covariance matrix that describes uncertainties and correlations of (unknown) observed data.

## Pantheon supernova data

We load the binned [Pantheon dataset](https://github.com/dscolnic/Pantheon/).
This includes redshifts and apparent magnitudes of over 1000 Type Ia supernovae.
```@example fit
using DataFrames, CSV, LinearAlgebra, PDMats

binned = true # use compressed dataset in the documentation
if binned # choose compressed dataset with 40 redshift bins
    data = "https://github.com/dscolnic/Pantheon/raw/master/Binned_data/lcparam_DS17f.txt"
    Csyst = "https://github.com/dscolnic/Pantheon/raw/master/Binned_data/sys_DS17f.txt"
else # choose full dataset with 1048 supernovae
    data = "https://github.com/dscolnic/Pantheon/raw/master/lcparam_full_long.txt"
    Csyst = "https://github.com/dscolnic/Pantheon/raw/master/sys_full_long.txt"
end

# Read data table
data = download(data)
data = CSV.read(data, DataFrame, delim = " ", silencewarnings = true)

# Read covariance matrix of apparent magnitudes (mb)
Csyst = download(Csyst)
Csyst = CSV.read(Csyst, DataFrame, header = false) # long vector
Csyst = collect(reshape(Csyst[2:end, 1], (Int(Csyst[1, 1]), Int(Csyst[1, 1])))) # to matrix
Cstat = Diagonal(data.dmb)^2 # TODO: should this be squared?
C = Csyst + Cstat

# Sort data and covariance matrix with decreasing redshift
is = sortperm(data, :zcmb, rev = true)
C = C[is, is]
C = PDMat(Symmetric(C)) # efficient sym-pos-def matrix with Cholesky factorization
data = data[is, :]
first(data, 10) # show 10 first rows
```
Let us plot the apparent magnitude as a function of redshift, and the covariance matrix:
```@example fit
using CairoMakie
fig = Figure(size = (600, 800))
ax1 = Axis(fig[1, 1:2], xlabel = "z", ylabel = "m", title = "Apparent brightness vs. redshift")
scatter!(ax1, data.zcmb, data.mb; markersize = 5, label = "data (Pantheon)")
errorbars!(ax1, data.zcmb, data.mb, data.dmb; linewidth = 1, whiskerwidth = 5)
ax2 = Axis(fig[2, 1]; xlabel = "z", ylabel = "z", title = "Covariance matrix", yreversed = true, aspect = 1)
hm = heatmap!(ax2, extrema(data.zcmb), extrema(data.zcmb), C; colormap = :balance, colorrange = (-0.001, +0.001))
Colorbar(fig[2, 2], hm)
fig
```

## Predicting luminosity distances

To predict luminosity distances
```math
d_L = \frac{r}{a} = \chi \, \mathrm{sinc} (\sqrt{k} \chi),
\qquad \text{where} \qquad
\chi = c \, (t_0 - t)
```
theoretically, we solve the w0waCDM model:
```@example fit
using SymBoltz, ModelingToolkit
g = SymBoltz.metric()
K = SymBoltz.curvature(g)
X = SymBoltz.w0wa(g; analytical = true)
M = RMΛ(K = K, Λ = X)
M = change_independent_variable(M, M.g.a; add_old_diff = true)
pars_fixed = Dict(M.t => 0.0, M.r.T₀ => NaN, M.X.cₛ² => NaN)
pars_varying = [M.r.Ω₀, M.m.Ω₀, M.K.Ω₀, M.X.Ω₀, M.g.h, M.X.w0, M.X.wa]

dL = SymBoltz.distance_luminosity_function(M, pars_fixed, pars_varying, data.zcmb)
μ(p) = 5 * log10.(dL(p)[begin:end-1] / (10*SymBoltz.pc)) # distance modulus

# Show example predictions
Mb = -19.3 # absolute supernova brightness (constant since SN-Ia are standard candles)
bgopts = (alg = SymBoltz.Tsit5(), reltol = 1e-5, maxiters = 1e3)
p0 = [9.3e-5, 0.3, 0.0, 0.7, 0.7, -1.0, 0.0] # fiducial parameters
μs = μ(p0)
mbs = μs .+ Mb
lines!(ax1, data.zcmb, mbs; color = :black, label = "theory (ΛCDM)")
axislegend(ax1, position = :rb)
fig
```

## Bayesian inference

To perform bayesian inference, we define a probabilistic model in [Turing.jl](https://turinglang.org/):
```@example fit
using Turing

@model function supernova(μ_pred, mbs, C; Mb = Mb, Ωr0 = 9.3e-5)
    # Parameter priors
    h ~ Uniform(0.1, 1.0)
    Ωm0 ~ Uniform(0.0, 1.0)
    Ωk0 ~ Uniform(-1.0, +1.0)
    w0 ~ Uniform(-2.0, 0.0)
    wa ~ Uniform(-1.0, +1.0)
    ΩX0 = 1 - Ωr0 - Ωm0 - Ωk0

    p = [Ωr0, Ωm0, Ωk0, ΩX0, h, w0, wa]
    μs_pred = μ_pred(p)
    if isempty(μs_pred)
        Turing.@addlogprob! -Inf
        return nothing
    end
    mbs_pred = μs_pred .+ Mb
    return mbs ~ MvNormal(mbs_pred, C) # read "measurements sampled from multivariate normal with predictions and covariance matrix"

    # equivalently:
    #Δmb = mbs .- mbs_pred
    #χ² = transpose(Δmb) * invC * Δmb
    #Turing.@addlogprob! -1/2 * χ²
    #return nothing
end

# https://github.com/JuliaStats/Distributions.jl/issues/1964 # TODO: get rid of? PR?
function MvNormal(μ::AbstractVector{<:Real}, Σ::AbstractPDMat{<:Real})
    R = Base.promote_eltype(μ, Σ)
    Distributions.MvNormal{R, typeof(Σ), typeof(μ)}(μ, Σ)
end
function MvNormal(μ, Σ)
    return Distributions.MvNormal(μ, Σ)
end

sn_w0waCDM = supernova(μ, data.mb, C);
sn_ΛCDM = fix(sn_w0waCDM, w0 = -1.0, wa = 0.0);
nothing # hide
```
We can now sample from the model to obtain a MCMC chain for the ΛCDM model:
```@example fit
chain = sample(sn_ΛCDM, NUTS(), 1000; initial_params = (h = 0.5, Ωm0 = 0.5, Ωk0 = 0.0))
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
We can easily repeat this for another model:
```@example fit
sn_w0CDM_flat = fix(sn_w0waCDM, Ωk0 = 0.0, wa = 0.0);
# TODO: describe Turing model more, e.g. loglikelihood(sn_fc, (h = 0.70, Ωm0 = 0.26, Ωk0 = 0.10, w0 = -1.01, wa = -0.07)) # hide
chain = sample(sn_w0CDM_flat, NUTS(), 1000; initial_params = (h = 0.5, Ωm0 = 0.5, w0 = -1.0))
pp = pairplot(chain => layout)
```

## Forecasting

Planning future cosmological surveys involves forecasting the precision of the constraints they are believed to place on parameters.
Here we show how one can perform forecasting by combining SymBoltz.jl with Turing.jl.

To start, create another probabilistic supernova model, but instead of observed luminosity distances, we now use simulated luminosity distances in a fiducial cosmology:
```@example fit
sn_fc_w0waCDM = supernova(μ, mbs, Diagonal(C));
sn_fc_w0CDM_flat = fix(sn_fc_w0waCDM, Ωk0 = 0.0, wa = 0.0);
nothing # hide
```

### MCMC-driven forecasting

A general assumption-free, but expensive method to perform forecasting is to explore the likelihood using MCMC (as before, only against simulated data):
```@example fit
chain_fc = sample(sn_fc_w0CDM_flat, NUTS(), 1000)
Plots.plot(chain_fc)
```
```@example fit
pars0 = Dict(pars_varying .=> p0)
truth = PairPlots.Truth((h = pars0[M.g.h], Ωm0 = pars0[M.m.Ω₀], w0 = pars0[M.X.w0]))
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
maxl_fc = maximum_likelihood(sn_fc_w0CDM_flat; initial_params = [0.5, 0.5, -1.0]) # TODO: or MAP?
@assert all(isapprox.(maxl_fc.values.array, [pars0[M.g.h], pars0[M.m.Ω₀], pars0[M.X.w0]]; atol = 1e-4)) # hide
maxl_fc # hide
```
As expected, the maximum likelihood corresponds to our chosen fiducial parameters.
Nevertheless, the returned mode estimate object offers convenience methods that greatly simplifies our following calculations.
Next, we calculate the Fisher matrix from the maximum likelihood, and invert it to get the covariance matrix:
```@example fit
using StatsBase
F_fc = informationmatrix(maxl_fc)
C_fc = inv(F_fc)
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

for i in eachindex(IndexCartesian(), C_fc)
    ix, iy = i[1], i[2]
    ix >= iy && continue
    μx = maxl_fc.values[ix]
    μy = maxl_fc.values[iy]
    for nstd in 1:2
        xs, ys = ellipse(C_fc.array, ix, iy, (μx, μy); nstd)
        lines!(pp_fc[iy,ix], xs, ys; color = :red)
    end
end

pp_fc
```

```@setup fit
#=
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
=#
nothing
```
