# Fitting parameter with Bayesian inference

This tutorial shows how to perform Bayesian parameter inference on a cosmological model by fitting it to data.

## Observed supernova data

We use the data from [Betoule+ (2014)](https://arxiv.org/abs/1401.4064) of redshifts and luminosity distances of supernovae.
Specifically, we use the values and errors in tables F.1 and F.2
(available from the [data files `data/jla_mub.txt` and `data/jla_mub_covmatrix.dat`](http://supernovae.in2p3.fr/sdss_snls_jla/jla_likelihood_v6.tgz)):
```@example 1
# TODO: do things properly
# TODO: use full covariance
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
data = (zs = data[:,1], dLs = data[:,2], ΔdLs = data[:,3], C = Diagonal(data[:,3] .^ 2))

using Plots
scatter(@. log10(data.zs+1), @. data.dLs / data.zs; yerror = data.ΔdLs ./ data.zs, xlabel = "lg(1+z)", ylabel = "dₗ/z/Gpc", label = "observations")
```

## Predicting luminosity distances

To predict luminosity distances
```math
d_L = \frac{r}{a} = \chi \, \mathrm{sinc} (\sqrt{k} \chi),
\qquad \text{where} \qquad
\chi = c \, (t_0 - t)
```
theoretically, we solve the standard ΛCDM model:
```@example 1
# TODO: generalize to non-flat
using SymBoltz, DataInterpolations

M = ΛCDM(recombination=false) # TODO: build m-only model

# TODO: improve performance of this
function dL(z, sol::CosmologySolution)
    as = sol[M.sys.g.a][begin:end-1]
    ts = sol[SymBoltz.t][begin:end-1] / sol.bg.ps[M.sys.g.H0] # s
    t_of_loga = CubicSpline(ts, log.(as); extrapolate=true) # TODO: do internally
    a0, a = 1, 1 ./ (z .+ 1)
    t0, t = t_of_loga(log(a0)), t_of_loga(log.(a))
    r = @. SymBoltz.c * (t0 - t)
    return @. r / a / SymBoltz.Gpc
end

function dL(z, M::CosmologyModel, Ωm0, h)
    return dL(z, solve(M, [
        M.sys.γ.Ω0 => 5e-5,
        M.sys.ν.Neff => 3.0,
        M.sys.b.Ω0 => 0.0,
        M.sys.c.Ω0 => Ωm0, # TODO: create matter-only model
        M.sys.g.h => h
    ]))
end

# Show example predictions
zs = data.zs
dLs = dL(zs, M, 0.3, 0.7)
scatter!(@. log10(zs+1), dLs ./ zs; label = "prediction")
```

## Bayesian inference

To perform bayesian inference, we use the [Turing.jl](https://turinglang.org/) package:
```@example 1
using Turing

@model function supernova(data, M)
    # Parameter priors
    Ωm0 ~ Uniform(0.2, 0.4)
    h ~ Uniform(0.6, 0.8)

    # Theoretical prediction
    dLs = dL(data.zs, M, Ωm0, h)

    # Compare predictions to data
    data.dLs ~ MvNormal(dLs, data.C) # multivariate Gaussian # TODO: full covariance
end

# TODO: speed up: https://discourse.julialang.org/t/modelingtoolkit-odesystem-in-turing/115700/
sn = supernova(data, M) # condition model on data
chain = sample(sn, MH(), MCMCSerial(), 500, 1) # TODO: NUTS() # TODO: MCMCThreads()
```
As we see above, the MCMC `chain` displays a summary with information about the fitted parameters, including their posterior means and standard deviations.
We can also plot the chains:
```@example 1
using StatsPlots
plot(chain)
```
