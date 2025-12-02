# Performance and benchmarks

Model setup and hardware information:

```@example bench
using MKL, SymBoltz, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSDIRK, OrdinaryDiffEqBDF, BenchmarkTools, Plots, BenchmarkPlots, StatsPlots, LinearSolve
M = ΛCDM(ν = nothing, K = nothing, h = nothing)
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
benchmarks = BenchmarkGroup()

using Dates, InteractiveUtils, LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)
println("Current time: ", now(), "\n")
println(sprint(InteractiveUtils.versioninfo)) # show computer information
println(LinearAlgebra.BLAS.get_config())
println("BLAS threads: ", LinearAlgebra.BLAS.get_num_threads())
nothing # hide
```

## Background: ODE solver

This plot shows the time to solve the background using different (implicit) Rosenbrock methods with a fixed tolerance.

```@example bench
bgalgs = [
    Rodas4()
    Rodas5()
    Rodas4P()
    Rodas5P()
    FBDF() # unstable for some tolerances
    QNDF() # unstable for some tolerances
]
for alg in bgalgs
    bgopts = (alg = alg, reltol = 1e-9)
    benchmarks["bgsolver"][nameof(typeof(alg))] = @benchmarkable $solve($prob; bgopts = $bgopts)
end
results = run(benchmarks["bgsolver"]; verbose = true)
plot(results; size = (800, 400))
```

## Background: precision-work diagram

This plot compares the time to solve the background vs. accuracy of the solution using different ODE solvers and tolerances.
Every solution is compared to a reference solution with very small tolerance.
The points on each curve correspond to a sequence of tolerances.

```@example bench
# following e.g. https://github.com/SciML/ModelingToolkit.jl/issues/2971#issuecomment-2310016590
using DiffEqDevTools

refalg = Rodas4P()
bgsol = solve(prob.bg, refalg; abstol = 1e-12, reltol = 1e-12) # reference solution (results are similar compared to Rodas4/4P/5P/FBDF)

abstols = 1 ./ 10 .^ (5:9)
reltols = 1 ./ 10 .^ (5:9)
setups = [Dict(:alg => alg) for alg in bgalgs]
wp = WorkPrecisionSet(prob.bg, abstols, reltols, setups; appxsol = bgsol, save_everystep = false, error_estimate = :l2)
plot(wp; title = "Reference: $(SymBoltz.algname(refalg))", size = (800, 400), margin = 5*Plots.mm)
```

## Perturbations: ODE solver

This plot shows the time to solve several perturbation $k$-modes using different implicit ODE solvers with fixed tolerance.

```@example bench
# TODO: test different nlsolve # hide
ks = 10 .^ range(-1, 4, length = 100)
ptalgs = [
    TRBDF2()
    KenCarp4()
    KenCarp47()
    Kvaerno5()
    Rodas5P()
    Rodas5()
    Rodas4()
    Rodas4P()
    QNDF()
    FBDF()
]
for alg in ptalgs
    ptopts = (alg = alg, reltol = 1e-5)
    benchmarks["ptsolver"][nameof(typeof(alg))] = @benchmarkable $solve($prob, $ks; ptopts = $ptopts)
end
results = run(benchmarks["ptsolver"]; verbose = true)
plot(results; size = (800, 400))
```

## Perturbations: precision-work diagram

This plot compares the time to solve a perturbation $k$-mode vs. accuracy of the solution using different ODE solvers and tolerances.
Each subplot corresponds to a different $k$-mode.
Every solution is compared to a reference solution with very small tolerance.
The points on each curve correspond to a sequence of tolerances.

```@example bench
ks = [1e0, 1e1, 1e2, 1e3]
ptprob0, ptprobgen = SymBoltz.setuppt(prob.pt, bgsol)

refalg = Rodas4P()
setups = [Dict(:alg => alg) for alg in ptalgs]
wps = []
for (i, k) in enumerate(ks)
    ptprob = ptprobgen(ptprob0, k)
    ptsol = solve(ptprob, refalg; reltol = 1e-12, abstol = 1e-12) # reference solution (results are similar compared to QNDF/FBDF/Rodas4/4P/5P/, somewhat different with KenCarp4/Kvaerno5; use Rodas4P which is also used in CLASS comparison)
    wp = WorkPrecisionSet(ptprob, abstols, reltols, setups; appxsol = ptsol, save_everystep = false, error_estimate = :l2)
    push!(wps, wp)
end
p = plot(layout = (length(ks), 1), size = (800, 500*length(ks)))
for (i, k) in enumerate(ks)
    wp = wps[i]
    plot!(p[i,1], wp; title = "Reference: $(SymBoltz.algname(refalg)), k = $k H₀/c", left_margin = 15*Plots.mm, bottom_margin = 5*Plots.mm)
end
p
```

## Perturbations: time per mode

This plot shows the time spent solving individual perturbation $k$-modes using different ODE solvers with fixed tolerance.

```@example bench
ptprob0, ptprobgen = SymBoltz.setuppt(prob.pt, bgsol)
solvemode(k, ptalg) = solve(ptprobgen(ptprob0, k); alg = ptalg, reltol = 1e-7, abstol = 1e-7)

ks = 10 .^ range(-2, 5, length = 50)
times = [[minimum(@elapsed solvemode(k, ptalg) for i in 1:3) for k in ks] for ptalg in ptalgs]

plot(
    log10.(ks), map(ts -> log10.(ts), times); marker = :auto, markersize = 2,
    xlabel = "lg(k)", ylabel = "lg(time / s)", xticks = range(log10(ks[begin]), log10(ks[end]), step=1),
    label = permutedims(SymBoltz.algname.(ptalgs)), legend_position = :topleft
)
```

## Perturbations: Jacobian method

The Jacobian of the perturbation ODEs can be computed in three ways:

1. explicitly from a symbolically generated function,
2. numerically using forward-mode dual numbers, or
3. numerically using finite differences.

This plot shows the time to solve several perturbation $k$-modes for each such method.

```@example bench
ks = 10 .^ range(-2, 5, length = 200)
prob_nojac = CosmologyProblem(M, pars; jac = false)

bgopts = (alg = Rodas4P(linsolve = RFLUFactorization(),), reltol = 1e-9)
ptopts = (alg = KenCarp4(linsolve = RFLUFactorization(),), reltol = 1e-8) # generate function for J symbolically
benchmarks["jacobian"]["symbolic"] = @benchmarkable $solve($prob, $ks; bgopts = $bgopts, ptopts = $ptopts)

bgopts = (alg = Rodas4P(linsolve = RFLUFactorization(), autodiff = true), reltol = 1e-9)
ptopts = (alg = KenCarp4(linsolve = RFLUFactorization(), autodiff = true), reltol = 1e-8) # compute J with forward-mode AD
benchmarks["jacobian"]["forward diff"] = @benchmarkable $solve($prob_nojac, $ks; bgopts = $bgopts, ptopts = $ptopts)

bgopts = (alg = Rodas4P(linsolve = RFLUFactorization(), autodiff = true), reltol = 1e-9) # fails with finite diff background J
ptopts = (alg = KenCarp4(linsolve = RFLUFactorization(), autodiff = false), reltol = 1e-8) # compute J with finite differences
benchmarks["jacobian"]["finite diff"] = @benchmarkable $solve($prob_nojac, $ks; bgopts = $bgopts, ptopts = $ptopts)

results = run(benchmarks["jacobian"]; verbose = true)
plot(results; size = (800, 400))
```

## Perturbations: linear system solver

At every time step, the implicit perturbation ODE solver solves a system of (nonlinear) equations with Newton's method.
In turn, Newton's method involves iteratively solving several linear systems $Ax = b$, where $A$ involves the ODE Jacobian matrix $J$.
This can be a bottleneck, so it is very important to use a linear solver that does this as fast as possible!

**Note that the optimal linear solver depends on both the model and your hardware.**
See [this tutorial on accelerating linear solves](https://docs.sciml.ai/LinearSolve/stable/tutorials/accelerating_choices/).
Also run Julia with the [optimal BLAS backend for your platform](https://docs.julialang.org/en/v1/manual/performance-tips/#man-backends-linear-algebra),
such as [MKL.jl](https://github.com/JuliaLinearAlgebra/MKL.jl) for Intel or [AppleAccelerate.jl](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl) for Macs,
instead of Julia's default OpenBLAS backend.

In particular, for large models the ODE Jacobian can have lots of zeros.
Here is an example for a model with many perturbation equations.

```@example bench
lmaxs = [4, 8, 16, 32, 64]
Ms = [ΛCDM(ν = nothing, K = nothing, h = nothing; lmax) for lmax in lmaxs]

probs_dense = [CosmologyProblem(M, pars; ptopts = (jac = true, sparse = false)) for M in Ms]
probs_sparse = [CosmologyProblem(M, pars; ptopts = (jac = true, sparse = true)) for M in Ms]

probs_sparse[end].pt.f.jac_prototype # example of sparse Jacobian
```

Sparse matrix methods are therefore very important to speed up the solution of large perturbation systems.
This plot compares the time to solve several perturbation $k$-modes with different dense and sparse linear matrix solvers.

```@example bench
#import Sparspak # SparspakFactorization crashes # hide
#import Pardiso # hide
#using LinearSolve: SparspakFactorization # hide
#using LinearSolve: MKLPardisoFactorize # hide
# hide
ks = 10 .^ range(-2, 5, length=200)
ptopts_dense1 = (alg = KenCarp4(linsolve = LUFactorization()),)
ptopts_dense2 = (alg = KenCarp4(linsolve = RFLUFactorization()),)
ptopts_sparse1 = (alg = KenCarp4(linsolve = LUFactorization()),) # Base dispatch on SparseMatrixCSC
ptopts_sparse2 = (alg = KenCarp4(linsolve = KLUFactorization()),)
ptopts_sparse3 = (alg = KenCarp4(linsolve = UMFPACKFactorization()),)
ts_dense1 = [@belapsed solve($prob, $ks; ptopts = $ptopts_dense1) samples=1 evals=3 for prob in probs_dense]
ts_dense2 = [@belapsed solve($prob, $ks; ptopts = $ptopts_dense2) samples=1 evals=3 for prob in probs_dense]
ts_sparse1 = [@belapsed solve($prob, $ks; ptopts = $ptopts_sparse1) samples=1 evals=3 for prob in probs_sparse]
ts_sparse2 = [@belapsed solve($prob, $ks; ptopts = $ptopts_sparse2) samples=1 evals=3 for prob in probs_sparse]
ts_sparse3 = [@belapsed solve($prob, $ks; ptopts = $ptopts_sparse3) samples=1 evals=3 for prob in probs_sparse]

p1 = plot(ylabel = "time / s", xticks = (lmaxs, ""), ylims = (0.0, ceil(max(maximum(ts_dense1), maximum(ts_dense2)))))
marker = :circle
plot!(p1, lmaxs, ts_sparse1; label = "sparse $(nameof(typeof(ptopts_sparse1.alg.linsolve))), $(length(ks))×k", marker)
plot!(p1, lmaxs, ts_sparse2; label = "sparse $(nameof(typeof(ptopts_sparse2.alg.linsolve))), $(length(ks))×k", marker)
plot!(p1, lmaxs, ts_sparse3; label = "sparse $(nameof(typeof(ptopts_sparse3.alg.linsolve))), $(length(ks))×k", marker)
plot!(p1, lmaxs, ts_dense1; label = "dense $(nameof(typeof(ptopts_dense1.alg.linsolve))), $(length(ks))×k", marker)
plot!(p1, lmaxs, ts_dense2; label = "dense $(nameof(typeof(ptopts_dense2.alg.linsolve))), $(length(ks))×k", marker)
text(prob::CosmologyProblem) = "$(length(prob.pt.u0)) eqs,\n$(round(SymBoltz.sparsity_fraction(prob.pt)*100, digits=1)) %\nsparse"
annotate!(p1, lmaxs, zeros(length(lmaxs)), [(text(prob), 5, :top) for prob in probs_sparse])

speedups1 = [ts_dense2[i]/ts_sparse1[i] for i in eachindex(lmaxs)]
speedups2 = [ts_dense2[i]/ts_sparse2[i] for i in eachindex(lmaxs)]
speedups3 = [ts_dense2[i]/ts_sparse3[i] for i in eachindex(lmaxs)]
speedups4 = [ts_dense2[i]/ts_dense1[i] for i in eachindex(lmaxs)]
speedups5 = [ts_dense2[i]/ts_dense2[i] for i in eachindex(lmaxs)]
ymax = Int(ceil(maximum(maximum.([speedups1, speedups2, speedups3, speedups4, speedups5]))))
ylims = (0, ymax)
yticks = 0:1:ymax
yticks = (collect(yticks), collect("$y×" for y in yticks))
p2 = plot(; xlabel = "ℓmax", ylabel = "speedup", xticks = lmaxs, yticks, ylims, marker)
plot!(p2, lmaxs, speedups1; marker, label = nothing)
plot!(p2, lmaxs, speedups2; marker, label = nothing)
plot!(p2, lmaxs, speedups3; marker, label = nothing)
plot!(p2, lmaxs, speedups4; marker, label = nothing)
plot!(p2, lmaxs, speedups5; marker, label = nothing)

plot(p1, p2; size = (800, 600), layout = grid(2, 1, heights=(3//4, 1//4)))
```

Except for models with a very small perturbation system, it is a good idea to generate the sparse Jacobian and use the sparse `KLUFactorization` linear solver.

## Perturbations: parallelization

SymBoltz parallelizes integration of different perturbation $k$-modes with multithreading by default.
Make sure you [run Julia with multiple threads](https://docs.julialang.org/en/v1/manual/multi-threading/).
This is a standard technique in Boltzmann solvers, as linear perturbation modes are mathematically independent.
It leads to a performance improvement depending on the number of threads available.
It [can be disabled](@ref "Solving models"), for example if your application permits parallelization at a higher level.

```@example bench
using Base.Threads
ks = 10 .^ range(-1, 4, length = 100)
for thread in [true, false]
    label = thread ? "$(nthreads()) threads" : "1 thread"
    benchmarks["thread"][label] = @benchmarkable $solve($prob, $ks; thread = $thread)
end
results = run(benchmarks["thread"]; verbose = true)
plot(results; size = (800, 400))
```

```@setup
# TODO: tune Krylov with verbose = 1, ILU, ..., atol, rtol # hide
# TODO: KenCarp47(linsolve, precs = incompletelu) # hide
#ptsol = @btime solvept(prob.pt, bgsol, ks; alg = KenCarp47(linsolve = KrylovJL_GMRES(rtol = 1e-3, atol = 1e-3)), reltol = 1e-8) # hide
# TODO: optimize prob.pt.f.f.f_iip !!! lots of unnecessary stuff?? try cse = false and cse = true
# TODO: why is it solvept() slower than solvept(; output_func = (sol, i) -> (sol, false) ???
nothing # hide
```
