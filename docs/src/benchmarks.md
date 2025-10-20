# Performance and benchmarks

!!! note
    This page is a work in progress.

```@example bench
using MKL, SymBoltz, OrdinaryDiffEq, BenchmarkTools, Plots, BenchmarkPlots, StatsPlots, InteractiveUtils
M = SymBoltz.ΛCDM(ν = nothing, K = nothing, h = nothing)
pars = SymBoltz.parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
ks = 10 .^ range(-1, 4, length = 100)
benchmarks = BenchmarkGroup()

using Dates
println("Current time: ", now())
print(sprint(InteractiveUtils.versioninfo)) # show computer information
nothing # hide
```

## Parallelize perturbations

By default, SymBoltz parallelizes integration of different perturbation modes $k$ with multithreading.
This is a standard technique in Boltzmann solvers, as linear perturbation modes are mathematically independent.
It leads to a performance improvement depending on the number of threads available, but [can be turned off](@ref "Solving models"), for example if your application permits parallelization at a higher level:
```@example bench
using Base.Threads
for thread in [true, false]
    label = thread ? "$(nthreads()) threads" : "1 thread"
    benchmarks["thread"][label] = @benchmarkable $solve($prob, $ks; thread = $thread)
end
results = run(benchmarks["thread"]; verbose = true)
plot(results; size = (800, 400))
```

## Linear algebra backend

By default, SymBoltz uses implicit ODE solvers to integrate approximation-free stiff equations.
Unlike explicit solvers (used in Boltzmann solvers with approximations), they are often bottlenecked by linear algebra operations on the Jacobian matrix at every time step.
[Choosing an optimal linear matrix solver](https://docs.sciml.ai/LinearSolve/stable/tutorials/accelerating_choices/) and/or [switching to an alternative BLAS library](https://docs.julialang.org/en/v1/manual/performance-tips/#man-backends-linear-algebra) (such as [MKL](https://github.com/JuliaLinearAlgebra/MKL.jl) or [AppleAccelerate](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl) over Julia's default OpenBLAS backend) can therefore improve performance.
**This is dependent on your hardware!**

Moreover, both BLAS and SymBoltz use multi-threading to parallellize linear algebra operations and integration of indendent perturbation modes, respectively.
These can conflict. For best multi-threading performance [it is recommended to restrict BLAS to a single thread](https://docs.julialang.org/en/v1/manual/performance-tips/#man-multithreading-linear-algebra) and parallellize only the integration of perturbations.

```@example bench
linsolves = [
    SymBoltz.LUFactorization()
    SymBoltz.RFLUFactorization()
    #SymBoltz.MKLLUFactorization() # fails/hangs/segfaults # TODO: restore
    SymBoltz.KrylovJL_GMRES(rtol = 1e-3, atol = 1e-3)
]
for linsolve in linsolves
    ptopts = (alg = KenCarp4(; linsolve), reltol = 1e-8)
    benchmarks["linsolve"][nameof(typeof(linsolve))] = @benchmarkable $solve($prob, $ks; ptopts = $ptopts)
end
results = run(benchmarks["linsolve"]; verbose = true)
plot(results; size = (800, 400))
```

## Background solver options

```@example bench
bgalgs = [
    Rodas4()
    Rodas5()
    Rodas4P()
    Rodas5P()
    #FBDF() # does not work with lower tolerances
    #QNDF()
]
for alg in bgalgs
    bgopts = (alg = alg, reltol = 1e-9)
    benchmarks["bgsolver"][nameof(typeof(alg))] = @benchmarkable $solve($prob; bgopts = $bgopts)
end
results = run(benchmarks["bgsolver"]; verbose = true)
plot(results; size = (800, 400))
```

## Perturbation ODE solver

```@example bench
# TODO: test accuracy; work-precision diagram? # hide
# TODO: test different nlsolve # hide
# TODO: use BenchmarkTools.BenchmarkGroup # hide
ptalgs = [
    TRBDF2()
    #Rosenbrock23() # disabled; causes MaxIters with default tolerances
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

## Perturbations Jacobian method

```@example bench
probs = [
    CosmologyProblem(M, pars; jac = false)
    CosmologyProblem(M, pars; jac = true)
]
for prob in probs
    for autodiff in (true, false)
        numerical = isnothing(prob.bg.f.jac)
        numerical && !autodiff && continue # finite difference Jacobian fails # TODO: make work?
        name = numerical ? "Numerical" : "Symbolic"
        name *= autodiff ? " (auto. diff.)" : " (fin. diff.)"
        bgopts = (alg = SymBoltz.Rodas4P(autodiff = autodiff), reltol = 1e-9)
        ptopts = (alg = SymBoltz.KenCarp4(autodiff = autodiff), reltol = 1e-8)
        benchmarks["jacobian"][name] = @benchmarkable $solve($prob, $ks; bgopts = $bgopts, ptopts = $ptopts)
    end
end
results = run(benchmarks["jacobian"]; verbose = true)
plot(results; size = (800, 400))
```

```@setup
# TODO: tune Krylov with verbose = 1, ILU, ..., atol, rtol # hide
# TODO: KenCarp47(linsolve, precs = incompletelu) # hide
#ptsol = @btime solvept(prob.pt, bgsol, ks, prob.bgspline; alg = SymBoltz.KenCarp47(linsolve = SymBoltz.KrylovJL_GMRES(rtol = 1e-3, atol = 1e-3)), reltol = 1e-8) # hide
# TODO: make linsolve work with sparse Jacobian (jac = true, sparse = true)
# TODO: thread = true is not helping as much as it should!
# TODO: optimize prob.pt.f.f.f_iip !!! lots of unnecessary stuff?? try cse = false and cse = true
# TODO: https://docs.sciml.ai/DiffEqDocs/stable/tutorials/advanced_ode_example/
# TODO: https://docs.sciml.ai/ModelingToolkit/stable/examples/sparse_jacobians/
# TODO: set 1 BLAS thread because using implicit ode method
# TODO: 100% CPU usage with threads = false?
# TODO: why is it solvept() slower than solvept(; output_func = (sol, i) -> (sol, false) ???
nothing # hide
```

## Background precision-work diagram

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

## Perturbations precision-work diagram

```@example bench
ks = [1e0, 1e1, 1e2, 1e3]
ptprob0, ptprobgen = SymBoltz.setuppt(prob.pt, bgsol, prob.bgspline)

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
    plot!(p[i,1], wp; title = "Reference: $(SymBoltz.algname(refalg)), k = $k", left_margin = 15*Plots.mm, bottom_margin = 5*Plots.mm)
end
p
```

## Perturbations time per mode

```@example bench
ptprob0, ptprobgen = SymBoltz.setuppt(prob.pt, bgsol, prob.bgspline)
solvemode(k) = solve(ptprobgen(ptprob0, k); alg = SymBoltz.DEFAULT_PTALG, reltol = 1e-8, abstol = 1e-8)

ks = 10 .^ range(-2, 5, length = 100)
times = [minimum(@elapsed solvemode(k) for i in 1:10) for k in ks]

plot(log10.(ks), times .* 1000; xlabel = "lg(k)", ylabel = "time / ms", xticks = range(log10(ks[begin]), log10(ks[end]), step=1), marker = :circle, label = nothing)
```
