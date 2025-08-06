# Performance and benchmarks

!!! note
    This page is a work in progress.

```@example bench
using MKL, SymBoltz, OrdinaryDiffEq, BenchmarkTools, Plots, BenchmarkPlots, StatsPlots, InteractiveUtils
M = SymBoltz.ΛCDM(ν = nothing, K = nothing, h = nothing)
pars = SymBoltz.parameters_Planck18(M) # TODO: faster if i set ΩΛ0 explicitly?
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
    SymBoltz.MKLLUFactorization()
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
    KenCarp4()
    KenCarp47()
    Kvaerno5()
    Rodas5P()
    Rodas5()
    Rodas4()
    Rodas4P()
    QNDF()
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
#ptsol = @btime solvept(prob.pt, bgsol, ks, prob.var2spl; alg = SymBoltz.KenCarp47(linsolve = SymBoltz.KrylovJL_GMRES(rtol = 1e-3, atol = 1e-3)), reltol = 1e-8) # hide
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

