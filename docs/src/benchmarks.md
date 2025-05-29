# Performance and benchmarks

!!! note
    This page is a work in progress.

## Linear algebra backend

By default, SymBoltz uses implicit ODE solvers to integrate approximation-free stiff equations.
Unlike explicit solvers (used in Boltzmann solvers with approximations), they are often bottlenecked by linear algebra operations on the Jacobian matrix at every time step.
[Choosing an optimal linear matrix solver](https://docs.sciml.ai/LinearSolve/stable/tutorials/accelerating_choices/) and/or [switching to an alternative BLAS library](https://docs.julialang.org/en/v1/manual/performance-tips/#man-backends-linear-algebra) (such as [MKL](https://github.com/JuliaLinearAlgebra/MKL.jl) or [AppleAccelerate](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl) over Julia's default OpenBLAS backend) can therefore improve performance.
**This is dependent on your hardware!**

Moreover, both BLAS and SymBoltz use multi-threading to parallellize linear algebra operations and integration of indendent perturbation modes, respectively.
These can conflict. For best multi-threading performance [it is recommended to restrict BLAS to a single thread](https://docs.julialang.org/en/v1/manual/performance-tips/#man-multithreading-linear-algebra) and parallellize only the integration of perturbations.

## ODE solver

```@example bench
using SymBoltz, OrdinaryDiffEq, BenchmarkTools, Plots
M = SymBoltz.ΛCDM(ν = nothing, K = nothing, h = nothing)
pars = SymBoltz.parameters_Planck18(M) # TODO: much faster if i set ΩΛ0 explicitly?
prob = CosmologyProblem(M, pars)
ks = 10 .^ range(-1, 4, length = 100)

# TODO: make proper; sort; test accuracy; work-precision diagram? # hide
# TODO: test different linsolve and nlsolve # hide
# TODO: use BenchmarkTools.BenchmarkGroup # hide
algs = [TRBDF2(), AutoTsit5(TRBDF2()), KenCarp4(), AutoTsit5(KenCarp4()), KenCarp47(), Kvaerno5(), Rodas5P(), Rodas5(), Rodas4(), Rodas4P(), QNDF()]
solve_with(alg) = solve(prob, ks; ptopts = (alg, reltol = 1e-5), thread = false)
timings = map(algs) do alg
    return @benchmark solve_with($alg)
end

# Sort by efficiency
idxs = sortperm(map(t -> mean(t.times), timings))
algs, timings = algs[idxs], timings[idxs]

bar(
    SymBoltz.algname.(algs),
    map(t -> mean(t.times/length(ks))/1e6, timings),
    yerror = map(t -> std(t.times/length(ks))/1e6, timings),
    ylabel = "time per k-mode / ms", label = false,
    permute = (:x, :y) # make bar plot horizontal
)
```

## ODE solver options

```@example bench
bgsol = @btime solvebg(prob.bg)
nothing # hide
```
```@example bench
ptsol = @btime solvept(prob.pt, bgsol, ks, prob.var2spl; alg = SymBoltz.KenCarp4(autodiff = true))
nothing # hide
```
```@example bench
ptsol = @btime solvept(prob.pt, bgsol, ks, prob.var2spl; alg = SymBoltz.KenCarp4(autodiff = false))
nothing # hide
```
```@example bench
ptsol = @btime solvept(prob.pt, bgsol, ks, prob.var2spl; alg = SymBoltz.KenCarp4(linsolve = SymBoltz.LUFactorization()))
nothing # hide
```
```@example bench
ptsol = @btime solvept(prob.pt, bgsol, ks, prob.var2spl; alg = SymBoltz.KenCarp4(linsolve = SymBoltz.RFLUFactorization()))
nothing # hide
```
```@example bench
ptsol = @btime solvept(prob.pt, bgsol, ks, prob.var2spl; alg = SymBoltz.KenCarp4(linsolve = SymBoltz.KrylovJL_GMRES()))
nothing # hide
```
```@example bench
ptsol = @btime solvept(prob.pt, bgsol, ks, prob.var2spl; alg = SymBoltz.KenCarp47())
nothing # hide
```
```@example bench
# TODO: tune Krylov with verbose = 1, ILU, ..., atol, rtol # hide
# TODO: KenCarp47(linsolve, precs = incompletelu) # hide
ptsol = @btime solvept(prob.pt, bgsol, ks, prob.var2spl; alg = SymBoltz.KenCarp47(linsolve = SymBoltz.KrylovJL_GMRES(rtol = 1e-3, atol = 1e-3)), reltol = 1e-8)
nothing # hide
```
```@setup
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

