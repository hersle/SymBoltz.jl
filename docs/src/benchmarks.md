# Benchmarks

!!! note
    This page is a work in progress.

```@example bench
using SymBoltz, Unitful, UnitfulAstro, OrdinaryDiffEq, LinearSolve, BenchmarkTools, Plots
M = SymBoltz.ΛCDM(h = nothing)
pars = SymBoltz.parameters_Planck18(M)
N = 20
ks = 10 .^ range(-5, 1, length=N) / u"Mpc"

# TODO: make proper; sort; test accuracy; work-precision diagram?
# TODO: test different linsolve and nlsolve
# TODO: use BenchmarkTools.BenchmarkGroup
solvers = [TRBDF2(), AutoTsit5(TRBDF2()), KenCarp4(), AutoTsit5(KenCarp4()), KenCarp47(), Kvaerno5(), Rodas5P(), QNDF()] # TODO: Rodas5?
timings = []
solve_with(solver; reltol = 1e-5) = solve(M, pars, ks; solver, reltol, thread=false)
for solver in solvers
    timing = @benchmark solve_with($solver)
    push!(timings, timing)
end
bar(
    SymBoltz.solvername.(solvers),
    map(t -> mean(t.times/N)/1e6, timings),
    yerror = map(t -> std(t.times/N)/1e6, timings),
    ylabel = "time per k-mode / ms", label = false,
    permute = (:x, :y) # make bar plot horizontal
)
```
