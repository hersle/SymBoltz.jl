# Benchmarks

!!! note
    This page is a work in progress.

```@example bench
using SymBoltz, Unitful, UnitfulAstro, OrdinaryDiffEq, BenchmarkTools, Plots
M = SymBoltz.ΛCDM(K = nothing)
pars = SymBoltz.parameters_Planck18(M) # TODO: much faster if i set ΩΛ0 explicitly?
prob = CosmologyProblem(M, pars)
N = 20
ks = 10 .^ range(-5, 1, length=N) / u"Mpc"

# TODO: make proper; sort; test accuracy; work-precision diagram?
# TODO: test different linsolve and nlsolve
# TODO: use BenchmarkTools.BenchmarkGroup
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
    map(t -> mean(t.times/N)/1e6, timings),
    yerror = map(t -> std(t.times/N)/1e6, timings),
    ylabel = "time per k-mode / ms", label = false,
    permute = (:x, :y) # make bar plot horizontal
)
```
