# Benchmarks

!!! note
    This page is a work in progress.

```@example bench
using SymBoltz, Unitful, UnitfulAstro, OrdinaryDiffEq, LinearSolve, BenchmarkTools
M = SymBoltz.ΛCDM(h = nothing)
pars = SymBoltz.parameters_Planck18(M)
ks = 10 .^ range(-5, 1, length=10) / u"Mpc"
```

```@example bench
# TODO: make proper; sort; test accuracy; work-precision diagram?
# TODO: test different linsolve and nlsolve
reltol = 1e-5
optss = [
    (solver = TRBDF2,),
    (solver = KenCarp4,),
    (solver = KenCarp47,),
    (solver = Kvaerno5,),
    (solver = Rodas5P,), # TODO: Rodas5?
    (solver = QNDF,),
    (solver = QBDF,),
]
for opts in optss
    print("$(opts.solver) (reltol = $(reltol)): ")
    solver = opts.solver()
    timing = @benchmark sol = solve($M, $pars, $ks; solver = $solver, thread=false)
    println("($(mean(timing).time / 1e9) ± $(std(timing).time)) s")
end
```