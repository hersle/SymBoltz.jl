import Base: nameof
import LinearAlgebra: issuccess, BLAS
import CommonSolve: solve
import SciMLBase: remake, successful_retcode
import PreallocationTools: DiffCache, get_tmp
import SciMLStructures
import SciMLStructures: canonicalize, Tunable
import OhMyThreads: TaskLocalValue
import SymbolicIndexingInterface
import SymbolicIndexingInterface: getsym, setsym_oop, parameter_values
using RecursiveFactorization # makes RFLUFactorization() available as linear solver: https://docs.sciml.ai/LinearSolve/stable/tutorials/accelerating_choices/
import NumericalIntegration: cumul_integrate
using SparseArrays

background(sys) = transform((sys, _) -> filter_system(isbackground, sys), sys)
perturbations(sys) = transform((sys, _) -> filter_system(isperturbation, sys), sys)

struct CosmologyProblem{Tbg <: ODEProblem, Tpt <: Union{ODEProblem, Nothing}, Tbgspline}
    M::System

    bg::Tbg
    pt::Tpt

    pars::Base.KeySet
    shoot::Base.KeySet
    conditions::AbstractArray

    bgspline::Tbgspline
end

struct CosmologySolution{Tbg <: ODESolution, Tpts <: Union{Nothing, EnsembleSolution, Vector{<:ODESolution}}, Tks <: Union{Nothing, AbstractVector}}
    prob::CosmologyProblem # problem which is solved
    bg::Tbg # background solution
    ks::Tks # perturbation wavenumbers
    pts::Tpts # perturbation solutions
end

algname(alg) = string(nameof(typeof(alg)))

function Base.show(io::IO, prob::CosmologyProblem; indent = "  ")
    print(io, "Cosmology problem for model ")
    printstyled(io, nameof(prob.bg.f.sys), '\n'; bold = true)

    printstyled(io, "Stages:"; bold = true)
    if !isnothing(prob.bg)
        print(io, '\n', indent, "Background")
        print(io, ": ", length(unknowns(prob.bg.f.sys)), " unknowns")
        print(io, ", ", 0, " splines")
        print(io, ", ", issparse(prob.bg) ? "$(round(sparsity_fraction(prob.bg)*100; digits=1)) % sparse" : "dense", " Jacobian")
    end
    if !isnothing(prob.pt)
        print(io, '\n', indent, "Perturbations")
        print(io, ": ", length(unknowns(prob.pt.f.sys)), " unknowns")
        print(io, ", ", isnothing(prob.bgspline) ? 0 : length(unknowns(prob.bg.f.sys)), " splines")
        print(io, ", ", issparse(prob.pt) ? "$(round(sparsity_fraction(prob.pt)*100; digits=1)) % sparse" : "dense", " Jacobian")
    end

    printstyled(io, "\nParameters & initial conditions:"; bold = true)
    for par in prob.pars
        val = getsym(prob, par)(prob)
        print(io, '\n', indent, "$par = $val", par in prob.shoot ? " (shooting)" : "")
    end

    !isempty(prob.conditions) && printstyled(io, "\nFinal conditions:"; bold = true) # TODO: initial vs final conditions
    for condition in prob.conditions
        print(io, '\n', indent, condition)
    end
end

function Base.show(io::IO, sol::CosmologySolution; indent = "  ")
    print(io, "Cosmology solution for model ")
    printstyled(io, nameof(sol.prob.bg.f.sys), '\n'; bold = true)

    retcode_color(retcode) = successful_retcode(retcode) ? :green : :red
    printstyled(io, "Stages:"; bold = true)
    if !isnothing(sol.bg)
        retcode = sol.bg.retcode
        print(io, '\n', indent, "Background: return code ")
        printstyled(io, retcode; color = retcode_color(retcode))
        print(io, "; solved with $(algname(sol.bg.alg)); $(length(sol.bg)) points")
    end
    if !isnothing(sol.pts)
        kmin, kmax = extrema(sol.ks)
        nmin, nmax = extrema(map(length, sol.pts))
        n = length(sol.pts)
        retcodes = unique(map(ptsol -> ptsol.retcode, sol.pts))
        print(io, '\n', indent, "Perturbations: return codes ")
        for (i, retcode) in enumerate(retcodes)
            printstyled(io, retcode; color = retcode_color(retcode))
            i < length(retcodes) && print(io, ", ")
        end
        print(io, "; solved with $(algname(sol.pts[1].alg)); $nmin-$nmax points; x$n k ∈ [$kmin, $kmax] H₀/c (interpolation in-between)")
    end
end

# Split parameters into DifferentialEquations' u0 and p convention
function split_vars_pars(M::System, x::Dict)
    pars = intersect(keys(x), parameters(M)) .|> ModelingToolkit.wrap # separate parameters from initial conditions # TODO: remove wrap
    vars = setdiff(keys(x), pars) # assume the rest are variables (do it without intersection to capture derivatives initial conditions)
    pars = Dict(par => x[par] for par in pars) # like p
    vars = Dict(var => x[var] for var in vars) # like u0
    return vars, pars
end

"""
    CosmologyProblem(
        M::System, pars::Dict, shoot_pars = Dict(), shoot_conditions = [];
        ivspan = (1e-6, 100.0), bg = true, pt = true, spline = true, debug = false, fully_determined = true, jac = true, sparse = false,
        bgopts = (), ptopts = (), kwargs...
    )

Create a numerical cosmological problem from the model `M` with parameters `pars`.
Optionally, the shooting method determines the parameters `shoot_pars` (mapped to initial guesses) such that the equations `shoot_conditions` are satisfied at the final time.

If `bg` and `pt`, the model is split into the background and perturbations stages.
If `spline` is a `Bool`, it decides whether all background unknowns in the perturbations system are replaced by splines.
If `spline` is a `Vector`, it rather decides which (unknown and observed) variables are splined.
If `jac`, analytic functions are generated for the ODE Jacobians; otherwise it is computed with forward-mode automatic differentiation by default.
If `sparse`, the perturbations ODE uses a sparse Jacobian matrix that is usually more efficient; otherwise a dense matrix is used.
"""
function CosmologyProblem(
    M::System, pars::Dict, shoot_pars = Dict(), shoot_conditions = [];
    ivspan = (1e-6, 100.0), bg = true, pt = true, spline = true, debug = false, fully_determined = true, jac = true, sparse = false,
    bgopts = (), ptopts = (), kwargs...
)
    pars = merge(pars, shoot_pars) # save full dictionary for constructor
    parsk = merge(pars, Dict(:k => NaN)) # k is unused, but must be set
    shoot_pars = keys(shoot_pars)

    if bg
        bg = background(M)
        bg = mtkcompile(bg)
        if debug
            bg = debug_system(bg)
        end

        # Set up callback for today # TODO: specify callbacks symbolically?
        terminate = !isnothing(term)
        iv = ModelingToolkit.get_iv(M)
        if Symbol(iv) == :τ
            aidx = ModelingToolkit.variable_index(bg, :a)
            f = (u, τ, integrator) -> u[aidx] - 1.0 # trigger callback when a = 1 (today)
        elseif Symbol(iv) == :a
            f = (u, a, integrator) -> a - 1.0 # a is independent variable
        else
            error("Don't know what to do when independent variable is $iv.")
        end
        parsymbols = Symbol.(parameters(bg))
        haveτ0 = Symbol(iv) == :τ && Symbol("τ0") in parsymbols
        haveκ0 = Symbol("b₊κ0") in parsymbols
        τ0idx = haveτ0 ? ModelingToolkit.parameter_index(bg, :τ0) : nothing
        _κidx = haveκ0 ? ModelingToolkit.variable_index(bg, M.b._κ) : nothing
        κ0idx = haveκ0 ? ModelingToolkit.parameter_index(bg, M.b.κ0) : nothing
        function affect!(integrator)
            if haveτ0
                integrator.ps[τ0idx] = integrator.t # set time today to time when a == 1 # TODO: what if τ is not iv
            end
            if haveκ0
                integrator.ps[κ0idx] = integrator.u[_κidx]
            end
            terminate && terminate!(integrator) # stop integration if desired
        end
        callback = ContinuousCallback(
            f, affect!;
            save_positions = (true, false), # don't duplicate final point
            rootfind = SciMLBase.RightRootFind # prefer right root, so a(τ₀) ≤ 1.0 and root finding algorithms get different signs also today (alternatively, try to enforce integrator.u[aidx] = 1.0 in affect! and set save_positions = (false, true), although this didn't work exactly last time)
        )

        bg = ODEProblem(bg, parsk, ivspan; fully_determined, callback, jac, bgopts..., kwargs...) # never sparse because small # TODO: hangs with jac = true, sparse = true
    else
        bg = nothing
    end

    if pt
        if spline == true
            spline = unknowns(bg.f.sys)
        elseif spline == false
            spline = []
        else
            # then spline should already be a vector of variables, so leave it unmodified
        end
        pt = perturbations(M)
        pt, bgspline = mtkcompile_spline(pt, spline)
        if debug
            pt = debug_system(pt)
        end
        # TODO: also remove_initial_conditions! from pt system (if initialization_eqs contain equations for splined variables)
        parsk = remove_initial_conditions!(parsk, spline) # must remove ICs of splined variables to avoid overdetermined initialization system
        pt = ODEProblem(pt, parsk, ivspan; fully_determined, jac, sparse, ptopts..., kwargs...)
    else
        pt = nothing
        bgspline = nothing
    end

    return CosmologyProblem(M, bg, pt, keys(pars), shoot_pars, shoot_conditions, bgspline)
end

"""
    function remake(
        prob::CosmologyProblem, pars::Dict;
        bg = true, pt = true, shoot = true,
        kwargs...
    )

Return an updated `CosmologyProblem` where parameters in `prob` are updated to values specified in `pars`.
Parameters that are not specified in `pars` keep their values from `prob`.
"""
function remake(
    prob::CosmologyProblem, pars::Dict;
    bg = true, pt = true, shoot = true,
    kwargs...
)
    vars, pars = split_vars_pars(prob.M, pars)
    vars = isempty(vars) ? missing : vars
    pars = isempty(pars) ? missing : pars
    bg = bg && !isnothing(prob.bg) ? remake(prob.bg; u0 = vars, p = pars, build_initializeprob = Val{!isnothing(prob.bg.f.initialization_data)}, kwargs...) : nothing
    if !ismissing(vars)
        vars = remove_initial_conditions!(vars, unknowns(prob.bg.f.sys)) # must filter ICs in remake, too
    end
    pt = pt && !isnothing(prob.pt) ? remake(prob.pt; u0 = vars, p = pars, build_initializeprob = Val{!isnothing(prob.pt.f.initialization_data)}, kwargs...) : nothing
    shoot_pars = shoot ? prob.shoot : keys(Dict())
    shoot_conditions = shoot ? prob.conditions : []
    return CosmologyProblem(prob.M, bg, pt, prob.pars, shoot_pars, shoot_conditions, prob.bgspline)
end

"""
    parameter_updater(prob::CosmologyProblem, idxs; kwargs...)

Create and return a function that updates the symbolic parameters `idxs` of the cosmological problem `prob`.
The returned function is called with numerical values (in the same order as `idxs`) and returns a new problem with the updated parameters.
"""
function parameter_updater(prob::CosmologyProblem, idxs; kwargs...)
    # define a closure based on https://docs.sciml.ai/ModelingToolkit/dev/examples/remake/#replace-and-remake
    # TODO: remove M, etc. for efficiency?

    @unpack bg, pt = prob

    bgsetsym = SymbolicIndexingInterface.setsym_oop(bg, idxs) # TODO: define setsym(::CosmologyProblem)?
    bgdiffcache = DiffCache(copy(canonicalize(Tunable(), parameter_values(bg))[1]))

    if !isnothing(pt)
        ptsetsym = setsym_oop(pt, idxs)
        ptdiffcache = DiffCache(copy(canonicalize(Tunable(), parameter_values(pt))[1]))
    end

    function updater(p)
        # Update background problem
        newu0, newp = bgsetsym(bg, p) # set new parameters
        bg_new = remake(bg; u0 = newu0, p = newp, kwargs...) # create updated problem (don't overwrite old)

        # Update perturbation problem
        if isnothing(pt)
            pt_new = pt
        else
            newu0, newp = ptsetsym(pt, p)
            pt_new = remake(pt; u0 = newu0, p = newp, kwargs...) # create updated problem (don't overwrite old)
        end

        return CosmologyProblem(prob.M, bg_new, pt_new, prob.pars, prob.shoot, prob.conditions, prob.bgspline)
    end
    function updater(p::Dict)
        p = [p[var] for var in idxs]
        return updater(p)
    end

    return updater
end

issparse(M::Nothing) = false
issparse(x) = SparseArrays.issparse(x)
issparse(prob::ODEProblem) = issparse(prob.f.jac_prototype)

function bgalg(prob::ODEProblem; stiff = true)
    if issparse(prob)
        linsolve = KLUFactorization()
    else
        linsolve = RFLUFactorization(throwerror = true)
    end
    if stiff
        return Rodas4P(; linsolve)
    else
        return Tsit5(; linsolve)
    end
end
bgalg(prob::CosmologyProblem; kwargs...) = bgalg(prob.bg; kwargs...)

function ptalg(prob::ODEProblem; accuracy = 1)
    if issparse(prob)
        linsolve = KLUFactorization()
    else
        linsolve = RFLUFactorization(throwerror = true)
    end
    nlsolve = NLNewton(fast_convergence_cutoff = 0, κ = 1)
    if accuracy == 0
        return TRBDF2(; linsolve, nlsolve)
    elseif accuracy == 1
        return KenCarp4(; linsolve, nlsolve)
    else
        return Rodas5P(; linsolve, nlsolve)
    end
end
ptalg(prob::CosmologyProblem; kwargs...) = ptalg(prob.pt; kwargs...)
shootalg(args...) = NewtonRaphson()

# TODO: want to use ODESolution's solver-specific interpolator instead of error-prone spline
"""
    function solve(
        prob::CosmologyProblem, ks::Union{Nothing, AbstractArray} = nothing;
        bgopts = (alg = bgalg(prob), reltol = 1e-9, abstol = 1e-9), bgextraopts = (),
        ptopts = (alg = ptalg(prob), reltol = 1e-8, abstol = 1e-8), ptextraopts = (),
        shootopts = (alg = shootalg(prob), abstol = 1e-5),
        thread = true, verbose = false, kwargs...
    )

Solve the cosmological problem `prob` up to the perturbative level with wavenumbers `ks` (or only to the background level if it is empty).
The options `bgopts` and `ptopts` are passed to the background and perturbations ODE `solve()` calls,
and `shootopts` to the shooting method nonlinear `solve()`.
If `threads`, integration over independent perturbation modes are parallellized.
"""
function solve(
    prob::CosmologyProblem, ks::Union{Nothing, AbstractArray} = nothing;
    bgopts = (alg = bgalg(prob), reltol = 1e-9, abstol = 1e-9), bgextraopts = (),
    ptopts = (alg = ptalg(prob), reltol = 1e-8, abstol = 1e-8), ptextraopts = (),
    shootopts = (alg = shootalg(prob), abstol = 1e-5),
    thread = true, verbose = false, kwargs...
)
    if !isempty(prob.shoot)
        bgsol = solvebg(prob.bg, collect(prob.shoot), prob.conditions; shootopts, verbose, bgopts..., bgextraopts..., kwargs...)
    else
        bgsol = solvebg(prob.bg; verbose, bgopts..., bgextraopts..., kwargs...)
    end

    if isnothing(ks) || isempty(ks)
        ks = nothing
        ptsol = nothing
    else
        ks = k_dimensionless.(ks, Ref(bgsol))
        ptsol = solvept(prob.pt, bgsol, ks, prob.bgspline; thread, verbose, ptopts..., ptextraopts..., kwargs...)
    end

    return CosmologySolution(prob, bgsol, ks, ptsol)
end
function solve(prob::CosmologyProblem, k::Number; kwargs...)
    return solve(prob, [k]; kwargs...)
end

function warning_failed_solution(sol::ODESolution, name = "ODE"; verbose = false)
    msg = "$name solution failed with return code $(sol.retcode)."
    if verbose
        t, u = sol.t[end], sol.u[end]
        msg *= " Final time and values:"
        msg *= "\n$(ModelingToolkit.get_iv(sol.prob.f.sys)) = $t"
        for (i, var) in enumerate(unknowns(sol.prob.f.sys))
            msg *= "\n$var = $(u[i])"
        end
    end
    msg *= "\nCheck the parameters and precision settings!"
    return msg
end

"""
    solvebg(bgprob::ODEProblem; alg = bgalg(bgprob), reltol = 1e-9, abstol = 1e-9, verbose = false, kwargs...)

Solve the background cosmology problem `bgprob`.
"""
function solvebg(bgprob::ODEProblem; alg = bgalg(bgprob), reltol = 1e-9, abstol = 1e-9, verbose = false, kwargs...)
    bgsol = solve(bgprob, alg; verbose, reltol, kwargs...)
    if !successful_retcode(bgsol)
        @warn warning_failed_solution(bgsol, "Background"; verbose)
    end

    τrecidx = ModelingToolkit.parameter_index(bgprob, :τrec)
    if !isnothing(τrecidx)
        bgsol.ps[τrecidx] = bgsol[:τ][argmax(bgsol[bgprob.f.sys.b.v])]
    end

    return bgsol
end
# TODO: more generic shooting method that can do anything (e.g. S8)
function solvebg(bgprob::ODEProblem, vars, conditions; alg = bgalg(bgprob), reltol = 1e-9, abstol = 1e-9, shootopts = (alg = shootalg(), reltol = 1e-3), verbose = false, build_initializeprob = Val{false}, kwargs...)
    length(vars) == length(conditions) || error("Different number of shooting parameters and conditions")

    setvars = SymbolicIndexingInterface.setsym_oop(bgprob, vars) # efficient setter
    getfuns = getsym(bgprob, map(eq -> eq.rhs - eq.lhs, conditions)) # efficient getter

    function f(vals, oldbgprob)
        #println("vals = $vals")

        # slow but "safe"
        #u0, p = SymBoltz.split_vars_pars(oldbgprob.f.sys, Dict(keys(vars) .=> vals))
        #newbgprob = remake(oldbgprob; u0, p)

        # fast but "unsafe"
        newu0, newp = setvars(oldbgprob, vals)
        newbgprob = remake(oldbgprob; u0 = newu0, p = newp, build_initializeprob)

        bgsol = solvebg(newbgprob; alg, reltol, abstol, kwargs..., save_everystep = false, save_start = false, save_end = true)
        return only(getfuns(bgsol)) # get final values
    end

    guess = map(var -> getsym(bgprob, var)(bgprob), vars)
    prob = NonlinearProblem(f, guess, bgprob)
    sol = solve(prob; show_trace = Val(verbose), shootopts...)

    u0, p = setvars(bgprob, sol.u)
    bgprob = remake(bgprob; u0, p, build_initializeprob)
    return solvebg(bgprob; alg, reltol, abstol, kwargs...)
end

function setuppt(ptprob::ODEProblem, bgsol::ODESolution, bgsplinepar)
    splset! = ModelingToolkit.setsym_oop(ptprob, [bgsplinepar])
    kset! = ModelingToolkit.setp(ptprob, k)

    ivspan = (bgsol.t[begin], bgsol.t[end])
    bgspline = spline(bgsol)
    newu0, newp = splset!(ptprob, Any[bgspline]) # TODO: why Vector{Any} needed to make solve() inferred?
    SciMLStructures.replace!(Tunable(), newp, canonicalize(Tunable(), parameter_values(bgsol))[1]) # copy background parameters to perturbations (e.g. τ0 and κ0)
    ptprob = remake(ptprob; tspan = ivspan, u0 = newu0, p = newp)

    #ptprob_tlv = TaskLocalValue{ODEProblem}(() -> remake(ptprob0; u0 = copy(ptprob0.u0) #= p is copied below =#)) # prevent conflicts where different tasks modify same problem: https://discourse.julialang.org/t/solving-ensembleproblem-efficiently-for-large-systems-memory-issues/116146/11 (alternatively copy just p and u0: https://github.com/SciML/ModelingToolkit.jl/issues/3056) # TODO: copy u0, p only?

    return ptprob, (ptprob, k) -> begin
        #ptprob = ptprob_tlv[]
        p = copy(newp) # newp specializes on spline types, while ptprob0.p does not; see https://github.com/SciML/ModelingToolkit.jl/issues/3715
        #p = copy(ptprob0.p) # see https://github.com/SciML/ModelingToolkit.jl/issues/3346 and https://github.com/SciML/ModelingToolkit.jl/issues/3056 # TODO: copy only Tunables
        kset!(p, k)
        ptprob = remake(ptprob; u0 = newu0, p = p, build_initializeprob = true) # solve for u0 # TODO: separate function?
        ptprob = remake(ptprob; u0 = ptprob.u0, p = p, build_initializeprob = false) # remake again with build_initializeprob = false makes following solve type-stable; https://github.com/SciML/ModelingToolkit.jl/issues/3715
        #println("Parameter type: ", typeof(ptprob.p))
        return ptprob
    end
end

"""
    solvept(ptprob::ODEProblem, bgsol::ODESolution, ks::AbstractArray, bgsplinepar; alg = ptalg(ptprob), reltol = 1e-8, abstol = 1e-8, output_func = (sol, i) -> sol, thread = true, verbose = false, kwargs...)

Solve the perturbation cosmology problem `ptprob` with wavenumbers `ks`.
A background solution `bgsol` must be passed (see `solvebg`), and a parameter `bgsplinepar` that refers to a spline in the perturbation problem of background unknowns.
If `thread` and Julia is running with multiple threads, the solution of independent wavenumbers is parallellized.
The return value is a vector with one `ODESolution` per wavenumber, or its mapping through `output_func` if a custom transformation is passed.
"""
function solvept(ptprob::ODEProblem, bgsol::ODESolution, ks::AbstractArray, bgsplinepar; alg = ptalg(ptprob), reltol = 1e-8, abstol = 1e-8, output_func = (sol, i) -> sol, thread = true, verbose = false, kwargs...)
    !issorted(ks) && throw(error("ks = $ks are not sorted in ascending order"))

    if thread && Threads.nthreads() == 1
        thread = false
        @warn "Multi-threading over perturbation modes was requested, but disabled, since Julia is running with only 1 thread. Restart Julia with more threads (e.g. `julia --threads=auto`) to enable multi-threading, or pass thread = false to explicitly disable it." maxlog=1
    end
    if thread && BLAS.get_num_threads() > 1
        @warn "Multi-threading over perturbation modes was requested, but BLAS is running with $(BLAS.get_num_threads()) threads.\nIt is recommended to restrict BLAS to one thread with `using LinearAlgebra: BLAS; BLAS.set_num_threads(1)`.\nFor more information, see https://docs.julialang.org/en/v1/manual/performance-tips/#man-multithreading-linear-algebra." maxlog=1
    end

    # TODO: can I exploit that the structure of the perturbation ODEs is ẏ = J * y with "constant" J?
    ptprob, ptprobgen = setuppt(ptprob, bgsol, bgsplinepar)

    function output_func_warn(sol, i)
        if !successful_retcode(sol)
            @warn warning_failed_solution(sol, "Perturbation (mode k = $(ks[i]))"; verbose)
        else
            verbose && print("\rSolved perturbations for wavenumber k = $(ks[i])")
        end
        return output_func(sol, i)
    end
    ptsols = [@spawnif output_func_warn(solve(ptprobgen(ptprob, ks[i]), alg; verbose, reltol, abstol, kwargs...), i) thread for i in eachindex(ks)]
    ptsols = fetch.(ptsols)
    verbose && println()
    return ptsols
end
"""
    solvept(ptprob::ODEProblem; alg = ptalg(ptprob), reltol = 1e-8, abstol = 1e-8, kwargs...)

Solve the perturbation problem `ptprob` and return the solution.
Its wavenumber and background spline must already be initialized, for example with `setuppt`.

# Examples

```julia
# ...
prob = CosmologyProblem(M, pars)
bgsol = solvebg(prob.bg)
ptprob0, ptprobgen = SymBoltz.setuppt(prob.pt, bgsol, prob.bgspline)
k = 1.0
ptprob = ptprobgen(ptprob0, k)
ptsol = solvept(ptprob)
```
"""
function solvept(ptprob::ODEProblem; alg = ptalg(ptprob), reltol = 1e-8, abstol = 1e-8, kwargs...)
    return solve(ptprob, alg; reltol, abstol, kwargs...)
end

function time_today(prob::CosmologyProblem)
    getτ0 = SymBoltz.getsym(prob.bg, :τ0)
    bgprob = prob.bg
    bgsol = solvebg(bgprob; save_everystep = false, save_start = false, save_end = true)
    return getτ0(bgsol)
end

"""
    issuccess(sol::CosmologySolution)

Returns whether the solution of a cosmological problem was successful (i.e. not failing due to instability or too many time steps).
"""
function issuccess(sol::CosmologySolution)
    return successful_retcode(sol.bg) && (isnothing(sol.pts) || all(successful_retcode(pt) for pt in sol.pts))
end

function integrate(xs, ys; integrator = TrapezoidalRule())
    prob = SampledIntegralProblem(ys, xs)
    sol = solve(prob, integrator)
    return sol.u
end
integrate_cumulative(sol::CosmologySolution, x, y) = cumul_integrate(sol[x], sol[y])
integrate_cumulative(sol::CosmologySolution, y) = integrate_cumulative(sol, sol.prob.M.τ, y)

# TODO: don't select time points as 2nd/3rd index, since these points will vary
const SymbolicIndex = Union{Num, AbstractArray{Num}}
function Base.getindex(sol::CosmologySolution, i::SymbolicIndex)
    if ModelingToolkit.isparameter(i) && !isequal(i, ModelingToolkit.get_iv(sol.prob.M)) # don't catch independent variable as parameter
        return sol.bg.ps[i] # assume all parameters are in background # TODO: index sol directly when this is fixed? https://github.com/SciML/ModelingToolkit.jl/issues/3267
    else
        return sol.bg[i]
    end
end
function Base.getindex(sol::CosmologySolution, i::SymbolicIndex, j)
    return stack(sol.bg[i, j])
end
Base.getindex(sol::CosmologySolution, i::Int, j::SymbolicIndex, k = :) = sol.pts[i][j, k]
Base.getindex(sol::CosmologySolution, i, j::SymbolicIndex, k = :) = [stack(sol[_i, j, k]) for _i in i]
Base.getindex(sol::CosmologySolution, i::Colon, j::SymbolicIndex, k = :) = sol[1:length(sol.pts), j, k]

"""
    express_derivatives(expr, prob::CosmologyProblem)

Express derivatives in the symbolic expression `expr` in terms of non-differentiated quantities in the system `sys`.
"""
function express_derivatives(expr, prob)
    # Create a map of as many derivatives-to-expressions that we know
    bg = prob.bg.f.sys
    dvarmap_bg = merge(
        Dict(var => eq.lhs for (var, eq) in map_variables_to_equations(bg) if is_derivative(var)), # observed dummy derivatives
        Dict(eq.lhs => eq.rhs for eq in equations(bg)) # ODE unknowns
    )
    pt = prob.pt.f.sys
    dvarmap_pt = merge(
        Dict(var => eq.lhs for (var, eq) in map_variables_to_equations(pt) if is_derivative(var)),
        Dict(eq.lhs => eq.rhs for eq in equations(pt))
    )
    while true
        expr0 = expr
        expr = expand_derivatives(expr) # isolate derivatives
        expr = substitute(expr, dvarmap_bg) # substitute derivatives
        expr = substitute(expr, dvarmap_pt) # substitute derivatives
        expr === expr0 && break # stop when expression doesn't change anymore
    end
    Symbolics.hasnode(is_derivative, expr) && error("Could not express derivative of $expr")
    return expr
end

"""
    get_is_deriv(prob::CosmologyProblem, is)

Make any transformation of symbolic indices `is` before querying an `ODESolution` with them and the best derivative order `deriv`.
"""
function get_is_deriv(prob::CosmologyProblem, is)
    arederivs = Symbolics.is_derivative.(unwrap.(is))
    if all(arederivs) && all(map(i -> Symbolics.contains_var(only(Symbolics.arguments(unwrap(i))), unknowns(prob)), is))
        is = map(i -> only(Symbolics.arguments(unwrap(i))), is) # peel off derivative operators, ...
        deriv = Val{1} # ... but request 1st derivative
    elseif any(arederivs)
        # expand derivatives in terms of non-differentiated variables
        is = map(i -> express_derivatives(i, prob), is)
        deriv = Val{0}
    else # no derivatives
        deriv = Val{0} # request 0th derivative; keep indices as-is (usual case)
    end
    return is, deriv
end

# TODO: match variable convention (i.e. δ(τ, k))
function (sol::CosmologySolution)(is::AbstractArray, ts::AbstractArray)
    #tmin, tmax = extrema(sol.bg.t[[begin, end]])
    #minimum(ts) >= tmin || minimum(ts) ≈ tmin || throw("Requested time $(minimum(ts)) is before initial time $tmin")
    #maximum(ts) <= tmax || maximum(ts) ≈ tmax || throw("Requested time $(maximum(ts)) is after final time $tmax")
    is, deriv = get_is_deriv(sol.prob, is)
    return sol.bg(ts, deriv; idxs=is)[:, :]
end
(sol::CosmologySolution)(i::Num, ts::AbstractArray) = sol([i], ts)[1, :]
(sol::CosmologySolution)(is::AbstractArray, t::Number) = sol(is, [t])[:, 1]
(sol::CosmologySolution)(i::Num, t::Number) = sol([i], [t])[1, 1]

# similar to https://github.com/SciML/SciMLBase.jl/blob/c568c0eb554ba78440a83792f058073c286a55d3/src/solutions/ode_solutions.jl#L277
function getfunc(sol::ODESolution, var; deriv = Val{0}, continuity = :left)
    ps = SciMLBase.parameter_values(sol)
    return t -> begin
        state = SciMLBase.ProblemState(; u = sol.interp(t, nothing, deriv, ps, continuity), p = ps, t)
        return getsym(sol, var)(state)
    end
end

function getsym(provider::Union{CosmologyProblem, CosmologySolution}, p)
    getsym_bg = SymbolicIndexingInterface.getsym(provider.bg, p)
    return provider -> getsym_bg(provider.bg)
end

function neighboring_modes_indices(sol::CosmologySolution, k)
    k = k_dimensionless.(k, Ref(sol.bg))
    if k == sol.ks[begin] # k == kmin
        i1 = i2 = 1
    elseif k == sol.ks[end] # k == kmax
        i1 = i2 = length(sol.ks)
    else
        i2 = searchsortedfirst(sol.ks, k) # index above target k
        i1 = i2 - 1 # index below target k
    end
    return i1, i2
end

Base.eltype(sol::CosmologySolution) = eltype(sol.bg)

function (sol::CosmologySolution)(out::AbstractArray, is::AbstractArray, ts::AbstractArray, ks::AbstractArray; smart = true, ktransform = log)
    if isnothing(sol.ks) || isempty(sol.ks)
        throw(error("No perturbations solved for. Pass ks to solve()."))
    end
    ks = k_dimensionless.(ks, Ref(sol.bg))
    kmin, kmax = extrema(sol.ks)
    minimum(ks) >= kmin || throw("Requested wavenumber k = $(minimum(ks)) is below the minimum solved wavenumber $kmin")
    maximum(ks) <= kmax || throw("Requested wavenumber k = $(maximum(ks)) is above the maximum solved wavenumber $kmax")
    # disabled; sensitive to exact endpoint values
    #tmin, tmax = extrema(sol.th.t)
    #minimum(ts) >= tmin || throw("Requested time $(minimum(ts)) is below minimum solved time $tmin")
    #maximum(ts) <= tmax || throw("Requested time $(maximum(ts)) is above maximum solved time $tmin")

    is, deriv = get_is_deriv(sol.prob, is)

    # Pre-allocate intermediate and output arrays
    v = similar(sol.bg, length(is), length(ts))
    v1 = similar(sol.bg, length(is), length(ts))
    v2 = similar(sol.bg, length(is), length(ts))

    i1_prev, i2_prev = -1, -1 # cache previous looked up solution and reuse it, if possible
    for ik in eachindex(ks) # TODO: multithreading leads to trouble; what about tmap?
        k = ks[ik]
        # Find two wavenumbers to interpolate between
        i1, i2 = neighboring_modes_indices(sol, k)

        # Evaluate solutions for neighboring wavenumbers,
        # but reuse those from the previous iteration if we are still between the same neighboring wavenumbers
        if i1 == i2_prev && smart
            v1 .= v2 # just set to v2 when incrementing i1 by 1
            i1_prev = i2_prev
        elseif i1 != i1_prev || !smart
            v1 .= sol.pts[i1](ts, deriv; idxs=is) # https://docs.sciml.ai/DiffEqDocs/latest/basics/solution/ # TODO: allocate less or make in-place (https://github.com/SciML/OrdinaryDiffEq.jl/issues/2562)
            i1_prev = i1
        end
        if i2 != i2_prev || !smart
            v2 .= sol.pts[i2](ts, deriv; idxs=is) # TODO: getu or similar for speed? possible while preserving interpolation?
            i2_prev = i2
        end
        v .= v1
        if i1 != i2
            # interpolate between solutions
            k1 = sol.ks[i1]
            k2 = sol.ks[i2]
            w = (ktransform(k) - ktransform(k1)) / (ktransform(k2) - ktransform(k1)) # interpolate between some function of the wavenumbers between ktransform(k) (e.g. k -> k or k -> log(k)) # TODO: cubic spline?
            @. v += (v2 - v1) * w # add to v1 from above
        end
        for ii in eachindex(is)
            out[ii, :, ik] .= v[ii, :]
        end
    end

    return out
end
function (sol::CosmologySolution)(is::AbstractArray, ts::AbstractArray, ks::AbstractArray; kwargs...)
    out = similar(sol.bg, length(is), length(ts), length(ks))
    return sol(out, is, ts, ks; kwargs...)
end
(sol::CosmologySolution)(is::AbstractArray, ts::AbstractArray, k::Number; kwargs...) = sol(is, ts, [k]; kwargs...)[:, :, 1]
(sol::CosmologySolution)(is::AbstractArray, t::Number, ks::AbstractArray; kwargs...) = sol(is, [t], ks; kwargs...)[:, 1, :]
(sol::CosmologySolution)(i::Num, ts::AbstractArray, ks::AbstractArray; kwargs...) = sol([i], ts, ks; kwargs...)[1, :, :]
(sol::CosmologySolution)(is::AbstractArray, t::Number, k::Number; kwargs...) = sol(is, [t], [k]; kwargs...)[:, 1, 1]
(sol::CosmologySolution)(i::Num, ts::AbstractArray, k::Number; kwargs...) = sol([i], ts, [k]; kwargs...)[1, :, 1]
(sol::CosmologySolution)(i::Num, t::Number, ks::AbstractArray; kwargs...) = sol([i], [t], ks; kwargs...)[1, 1, :]
(sol::CosmologySolution)(i::Num, t::Number, k::Number; kwargs...) = sol([i], [t], [k]; kwargs...)[1, 1, 1]

function (sol::CosmologySolution)(is, tmap::Pair)
    tvar, ts = tmap
    ts = timeseries(sol, tvar, ts)
    return sol(is, ts)
end

function (sol::CosmologySolution)(is, tmap::Pair, ks)
    tvar, ts = tmap
    ts = timeseries(sol, tvar, ts)
    return sol(is, ts, ks)
end

function timeseries(sol::CosmologySolution; kwargs...)
    ts = sol.bg.t
    return timeseries(ts; kwargs...)
end
function timeseries(sol::CosmologySolution, k; kwargs...)
    i1, i2 = neighboring_modes_indices(sol, k)
    t1s = sol.pts[i1].t
    t2s = sol.pts[i2].t
    ts = sort!(unique!([t1s; t2s])) # average or interleave?
    return timeseries(ts; kwargs...)
end
function timeseries(ts::AbstractArray; Nextra = 0)
    if Nextra > 0
        ts = extend_array(ts, Nextra)
    end
    return ts
end
"""
    timeseries(sol::CosmologySolution, var, vals; kwargs...)

Find the times when some variable `var` equals some values `vals` with a spline.
"""
function timeseries(sol::CosmologySolution, var, vals; alg = ITP(), kwargs...)
    allequal(sign.(diff(sol[var]))) || error("$var is not monotonic")
    varfunc = getfunc(sol.bg, var)
    f(t, p) = varfunc(t) - p # var(t) == val when f(t) == 0
    ivspan = extrema(sol.bg.t)
    prob = IntervalNonlinearProblem(f, ivspan, vals[1]; kwargs...)
    return map(val -> solve(remake(prob; p = val); alg).u, vals)
end
"""
    timeseries(sol::CosmologySolution, var, dvar, vals::AbstractArray)

Find the times when some variable `var` equals some values `vals` with a Hermite spline, also taking the derivative `dvar` of `var` into account.
"""
function timeseries(sol::CosmologySolution, var, dvar, vals::AbstractArray; kwargs...)
    ts = timeseries(sol; kwargs...)
    xs = sol(var, ts)
    ẋs = sol(dvar, ts) # dx/dt
    ṫs = 1 ./ ẋs # dt/dx
    spl = spline(ts, ṫs, xs)
    ts = spl(vals)
    return ts
end

"""
    unknowns(prob::CosmologyProblem)

Get all unknown variables from the background and perturbations of the cosmological problem `prob`.
"""
function unknowns(prob::CosmologyProblem)
    return [unknowns(prob.bg.f.sys); unknowns(prob.pt.f.sys)]
end

"""
    parameters(prob::CosmologyProblem; bg = true, pt = true, spline = false)

Get all parameter values of the cosmological problem `prob`.
"""
function parameters(prob::CosmologyProblem; bg = true, pt = true, spline = false)
    bg = bg && !isnothing(prob.bg) ? parameters(prob.bg) : Dict()
    pt = pt && !isnothing(prob.pt) ? parameters(prob.pt) : Dict()
    pars = merge(bg, pt)
    !spline && delete!.(Ref(pars), values(prob.bgspline))
    return pars
end
function parameters(prob::ODEProblem)
    pars = parameters(prob.f.sys)
    return Dict(pars .=> prob.ps[pars])
end
function parameters(sol::CosmologySolution; kwargs...)
    return parameters(sol.prob; kwargs...)
end

# Fix model/solution under broadcasted calls
Base.broadcastable(sys::System) = Ref(sys)
Base.broadcastable(sol::CosmologySolution) = Ref(sol)

# Statistics for solution of background
function statsbg(sol::CosmologySolution)
    return sol.bg.stats
end

# Summarized statistics for solution of all perturbation modes
function statspt(sol::CosmologySolution)
    stats = SciMLBase.DEStats(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NaN)
    for field in fieldnames(typeof(stats))
        @eval $stats.$field = sum(ptsol.stats.$field for ptsol in $sol.pts)
    end
    return stats
end

function sparsity_fraction(J::SparseMatrixCSC)
    nall = length(J)
    nnonzeros = nnz(J)
    nzeros = nall - nnonzeros
    return nzeros / nall
end
function sparsity_fraction(prob::ODEProblem)
    J = prob.f.jac_prototype
    if isnothing(J)
        return 0.0 # matrix is dense
    else
        return sparsity_fraction(J)
    end
end
