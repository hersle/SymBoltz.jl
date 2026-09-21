import Base: nameof
import LinearAlgebra: issuccess, BLAS
import CommonSolve: solve
import SciMLBase: remake, successful_retcode
import SciMLLogging
import PreallocationTools: DiffCache, get_tmp
import SciMLStructures
import SciMLStructures: canonicalize, Tunable
import OhMyThreads: TaskLocalValue
import SymbolicIndexingInterface
import SymbolicIndexingInterface: getsym, setsym_oop, parameter_values, variable_index, parameter_index
using RecursiveFactorization # makes RFLUFactorization() available as linear solver: https://docs.sciml.ai/LinearSolve/stable/tutorials/accelerating_choices/
import NumericalIntegration: cumul_integrate
using SparseArrays
import NonlinearSolve.BracketingNonlinearSolve: AbstractBracketingAlgorithm

background(sys) = transform((sys, _) -> filter_system(isbackground, sys), sys)
perturbations(sys) = transform((sys, _) -> filter_system(isperturbation, sys), sys)

struct CosmologyProblem{Tbg <: Tuple{Vararg{ODEProblem}}, Tpt <: Union{ODEProblem, Nothing}}
    M::System

    bg::Tbg # background stages solved in order, each with the unknowns of the previous stages splined in (the last has all background variables)
    pt::Tpt

    pars::Vector{Symbolics.SymbolicT}
    shoot::Dict
    conditions::Vector{Equation} # shooting conditions in the form lhs - rhs ~ 0
    terminate::Union{Nothing, Equation} # event that stops the first forwards integration stage
end

struct CosmologySolution{Tbg <: Tuple{Vararg{ODESolution}}, Tpts <: Union{Nothing, EnsembleSolution, Vector{<:ODESolution}}, Tks <: Union{Nothing, AbstractVector}}
    prob::CosmologyProblem # problem which is solved
    bg::Tbg # background stage solutions
    ks::Tks # perturbation wavenumbers
    pts::Tpts # perturbation solutions
end

algname(alg) = string(nameof(typeof(alg)))

isbackwards(prob::ODEProblem) = prob.tspan[end] < prob.tspan[begin] # TODO: assumes iv is e.g. τ or log(a), but not e.g. z
isbackwards(sol::ODESolution) = sol.t[end] < sol.t[begin]

# Print the unknowns of a stage with n shortest names (usually the most fundamental variables)
function show_unknowns(io::IO, prob::ODEProblem; n = 3)
    vars = sort(string.(unknowns(prob.f.sys)); by = length)
    print(io, " (", join(first(vars, n), ", ")) # restrict to n variables
    length(vars) > n && print(io, ", …")
    print(io, ")")
end

# Print the pairs of a mapping on separate lines like a Dict (with aligned keys)
function show_mapping(io::IO, mapping; color = false)
    isempty(mapping) && return
    context = IOContext(io, :limit => true, :displaysize => (typemax(Int), displaysize(io)[2]), :color => color) # limit makes "=>" align; typemax(Int) shows all rows
    str = sprint(show, "text/plain", Dict(mapping); context) # print Dict mapping to a string
    print(io, str[findfirst('\n', str):end]) # remove first line with "Dict{...} with N entries:"
end

function Base.show(io::IO, prob::CosmologyProblem; indent = "  ", compact = true, bold = true, color = false)
    io = IOContext(io, :compact => compact) # print numbers compactly, like in arrays
    symio = IOContext(io, :color => color) # print symbolic expressions without syntax highlighting
    print(io, "Cosmology problem for model ")
    printstyled(io, nameof(prob.M), '\n'; bold)

    iv = ModelingToolkit.get_iv(prob.M)
    tmin, tmax = extrema(prob.bg[1].tspan) # in either direction
    printstyled(io, "Timespan:"; bold)
    print(io, " from ", iv, " = ", tmin, " till ", iv, " = ", tmax)
    !isnothing(prob.terminate) && print(symio, " or ", prob.terminate)

    printstyled(io, "\nStages:"; bold)
    for (i, stage) in enumerate(prob.bg)
        print(io, '\n', indent, "Background $i: ", isbackwards(stage) ? "backwards" : "forwards")
        nvars = length(unknowns(stage.f.sys))
        print(io, ", ", nvars, " unknowns")
        nvars > 0 && show_unknowns(io, stage)
        print(io, ", ", issparse(stage) ? "$(round(sparsity_fraction(stage)*100; digits=1)) % sparse" : "dense", " Jacobian")
    end
    if !isnothing(prob.pt)
        print(io, '\n', indent, "Perturbations")
        nvars = length(unknowns(prob.pt.f.sys))
        print(io, ": ", nvars, " unknowns")
        nvars > 0 && show_unknowns(io, prob.pt)
        print(io, ", ", issparse(prob.pt) ? "$(round(sparsity_fraction(prob.pt)*100; digits=1)) % sparse" : "dense", " Jacobian")
    end

    printstyled(io, "\nIndependent parameters:"; bold)
    show_mapping(io, Dict(par => getsym(prob, par)(prob) for par in prob.pars if !(par in keys(prob.shoot)))) # shooting parameters are printed below

    !isempty(prob.shoot) && printstyled(io, "\nShooting guesses:"; bold)
    show_mapping(io, prob.shoot)

    !isempty(prob.conditions) && printstyled(io, "\nShooting conditions:"; bold)
    for condition in prob.conditions
        print(symio, '\n', indent, condition)
    end
end

function Base.show(io::IO, sol::CosmologySolution; indent = "  ", bold = true)
    print(io, "Cosmology solution for model ")
    printstyled(io, nameof(sol.prob.M), '\n'; bold)

    retcode_color(retcode) = successful_retcode(retcode) ? :green : :red
    printstyled(io, "Stages:"; bold)
    for (i, bgsol) in enumerate(sol.bg)
        retcode = bgsol.retcode
        print(io, '\n', indent, "Background $i: return code ")
        printstyled(io, retcode; color = retcode_color(retcode))
        print(io, "; solved with $(algname(bgsol.alg)); $(length(bgsol.u)) points")
    end
    if !isnothing(sol.pts)
        kmin, kmax = extrema(sol.ks)
        nmin, nmax = extrema(map(ptsol -> length(ptsol.u), sol.pts))
        n = length(sol.pts)
        retcodes = unique(map(ptsol -> ptsol.retcode, sol.pts))
        print(io, '\n', indent, "Perturbations: return codes ")
        for (i, retcode) in enumerate(retcodes)
            printstyled(io, retcode; color = retcode_color(retcode))
            i < length(retcodes) && print(io, ", ")
        end
        print(io, "; solved with $(algname(sol.pts[1].alg)); $nmin-$nmax points; x$n k ∈ [$kmin, $kmax]")
    end
end

# Require that only parameters are specified (initial conditions must be set through parameters declared in the model)
function check_parameters(pars::Dict)
    nonpars = filter(!ModelingToolkit.isparameter, collect(keys(pars)))
    isempty(nonpars) || error("Only parameters can be specified, but got $(join(nonpars, ", ")). To specify an initial condition, declare a parameter for it in the model with initial_conditions = [x => xini].")
end

# Select the value of an option for background stage i of n, which is either shared by all stages or given per stage in a Tuple
function stageopt(opt, i, n)
    opt isa Tuple && !isempty(opt) || return opt
    length(opt) == n || error("Got $(length(opt)) values $opt for $n background stages")
    return opt[i]
end
stageopts(opts, i, n) = map(opt -> stageopt(opt, i, n), NamedTuple(opts))

"""
    CosmologyProblem(
        M::System, pars::Dict, shoot_pars = Dict(), shoot_conditions = [];
        tspan = (1e-6, 100.0), terminate = M.a ~ 1,
        bg = true, pt = true, spline = true, debug = false, fully_determined = true,
        bgjac = true, bgsparse = false, bgopts = (),
        ptjac = true, ptsparse = true, ptopts = (),
        iip = true, specialize = SciMLBase.AutoSpecialize,
        kwargs...
    )

Create a numerical cosmological problem from the model `M` with parameters `pars`.

Optionally, the shooting method determines the parameters `shoot_pars` (mapped to initial guesses) such that the equations `shoot_conditions` are satisfied at the final time.
Shooting parameters and conditions declared in `M` are included automatically, and guesses in `shoot_pars` override those in `M`.

The background is solved in stages given by the Tuple `bg` of variable vectors, each interpolating previous stages with splines.
If `bg = true`, the stages are detected from the dependencies between background variables with `SymBoltz.split_stages`.
Stages with variables declared with `[backwards = true]` (like `χ` and `κ`) are integrated backwards from today, and other stages forwards.
Each option in `bgopts` is a single value for all stages, or a Tuple with one value per stage.

The first stage is integrated over `tspan`, and later stages over the span of the previous stage.
The first forwards stage terminates at the event `terminate` (default today when ``a = 1``); pass `terminate = nothing` to integrate over all of `tspan`.

If `pt = false`, or if `M` has no wavenumber parameter `k`, the perturbations are not created.
The extra options `bgopts` and `ptopts` are passed to the `ODEProblem` constructors of the background stages and perturbations, and override the prefixed options above.
The unprefixed options in `kwargs` (like `jac` or `sparse`) are passed to all of them last, and override all of them.

If `spline` is a `Bool`, it decides whether all background unknowns in the perturbations system are replaced by splines.
If `spline` is a `Vector`, it rather decides which (unknown and observed) variables are splined.

If `bgjac`/`ptjac`, analytic functions are generated for the background/perturbation ODE Jacobians; otherwise they are computed with forward-mode automatic differentiation by default.
If `bgsparse`/`ptsparse`, the background/perturbation ODEs use sparse Jacobian matrices that are usually more efficient for large systems; otherwise dense matrices are used.
By default the perturbations are sparse, while the smaller background stages are dense.

If `fully_determined`, the initialization system of every stage must have as many equations as unknowns.
If `debug`, the system of every stage is wrapped with `ModelingToolkit.debug_system` to help locate errors in the equations.

The [SciMLBase type parameters](https://docs.sciml.ai/SciMLBase/stable/interfaces/Problems/) `iip` and `specialize` are forwarded to internal `ODEProblem{iip, specialize}(...)` constructors.
"""
function CosmologyProblem(
    M::System, pars::Dict, shoot_pars = Dict(), shoot_conditions = [];
    tspan = (1e-6, 100.0), terminate = M.a ~ 1,
    bg = true, pt = true, spline = true, debug = false, fully_determined = true,
    bgjac = true, bgsparse = false, bgopts = (),
    ptjac = true, ptsparse = true, ptopts = (),
    iip = true, specialize = SciMLBase.AutoSpecialize,
    kwargs...
)
    p_constructor(buf) = convert(Vector{isempty(buf) ? eltype(buf) : typeof(first(buf))}, buf) # converts nonnumeric Any vector to vector of concrete spline type
    check_parameters(pars)
    shoot_pars_sys = shootvars(M)
    conditions_sys = ModelingToolkit.get_constraints(M)
    shoot_pars = mergesafe(shoot_pars_sys, shoot_pars) # read from system, but let passed guesses override them
    shoot_conditions = unionsafe(shoot_conditions, conditions_sys) # read from system

    length(shoot_pars) != length(shoot_conditions) && error("Got $(length(shoot_pars)) shooting parameters ($(join(keys(shoot_pars), ", "))), but $(length(shoot_conditions)) conditions ($(join(shoot_conditions, ", ")))")
    length(shoot_pars) > 1 && any(x -> x[2] isa Tuple, shoot_pars) && error("Shooting with multiple parameters requires scalar guesses, but got interval guesses for $(join(filter(x -> x[2] isa Tuple, shoot_pars), ", ")).")
    all(ModelingToolkit.isparameter, keys(shoot_pars)) || error("Shooting parameters $(join(filter(!ModelingToolkit.isparameter, keys(shoot_pars)), ", ")) must be declared with @parameters (not @variables). If the shooting parameter represents an initial value, add a trivial parameter and equation for the initial condition.")

    parsk = mergesafe(pars, Dict(par => first(guess) for (par, guess) in shoot_pars)) # if guess is a tuple (x1, x2) for bracketing solvers, then use just x1 for setting up the problem

    pt = pt && k in Set(ModelingToolkit.get_ps(M))
    if pt
        parsk[k] = NaN
    end

    sys = ModelingToolkit.flatten(background(M)) # flatten once, so it can be split and compiled below

    if bg === true
        bg, backwards = split_stages(sys)
    else
        bg isa Tuple && !isempty(bg) || error("bg must be true or a non-empty Tuple of background stages, but got $bg")
        backwards = map(stagebackwards, bg)
    end
    nbg = length(bg)
    iterminate = findfirst(!, backwards) # first forwards stage
    isnothing(terminate) || !isnothing(iterminate) || error("All background stages are integrated backwards, so none can terminate at the event $terminate; pass terminate = nothing.")

    bgopts = (jac = bgjac, sparse = bgsparse, bgopts...) # extra options take precedence
    bgprobs = ODEProblem[]
    splitvars = [] # variables of this and all previous stages
    splvars = Symbolics.SymbolicT[] # unknowns of all previous stages, which are splined into the next ones
    for i in 1:nbg
        append!(splitvars, bg[i])
        stagesys = i == nbg ? sys : split_system(sys, splitvars) # the last stage has the complete background
        if isempty(splvars)
            stagesys = mtkcompile(stagesys)
            stageparsk = parsk
        else
            splvarset = Set(basevar.(splvars))
            stagesys, splpar = mtkcompile_spline(stagesys, splvars; removeics! =ics -> remove_initial_conditions!(ics, splvarset))
            stageparsk = merge(parsk, Dict(splpar => dummyspline(length(splvars)))) # set dummy spline parameter (only known to this stage, so keep parsk clean for the perturbations)
        end
        if debug
            stagesys = debug_system(stagesys)
        end

        # Stop the first forwards stage at the event given by the symbolic `terminate` event (default is today: a ~ 1)
        callback = if isnothing(terminate) || i != iterminate
            nothing # no event; integrate the whole span
        else
            eventfunc = ModelingToolkit.build_explicit_observed_function(stagesys, terminate.lhs - terminate.rhs) # works whether the event's variables are independent, unknown or observed
            ContinuousCallback(
                (u, t, integrator) -> eventfunc(u, integrator.p, t), terminate!;
                save_positions = (true, false), # don't duplicate final point
                rootfind = SciMLBase.RightRootFind # prefer right root, so a(τ₀) ≤ 1.0 and root finding algorithms get different signs also today (alternatively, try to enforce integrator.u[aidx] = 1.0 in affect! and set save_positions = (false, true), although this didn't work exactly last time)
            )
        end

        ts = ModelingToolkit.get_tearing_state(stagesys)
        if !isempty(splvars)
            @set! stagesys.tearing_state = nothing # splining reorders variables and gives an incorrect Jacobian; see comment in the perturbations stage below
        end
        stagetspan = backwards[i] ? reverse(tspan) : tspan # later stages get their span from the previous stages when solved
        stageparsk = i == nbg ? stageparsk : restrict(stageparsk, stagesys) # the split-off stage only knows about some variables
        stage = ODEProblem{iip, specialize}(stagesys, stageparsk, stagetspan; fully_determined, callback, p_constructor, stageopts(bgopts, i, nbg)..., kwargs...)
        if !isempty(splvars)
            newsys = stage.f.sys
            @set! newsys.tearing_state = ts
            @set! stage.f = remake(stage.f; sys = newsys)
        end

        push!(bgprobs, stage)
        append!(splvars, unknowns(stage.f.sys))
    end
    bg = Tuple(bgprobs)

    if pt
        pt = perturbations(M)
        if spline == true
            spline = splvars # all background unknowns
        end
        if spline isa AbstractVector && !isempty(spline)
            pt, splpar = mtkcompile_spline(pt, spline)
            parsk = merge(parsk, Dict(splpar => dummyspline(length(spline)))) # set dummy spline parameter
        else
            pt = mtkcompile(pt)
        end
        if debug
            pt = debug_system(pt)
        end
        ts = ModelingToolkit.get_tearing_state(pt)
        @set! pt.tearing_state = nothing # additional pass in mtkcompile_spline modifies variable ordering and leads to an incorrect Jacobian; reset tearing state to nothing to trigger "manual" computation of the Jacobian
        pt = ODEProblem{iip, specialize}(pt, parsk, tspan; fully_determined, jac = ptjac, sparse = ptsparse, p_constructor, ptopts..., kwargs...)
        # restore tearing state via remake (not @set!) on pt.f while preserving the specialize level
        # (@set!-ing into a nested AbstractSciMLFunction field reconstructs it through ConstructionBase,
        # whose constructorof for SciML function types hardcodes SciMLBase.DEFAULT_SPECIALIZATION (i.e. AutoSpecialize))
        newsys = pt.f.sys
        @set! newsys.tearing_state = ts
        @set! pt.f = remake(pt.f; sys = newsys)
    else
        pt = nothing
    end

    pars = [unwrap(par) for (par, val) in pars]
    shoot_conditions = Equation[eq.lhs - eq.rhs ~ 0 for eq in shoot_conditions]
    return CosmologyProblem(M, bg, pt, pars, shoot_pars, shoot_conditions, terminate)
end

# restrict a variable/parameter map (or list) to those that are part of the problem, which can be a reduced subsystem
problem_symbols(sys::System) = Set(unwrap.([unknowns(sys); parameters(sys)]))
problem_symbols(prob::ODEProblem) = problem_symbols(prob.f.sys)
restrict(x::Missing, _) = x
restrict(x::Dict, target) = (syms = problem_symbols(target); filter(kv -> basevars(first(kv)) ⊆ syms, x))

# Create a function that returns the background stages bg with the parameters pars set to new values in all stages that have them
function bgsetter(bg::Tuple, pars; kwargs...)
    stageis = map(stage -> [i for (i, par) in enumerate(pars) if unwrap(par) in problem_symbols(stage)], bg) # each stage only knows about the parameters it needs
    setters = map((stage, is) -> isempty(is) ? nothing : setsym_oop(stage, pars[is]), bg, stageis)
    return vals -> map(bg, stageis, setters) do stage, is, setter
        isempty(is) && return stage
        u0, p = setter(stage, vals[is])
        return remake(stage; u0, p, kwargs...) # create updated problem (don't overwrite old)
    end
end

"""
    remake_function(prob::CosmologyProblem, pars; kwargs...)

Create an efficient function `f` for updating the values of the independent parameters `pars` in `prob`.
It is called like `newprob = f(vals)`, where `vals` are the new numerical values in the same order as in `pars`.
The symbolic parameters `pars` can be a single parameter or a vector or tuple of parameters.

# Examples
```julia
probf = remake_function(prob, M.c.Ω₀)
newprob = probf(0.3)

probf = remake_function(prob, [M.c.Ω₀, M.g.h])
newprob = probf([0.3, 0.7])
```
"""
function remake_function(prob::CosmologyProblem, pars; kwargs...)
    @unpack bg, pt = prob
    scalar = !(pars isa Union{AbstractArray, Tuple})
    parsvec = scalar ? [pars] : collect(pars)

    updatable = Set(prob.pars)
    nonupdatable = filter(par -> !(unwrap(par) in updatable), parsvec)
    isempty(nonupdatable) || error("Cannot update $(join(nonupdatable, ", ")) because they are not among the independent parameters $(join(updatable, ", ")).")

    bgset = bgsetter(bg, parsvec; kwargs...)
    ptset = isnothing(pt) ? nothing : setsym_oop(pt, parsvec)
    return vals -> begin
        vals = scalar ? [vals] : vals
        bgnew = bgset(vals)
        ptnew = isnothing(pt) ? nothing : remake(pt; zip((:u0, :p), ptset(pt, vals))..., kwargs...)
        return CosmologyProblem(prob.M, bgnew, ptnew, prob.pars, prob.shoot, prob.conditions, prob.terminate)
    end
end

"""
    remake(prob::CosmologyProblem, pars; kwargs...)

Return a new problem with updated independent parameter values `pars` (a `Dict`, pair or vector of pairs in the form `par => val`).
Unspecified parameters keep their values in `prob`.

For repeated updates, prefer [`remake_function`](@ref) instead.

# Examples
```julia
newprob = remake(prob, M.c.Ω₀ => 0.3)
newprob = remake(prob, [M.c.Ω₀ => 0.3, M.g.h => 0.7])
newprob = remake(prob, Dict(M.c.Ω₀ => 0.3, M.g.h => 0.7))
```
"""
function remake(prob::CosmologyProblem, pars; kwargs...)
    pars = pars isa Pair ? [pars] : collect(pars)
    return remake_function(prob, first.(pars); kwargs...)(last.(pars))
end

issparse(M::Nothing) = false
issparse(x) = SparseArrays.issparse(x)
issparse(prob::ODEProblem) = issparse(prob.f.jac_prototype)

function default_bgalg(prob::ODEProblem; stiff = true)
    if issparse(prob)
        linsolve = PureKLUFactorization()
    else
        linsolve = RFLUFactorization(throwerror = true)
    end
    if stiff
        return Rodas5P(; linsolve)
    else
        return Tsit5(; linsolve)
    end
end
default_bgalg(prob::CosmologyProblem; kwargs...) = map(stage -> default_bgalg(stage; kwargs...), prob.bg) # one per stage

function default_ptalg(prob::ODEProblem; accuracy = 2)
    if issparse(prob)
        linsolve = PureKLUFactorization()
    else
        linsolve = RFLUFactorization(throwerror = true)
    end
    nlsolve = NLNewton(fast_convergence_cutoff = 0, κ = 1)
    if accuracy == 0
        return TRBDF2(; linsolve, nlsolve)
    elseif accuracy == 1
        return KenCarp4(; linsolve, nlsolve)
    else
        return Rodas5P(; linsolve) # does not do nonlinear solve
    end
end
default_ptalg(prob::CosmologyProblem; kwargs...) = default_ptalg(prob.pt; kwargs...)
default_ptalg(prob::Nothing) = nothing

function default_shootalg(prob::CosmologyProblem; accuracy = 2)
    if length(prob.shoot) == 1 && only(values(prob.shoot)) isa Tuple
        return ITP() # bracketing solver for interval guesses, regardless of accuracy
    elseif accuracy == 0
        return NewtonRaphson(linesearch = BackTracking())
    else # accuracy >= 1
        return TrustRegion()
    end
end
default_shootalg() = nothing

function check_solve_args(prob::ODEProblem, alg)
    if hasproperty(alg, :linsolve) && !isnothing(alg.linsolve) # if nothing, OrdinaryDiffEq automatically finds a compatible linear solver
        issparse(prob) && !(alg.linsolve isa LinearSolve.AbstractSparseFactorization) && error("ODE with sparse Jacobian must be solved with sparse linear solver")
        !issparse(prob) && alg.linsolve isa LinearSolve.AbstractSparseFactorization && error("ODE with dense Jacobian must be solved with dense linear solver")
    end
end

# TODO: want to use ODESolution's solver-specific interpolator instead of error-prone spline
"""
    solve(
        prob::CosmologyProblem, ks::Union{Nothing, AbstractArray} = nothing;
        bgalg = default_bgalg(prob), bgreltol = 1e-7, bgabstol = 1e-7, bgopts = (),
        ptalg = default_ptalg(prob), ptreltol = 1e-5, ptabstol = 1e-5, ptopts = (),
        shootalg = default_shootalg(prob), shootabstol = 1e-5, shootopts = (),
        thread = true, verbose = false, kwargs...
    )

Solve the cosmological problem `prob` up to the perturbative level with wavenumbers `ks` (or only to the background level if it is empty).
The background stages, perturbations and shooting method are solved with the algorithms `bgalg`, `ptalg` and `shootalg`
and the tolerances `bgreltol`/`bgabstol`, `ptreltol`/`ptabstol` and `shootabstol`.
The extra options `bgopts`, `ptopts` and `shootopts` are passed to the same `solve()` calls, and override the prefixed options above.
The unprefixed options in `kwargs` (like `reltol` or `maxiters`) are applied to both the background and perturbations last, and override all of them.
Each background option can be a single value for all stages, or a Tuple with one value per stage.
If `threads`, integration over independent perturbation modes are parallellized.

See also [`solvebg`](@ref) and [`solvept`](@ref).
"""
function solve(
    prob::CosmologyProblem, ks::Union{Nothing, AbstractArray} = nothing;
    bgalg = default_bgalg(prob), bgreltol = 1e-7, bgabstol = 1e-7, bgopts = (),
    ptalg = default_ptalg(prob), ptreltol = 1e-5, ptabstol = 1e-5, ptopts = (),
    shootalg = default_shootalg(prob), shootabstol = 1e-5, shootopts = (),
    thread = true, verbose = false, kwargs...
)
    bgopts = (alg = bgalg, reltol = bgreltol, abstol = bgabstol, bgopts...)
    ptopts = (alg = ptalg, reltol = ptreltol, abstol = ptabstol, ptopts...)
    shootopts = (alg = shootalg, abstol = shootabstol, shootopts...)

    bgsols = solvebg(prob; shootopts, verbose, bgopts..., kwargs...)

    if isnothing(ks) || isempty(ks) || !all(successful_retcode, bgsols) # no perturbations requested, or they cannot be set up on a failed background
        ks = nothing
        ptsol = nothing
    else
        ptsol = solvept(prob.pt, bgsols, ks; thread, verbose, ptopts..., kwargs...)
    end

    return CosmologySolution(prob, bgsols, ks, ptsol)
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
    solvebg(bgprob::ODEProblem; alg = default_bgalg(bgprob), reltol = 1e-7, abstol = 1e-7, verbose = false, name = "Background", kwargs...)

Solve the background stage `bgprob` and return its solution.
Its splines from previous stages must already be set, for example with [`setupbg`](@ref).
"""
function solvebg(bgprob::ODEProblem; alg = default_bgalg(bgprob), reltol = 1e-7, abstol = 1e-7, verbose = false, name = "Background", kwargs...)
    check_solve_args(bgprob, alg)
    sol = solve(bgprob, alg; verbose = verbosity(verbose), reltol, abstol, kwargs...)
    if !successful_retcode(sol)
        @warn warning_failed_solution(sol, name; verbose)
    end
    return sol
end

"""
    solvebg(bg::Tuple; verbose = false, kwargs...)

Solve the background stages `bg` in order, each set up with the solutions of the previous stages (see [`setupbg`](@ref)), and return a Tuple with their solutions.
If a stage fails, the later stages are not solved.
Each option in `kwargs` can be a single value for all stages, or a Tuple with one value per stage.
"""
function solvebg(bg::Tuple; verbose = false, kwargs...)
    n = length(bg)
    bgsols = ()
    for i in 1:n
        bgprob = setupbg(bg[i], bgsols)
        opts = stageopts(kwargs, i, n)
        if i < n
            opts = (; opts..., save_everystep = true, save_start = true, save_end = true, dense = true) # spline the whole solution into later stages
        end
        bgsol = solvebg(bgprob; verbose, name = "Background stage $i", opts...)
        bgsols = (bgsols..., bgsol)
        successful_retcode(bgsol) || break # cannot set up later stages
    end
    return bgsols
end

# TODO: more generic shooting method that can do anything (e.g. S8)
function _solvebg_shoot_f(x, p)
    n, setvars, getconds, scale, kwargs, verbose, varstrs, constrs = p # unpack
    u = x .* scale
    bgsols = solvebg(setvars(u isa Number ? [u] : u); kwargs..., save_everystep = false, save_start = true, save_end = true, verbose)
    if length(bgsols) < n || !successful_retcode(bgsols[end])
        verbose && eltype(u) <: AbstractFloat && println("Shooting: ODE failed with ", varvalstr(varstrs, u), " (returning NaN)")
        return u .* NaN # return NaN instead of erroring, so solvers can use this information to backtrack/retry into valid regions
    end
    bgsol = bgsols[end] # the complete background
    result = getconds(bgsol, argmax(bgsol.t)) # today is the last or first saved step, depending on direction
    result = u isa Number ? only(result) : result
    verbose && eltype(u) <: AbstractFloat && println("Shooting: ", varvalstr(varstrs, u), " -> ", varvalstr(constrs, result))
    return result
end

# Solve the background stages with the shooting method for the parameters `vars` (mapped to initial guesses), so that the equations `conditions` hold at the final time
function solvebg(bg::Tuple, vars, conditions; shootopts = (alg = default_shootalg(), abstol = 1e-5), verbose = false, kwargs...)
    length(vars) == length(conditions) || error("Different number of shooting parameters and conditions")

    guess = collect(values(vars))
    vars = collect(keys(vars))
    conditions = map(eq -> eq.lhs - eq.rhs, conditions)
    varstrs = string.(vars)
    constrs = string.(conditions)
    guess = map(g -> issymbolic(g) ? bg[end][g] : g, guess) # evaluate symbolic guesses, keep numerical ones as they are
    if length(vars) == 1 # work with scalars instead of vectors to support interval methods
        guess = only(guess)
    end
    scale = guess isa Tuple ? 1 : map(g -> max(abs(g), one(g)), guess) # solve for large parameters relative to their guesses, so they are of order unity
    setvars = bgsetter(bg, vars) # efficient setter
    getconds = getsym(bg[end], conditions) # efficient getter in the complete background

    if guess isa Tuple
        if shootopts.alg isa AbstractBracketingAlgorithm
            NonlinearProblemT = IntervalNonlinearProblem
        else
            error("Shooting with interval guess requires bracketing nonlinear solver")
        end
    else # guess isa Vector
        if shootopts.alg isa AbstractBracketingAlgorithm
            error("Shooting with scalar guesses requires nonbracketing nonlinear solver")
        else
            NonlinearProblemT = NonlinearProblem
        end
    end
    prob = NonlinearProblemT(_solvebg_shoot_f, guess ./ scale, (length(bg), setvars, getconds, scale, kwargs, verbose, varstrs, constrs))
    sol = solve(prob; shootopts...)
    u = sol.u .* scale

    if !successful_retcode(sol)
        error("Shooting failed to converge. Last result was $(varvalstr(varstrs, u)). Run with `verbose = true` for more output. Change the initial shooting guesses.")
    end

    return solvebg(setvars(u isa Number ? [u] : u); verbose, kwargs...)
end

# Set the nonnumeric background spline parameter of a problem (and optionally copy tunable parameters into it)
function setspline(prob::ODEProblem, bgspline, bgtunables = nothing)
    function concretize(p)
        # copy background tunable params into p and concretize spline type
        isnothing(bgtunables) || SciMLStructures.replace!(Tunable(), p, bgtunables)
        if !isempty(p.nonnumeric)
            @set! p.nonnumeric = ([bgspline],)
        end
        return p
    end

    newp = concretize(prob.p)
    @set! prob.p = newp

    # same for initialization problem
    # prevent type assert error, see https://github.com/hersle/SymBoltz.jl/pull/96
    if !isnothing(prob.f.initialization_data)
        newinitp = concretize(prob.f.initialization_data.initializeprob.p)
        new_initprob = remake(prob.f.initialization_data.initializeprob; p = newinitp)
        # remake (not @set!) on prob.f (see comment in CosmologyProblem)
        newid = prob.f.initialization_data
        @set! newid.initializeprob = new_initprob
        @set! prob.f = remake(prob.f; initialization_data = newid)
    end

    return prob, newp
end

"""
    setupbg(bgprob::ODEProblem, bgsols::Tuple)

Prepare the background stage `bgprob` to be solved on top of the solutions `bgsols` of all previous stages:
spline their unknowns into it, and integrate over the span of the last previous stage.
"""
function setupbg(bgprob::ODEProblem, bgsols::Tuple)
    isempty(bgsols) && return bgprob
    tspan = extrema(bgsols[end].t)
    tspan = isbackwards(bgprob) ? reverse(tspan) : tspan
    p = bgprob.p
    if !isempty(p.nonnumeric)
        bgprob, p = setspline(bgprob, spline(bgsols...))
    end
    return remake(bgprob; u0 = bgprob.u0, p, tspan)
end

"""
    solvebg(prob::CosmologyProblem; shootopts = (alg = default_shootalg(prob), abstol = 1e-5), verbose = false, kwargs...)

Solve all background stages of the cosmological problem `prob` in order, and return a Tuple with their solutions.
If the problem requires shooting, all stages are solved repeatedly until the shooting conditions hold at the final time (today) of the last stage.
Each option in `kwargs` can be a single value for all stages, or a Tuple with one value per stage.
"""
function solvebg(prob::CosmologyProblem; shootopts = (alg = default_shootalg(prob), abstol = 1e-5), verbose = false, kwargs...)
    isempty(prob.shoot) && return solvebg(prob.bg; verbose, kwargs...)
    return solvebg(prob.bg, prob.shoot, prob.conditions; shootopts, verbose, kwargs...)
end

function setuppt(ptprob::ODEProblem, bgsols::Tuple)
    tspanbg = extrema(bgsols[end].t) # e.g. until the background terminates
    bgtunables = canonicalize(Tunable(), parameter_values(bgsols[end]))[1] # tunable parameters from the complete background (e.g. set by shooting)

    # copy parameters from background solution to perturbations problem, and spline all background unknowns into it
    ptprob, newp = setspline(ptprob, spline(bgsols...), bgtunables)

    kset! = ModelingToolkit.setp(ptprob, k)
    return k -> begin
        p = copy(newp) # newp specializes on spline types, while ptprob0.p does not; see https://github.com/SciML/ModelingToolkit.jl/issues/3715
        kset!(p, k)
        newptprob = remake(ptprob; u0 = ptprob.u0, p = p, tspan = tspanbg)
        return newptprob
    end
end

"""
    solvept(ptprob::ODEProblem, bgsols::Tuple, ks::AbstractArray; alg = default_ptalg(ptprob), reltol = 1e-5, abstol = 1e-5, output_func = (sol, i) -> sol, thread = true, verbose = false, kwargs...)

Solve the perturbation cosmology problem `ptprob` with wavenumbers `ks` on top of the solutions `bgsols` of all background stages (see [`solvebg`](@ref)).
If `thread` and Julia is running with multiple threads, the solution of independent wavenumbers is parallellized.
The return value is a vector with one `ODESolution` per wavenumber, or its mapping through `output_func` if a custom transformation is passed.
"""
function solvept(ptprob::ODEProblem, bgsols::Tuple, ks::AbstractArray; alg = default_ptalg(ptprob), reltol = 1e-5, abstol = 1e-5, output_func = (sol, i) -> sol, callback = (i -> nothing), thread = true, verbose = false, kwargs...)
    check_solve_args(ptprob, alg)

    #= # do not show threading warnings; these are Julia runtime options that the user is reponsible for setting
    if thread && Threads.nthreads() == 1
        thread = false
        @warn "Multi-threading over perturbation modes was requested, but disabled, since Julia is running with only 1 thread. Restart Julia with more threads (e.g. `julia --threads=auto`) to enable multi-threading, or pass thread = false to explicitly disable it." maxlog=1
    end
    if thread && BLAS.get_num_threads() > 1
        @warn "Multi-threading over perturbation modes was requested, but BLAS is running with $(BLAS.get_num_threads()) threads.\nIt is recommended to restrict BLAS to one thread with `using LinearAlgebra: BLAS; BLAS.set_num_threads(1)`.\nFor more information, see https://docs.julialang.org/en/v1/manual/performance-tips/#man-multithreading-linear-algebra." maxlog=1
    end
    =#

    # TODO: can I exploit that the structure of the perturbation ODEs is ẏ = J * y with "constant" J?
    ptprobf = setuppt(ptprob, bgsols)

    function output_func_warn(sol, i)
        if !successful_retcode(sol)
            @warn warning_failed_solution(sol, "Perturbation (mode k = $(ks[i]))"; verbose)
        elseif verbose
            print("\rSolved perturbations for wavenumber k = $(ks[i])")
        end
        return output_func(sol, i)
    end

    ptsols = fetch.(@spawnif output_func_warn(solve(ptprobf(ks[i]), alg; verbose = verbosity(verbose), reltol, abstol, callback = callback(i), kwargs...), i) thread for i in eachindex(ks)) # wait for all tasks to finish and get the returned solutions
    verbose && println()
    return ptsols
end
"""
    solvept(ptprob::ODEProblem; alg = default_ptalg(ptprob), reltol = 1e-5, abstol = 1e-5, kwargs...)

Solve the perturbation problem `ptprob` and return the solution.
Its wavenumber and background spline must already be initialized, for example with `setuppt`.

# Examples

```julia
# ...
prob = CosmologyProblem(M, pars)
bgsols = solvebg(prob)
ptprobf = SymBoltz.setuppt(prob.pt, bgsols)
k = 1.0
ptprob = ptprobf(k)
ptsol = solvept(ptprob)
```
"""
function solvept(ptprob::ODEProblem; alg = default_ptalg(ptprob), reltol = 1e-5, abstol = 1e-5, kwargs...)
    return solve(ptprob, alg; reltol, abstol, kwargs...)
end

function time_today(prob::CosmologyProblem)
    return maximum(solvebg(prob)[end].t)
end

"""
    issuccess(sol::CosmologySolution)

Returns whether the solution of a cosmological problem was successful (i.e. not failing due to instability or too many time steps).
"""
function issuccess(sol::CosmologySolution)
    bgok = all(successful_retcode, sol.bg)
    ptok = isnothing(sol.pts) || all(successful_retcode, sol.pts)
    return bgok && ptok
end

function integrate(xs, ys; integrator = Trapezoidal())
    return NumericalIntegration.integrate(xs, ys, integrator)
end
integrate_cumulative(sol::CosmologySolution, x, y) = cumul_integrate(sol[x], sol[y])
integrate_cumulative(sol::CosmologySolution, y) = integrate_cumulative(sol, sol.prob.M.τ, y)

# TODO: don't select time points as 2nd/3rd index, since these points will vary
const SymbolicIndex = Union{Num, AbstractArray{Num}}

function Base.getindex(sol::CosmologySolution, i::SymbolicIndex)
    iv = ModelingToolkit.get_iv(sol.prob.M)
    if all(var -> ModelingToolkit.isparameter(var) && !isequal(var, iv), basevars(i)) # expression of only parameters (but don't catch independent variable as parameter)
        return sol.bg[end].ps[i] # assume all parameters are in background # TODO: index sol directly when this is fixed? https://github.com/SciML/ModelingToolkit.jl/issues/3267
    else
        return sol.bg[end](timeseries(sol); idxs = i).u # the last stage has all background variables
    end
end
function Base.getindex(sol::CosmologySolution, i::SymbolicIndex, j)
    return stack(sol[i][j])
end
Base.getindex(sol::CosmologySolution, i::Int, j::SymbolicIndex, k = :) = sol.pts[i][j, k]
Base.getindex(sol::CosmologySolution, i, j::SymbolicIndex, k = :) = [stack(sol[_i, j, k]) for _i in i]
Base.getindex(sol::CosmologySolution, i::Colon, j::SymbolicIndex, k = :) = sol[1:length(sol.pts), j, k]

# TODO: match variable convention (i.e. δ(τ, k))
function (sol::CosmologySolution)(is::AbstractArray, ts::AbstractArray)
    return sol.bg[end](ts; idxs = is)[:, :]
end
(sol::CosmologySolution)(i::Num, ts::AbstractArray) = sol([i], ts)[1, :]
(sol::CosmologySolution)(is::AbstractArray, t::Number) = sol(is, [t])[:, 1]
(sol::CosmologySolution)(i::Num, t::Number) = sol([i], [t])[1, 1]

# similar to https://github.com/SciML/SciMLBase.jl/blob/c568c0eb554ba78440a83792f058073c286a55d3/src/solutions/ode_solutions.jl#L277
function getfunc(sol::ODESolution, var; continuity = :left)
    ps = SciMLBase.parameter_values(sol)
    return t -> begin
        state = SciMLBase.ProblemState(; u = sol.interp(t, nothing, Val{0}, ps, continuity), p = ps, t)
        return getsym(sol, var)(state)
    end
end

function getsym(provider::CosmologyProblem, p)
    getsym_bg = SymbolicIndexingInterface.getsym(provider.bg[end], p)
    return provider -> getsym_bg(provider.bg[end])
end
function getsym(provider::CosmologySolution, p)
    getsym_bg = SymbolicIndexingInterface.getsym(provider.bg[end], p)
    return provider -> getsym_bg(provider.bg[end])
end

function neighboring_modes_indices(sol::CosmologySolution, k)
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

Base.eltype(sol::CosmologySolution) = eltype(sol.bg[end])

function (sol::CosmologySolution)(out::AbstractArray, is::AbstractArray, ts::AbstractArray, ks::AbstractArray; smart = true, ktransform = log)
    if isnothing(sol.ks) || isempty(sol.ks)
        throw(error("No perturbations solved for. Pass ks to solve()."))
    end
    if !issorted(sol.ks)
        throw(error("Solution wavenumbers are not sorted in ascending order"))
    end
    kmin, kmax = extrema(sol.ks)
    minimum(ks) >= kmin || throw("Requested wavenumber k = $(minimum(ks)) is below the minimum solved wavenumber $kmin")
    maximum(ks) <= kmax || throw("Requested wavenumber k = $(maximum(ks)) is above the maximum solved wavenumber $kmax")

    # Pre-allocate intermediate and output arrays
    v = similar(sol.bg[end], length(is), length(ts))
    v1 = similar(sol.bg[end], length(is), length(ts))
    v2 = similar(sol.bg[end], length(is), length(ts))

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
            v1 .= sol.pts[i1](ts; idxs=is) # https://docs.sciml.ai/DiffEqDocs/latest/basics/solution/ # TODO: allocate less or make in-place (https://github.com/SciML/OrdinaryDiffEq.jl/issues/2562)
            i1_prev = i1
        end
        if i2 != i2_prev || !smart
            v2 .= sol.pts[i2](ts; idxs=is) # TODO: getu or similar for speed? possible while preserving interpolation?
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
    out = similar(sol.bg[end], length(is), length(ts), length(ks))
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

"""
    timeseries(bgsols::Tuple; kwargs...)
    timeseries(sol::CosmologySolution; kwargs...)

Return the time steps of all background stage solutions within the span of the last stage, in ascending order.
"""
function timeseries(bgsols::Tuple{Vararg{ODESolution}}; kwargs...)
    ts = sort!(mapreduce(sol -> collect(sol.t), vcat, bgsols)) # collect, so no solution's own time steps are mutated
    tmin, tmax = extrema(bgsols[end].t) # each stage integrates within the span of the previous one
    ts = [t for (i, t) in enumerate(ts) if tmin ≤ t ≤ tmax && (i == 1 || t != ts[i-1])] # remove duplicates with ==, which (unlike unique) ignores ForwardDiff partials
    return timeseries(ts; kwargs...)
end
timeseries(sol::CosmologySolution; kwargs...) = timeseries(sol.bg; kwargs...)
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
    varfunc = getfunc(sol.bg[end], var)
    f(t, p) = varfunc(t) - p # var(t) == val when f(t) == 0
    tspan = extrema(sol.bg[end].t)
    prob = IntervalNonlinearProblem(f, tspan, vals[1]; kwargs...)
    return map(val -> solve(remake(prob; p = val); alg).u, vals)
end
"""
    timeseries(sol::CosmologySolution, var, dvar, vals::AbstractArray)

Find the times when some variable `var` equals some values `vals` with a Hermite spline, also taking the derivative `dvar` of `var` into account.
"""
function timeseries(sol::CosmologySolution, var, dvar, vals::AbstractArray; kwargs...)
    ts = timeseries(sol; kwargs...)
    xs = sol(var, ts)
    ẋs = sol(dvar, ts) # dx/dt
    ṫs = 1 ./ ẋs # dt/dx
    spl = spline(ts, ṫs, xs)
    ts = spl(vals)
    return ts
end

"""
    unknowns(prob::CosmologyProblem)

Get all unknown variables from the background stages and perturbations of the cosmological problem `prob`.
"""
function unknowns(prob::CosmologyProblem)
    bg = reduce(vcat, (unknowns(stage.f.sys) for stage in prob.bg))
    return [bg; unknowns(prob.pt.f.sys)]
end

"""
    parameters(prob::CosmologyProblem; nonnumeric = false)

Get all parameter values of the cosmological problem `prob`.
"""
function parameters(prob::CosmologyProblem; nonnumeric = false)
    ptpars = isnothing(prob.pt) ? Dict() : parameters(prob.pt)
    pars = merge(map(parameters, prob.bg)..., ptpars)
    !nonnumeric && filter!(par_val -> par_val[2] isa Number, pars)
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

# Statistics for solutions of the background stages
function statsbg(sol::CosmologySolution)
    return map(bgsol -> bgsol.stats, sol.bg)
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
