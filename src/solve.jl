import Base: nameof
import LinearAlgebra: issuccess
import CommonSolve: solve
import SciMLBase: remake, successful_retcode
import PreallocationTools: DiffCache, get_tmp
import SciMLStructures
import SciMLStructures: canonicalize, Tunable
import OhMyThreads: TaskLocalValue
import SymbolicIndexingInterface: getsym, setp, parameter_values

background(sys) = transform((sys, _) -> taylor(sys, ϵ, 0:0; fold = false), sys)
perturbations(sys) = transform((sys, _) -> taylor(sys, ϵ, 0:1; fold = false), sys)

struct CosmologyProblem
    M::ODESystem

    bg::ODEProblem
    pt::Union{ODEProblem, Nothing}

    pars::Dict

    shoot::Base.KeySet
    conditions::AbstractArray

    var2spl::Dict
end

struct CosmologySolution{Tbg, Tks, Tpts, Th}
    prob::CosmologyProblem
    bg::Tbg
    ks::Tks
    pts::Tpts
    h::Th
end

algname(alg) = string(nameof(typeof(alg)))
algname(alg::CompositeAlgorithm) = join(algname.(alg.algs), "+")

function Base.show(io::IO, prob::CosmologyProblem; indent = "  ")
    print(io, "Cosmology problem for model ")
    printstyled(io, nameof(prob.bg.f.sys), '\n'; bold = true)

    printstyled(io, "Stages:"; bold = true)
    if !isnothing(prob.bg)
        print(io, '\n', indent, "Background")
        print(io, ": ", length(unknowns(prob.bg.f.sys)), " unknowns")
        print(io, ", ", 0, " splines")
    end
    if !isnothing(prob.pt)
        print(io, '\n', indent, "Perturbations")
        print(io, ": ", length(unknowns(prob.pt.f.sys)), " unknowns")
        print(io, ", ", length(prob.var2spl), " splines")
    end

    printstyled(io, "\nParameters:"; bold = true)
    for (par, val) in prob.pars
        print(io, '\n', indent, "$par = $val", par in prob.shoot ? " (shooting guess)" : "")
    end

    !isempty(prob.conditions) && printstyled(io, "\nFinal conditions:"; bold = true) # TODO: initial vs final conditions
    for condition in prob.conditions
        print(io, '\n', indent, condition)
    end
end

function Base.show(io::IO, sol::CosmologySolution; indent = "  ")
    print(io, "Cosmology solution for model ")
    printstyled(io, nameof(sol.prob.bg.f.sys), '\n'; bold = true)

    printstyled(io, "Stages:"; bold = true)
    if !isnothing(sol.bg)
        print(io, '\n', indent, "Background: solved with $(algname(sol.bg.alg)), $(length(sol.bg)) points")
    end
    if !isnothing(sol.pts)
        kmin, kmax = extrema(sol.ks)
        nmin, nmax = extrema(map(length, sol.pts))
        n = length(sol.pts)
        print(io, '\n', indent, "Perturbations: solved with $(algname(sol.pts[1].alg)), $nmin-$nmax points, x$n k ∈ [$kmin, $kmax] H₀/c (linear interpolation between log(k))")
    end

    printstyled(io, "\nParameters:"; bold = true)
    for (par, val) in sol.prob.pars
        print(io, '\n', indent, "$par = $val", par in sol.prob.shoot ? " (shooting solution)" : "")
    end
end

# Split parameters into DifferentialEquations' u0 and p convention
function split_vars_pars(M::ODESystem, x::Dict)
    pars = intersect(keys(x), parameters(M)) .|> ModelingToolkit.wrap # separate parameters from initial conditions # TODO: remove wrap
    vars = setdiff(keys(x), pars) # assume the rest are variables (do it without intersection to capture derivatives initial conditions)
    pars = Dict(par => x[par] for par in pars) # like p
    vars = Dict(var => x[var] for var in vars) # like u0
    return vars, pars
end

"""
    CosmologyProblem(
        M::ODESystem, pars::Dict, shoot_pars = Dict(), shoot_conditions = [];
        bg = true, pt = true, spline = true, debug = false, kwargs...
    )

Create a numerical cosmological problem from the model `M` with parameters `pars`.
Optionally, the shooting method determines the parameters `shoot_pars` (mapped to initial guesses) such that the equations `shoot_conditions` are satisfied at the final time.

If `bg` and `pt`, the model is split into the background and perturbations stages.
If `spline` is a `Bool`, it decides whether all background unknowns in the perturbations system are replaced by splines.
If `spline` is a `Vector`, it rather decides which (unknown and observed) variables are splined.
"""
function CosmologyProblem(
    M::ODESystem, pars::Dict, shoot_pars = Dict(), shoot_conditions = [];
    ivspan = (0.0, 100.0), bg = true, pt = true, spline = true, debug = false, fully_determined = true, kwargs...
)
    pars_full = merge(pars, shoot_pars) # save full dictionary for constructor
    vars, pars = split_vars_pars(M, pars_full)
    parsk = merge(pars, Dict(M.k => NaN)) # k is unused, but must be set
    shoot_pars = keys(shoot_pars)

    if bg
        bg = background(M)
        bg = structural_simplify(bg)
        if debug
            bg = debug_system(bg)
        end
        bg = ODEProblem(bg, vars, ivspan, parsk; fully_determined, kwargs...)
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
        pt, var2spl = structural_simplify_spline(pt, spline)
        if debug
            pt = debug_system(pt)
        end
        vars = remove_initial_conditions!(vars, keys(var2spl)) # must remove ICs of splined variables to avoid overdetermined initialization system
        pt = ODEProblem(pt, vars, ivspan, parsk; fully_determined, kwargs...)
    else
        pt = nothing
        var2spl = Dict()
    end

    return CosmologyProblem(M, bg, pt, pars_full, shoot_pars, shoot_conditions, var2spl)
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
    pars_full = merge(prob.pars, pars) # save full dictionary for constructor
    vars, pars = split_vars_pars(prob.M, pars)
    vars = isempty(vars) ? missing : vars
    pars = isempty(pars) ? missing : pars
    bg = bg && !isnothing(prob.bg) ? remake(prob.bg; u0 = vars, p = pars, build_initializeprob = Val{!isnothing(prob.bg.f.initialization_data)}, kwargs...) : nothing
    if !ismissing(vars)
        vars = remove_initial_conditions!(vars, keys(prob.var2spl)) # must filter ICs in remake, too
    end
    pt = pt && !isnothing(prob.pt) ? remake(prob.pt; u0 = vars, p = pars, build_initializeprob = Val{!isnothing(prob.pt.f.initialization_data)}, kwargs...) : nothing
    shoot_pars = shoot ? prob.shoot : keys(Dict())
    shoot_conditions = shoot ? prob.conditions : []
    return CosmologyProblem(prob.M, bg, pt, pars_full, shoot_pars, shoot_conditions, prob.var2spl)
end

function parameter_updater(prob::CosmologyProblem, idxs; kwargs...)
    # define a closure based on https://docs.sciml.ai/ModelingToolkit/dev/examples/remake/#replace-and-remake # hide
    # TODO: perturbations, remove M, pars, etc. for efficiency?

    @unpack M, bg, pars, shoot, conditions, var2spl = prob
    bgps = parameter_values(bg)
    bgsetp! = setp(bg, idxs)
    diffcache = DiffCache(copy(canonicalize(Tunable(), bgps)[1])) # TODO: separate bg/pt?

    return p -> begin
        buffer = get_tmp(diffcache, p) # get newly typed buffer
        copyto!(buffer, canonicalize(Tunable(), bgps)[1]) # copy all parameters to buffer
        bgps_new = SciMLStructures.replace(Tunable(), bgps, buffer) # get newly typed parameter object
        bgsetp!(bgps_new, p) # set new parameters
        bg_new = remake(bg; p = bgps_new, kwargs...) # create updated problem
        return CosmologyProblem(M, bg_new, nothing, pars, shoot, conditions, var2spl) # TODO: update pars? or remove that field and read parameters from subproblems?
    end
end

# TODO: want to use ODESolution's solver-specific interpolator instead of error-prone spline
"""
    function solve(
        prob::CosmologyProblem, ks::AbstractArray = [];
        aterm = 1.0,
        bgopts = (alg = Rodas4P(), reltol = 1e-8,),
        ptopts = (alg = KenCarp4(), reltol = 1e-8,),
        shootopts = (alg = NewtonRaphson(), reltol = 1e-3, pt = false),
        thread = true, verbose = false, kwargs...
    )

Solve the cosmological problem `prob` up to the perturbative level with wavenumbers `ks` (or only to the background level if it is empty).
The options `bgopts` and `ptopts` are passed to the background and perturbations ODE `solve()` calls,
and `shootopts` to the shooting method nonlinear `solve()`.
If `threads`, integration over independent perturbation modes are parallellized.
"""
function solve(
    prob::CosmologyProblem, ks::AbstractArray;
    term = prob.M.g.a => 1.0,
    bgopts = (alg = Rodas4P(), reltol = 1e-8,),
    ptopts = (alg = KenCarp4(), reltol = 1e-8,),
    shootopts = (alg = NewtonRaphson(), reltol = 1e-3, pt = false),
    thread = true, verbose = false, kwargs...
)
    if !isempty(prob.shoot)
        length(prob.shoot) == length(prob.conditions) || error("Different number of shooting parameters and conditions")
        pars = Dict(par => prob.pars[par] for par in prob.shoot)
        pars = shoot(prob, pars, prob.conditions; verbose, shootopts...)
        prob = remake(prob, pars)
    end

    M = prob.M
    callback = isnothing(term) ? nothing : callback_terminator(prob.bg, term[1], term[2])

    if !isnothing(prob.bg)
        bgprob = prob.bg
        if !haskey(bgopts, :alg)
            bgopts = merge(bgopts, (alg = Rodas4P(),)) # default if unspecified
        end
        bg = solve(bgprob; callback, verbose, bgopts..., kwargs...)

        # Offset optical depth, so it's 0 today
        if have(M, :b)
            idx_τ = ModelingToolkit.variable_index(prob.bg, M.b.rec.τ)
            for i in 1:length(bg.u)
                bg.u[i][idx_τ] -= bg.u[end][idx_τ]
            end
        end
    else
        bg = nothing
    end

    ivspan = extrema(bg.t)

    h = prob.pars[prob.bg.f.sys.g.h]
    if !isnothing(prob.pt) && !isempty(ks)
        ks = k_dimensionless.(ks, h)
        !issorted(ks) && throw(error("ks = $ks are not sorted in ascending order"))

        if Threads.nthreads() == 1 && thread
            @warn "Multi-threading was requested, but disabled, since Julia is running with only 1 thread. Restart Julia with more threads (e.g. `julia --threads=auto`) to enable multi-threading, or pass thread = false to explicitly disable it."
            thread = false
        end

        # TODO: can I exploit that the structure of the perturbation ODEs is ẏ = J * y with "constant" J?
        ptprob0 = remake(prob.pt; tspan = ivspan)
        update_vars = Dict()
        is_unknown = Set(nameof.(Symbolics.operation.(unknowns(prob.bg.f.sys)))) # set for more efficient lookup # TODO: process in CosmologyProblem
        for (var, spl) in prob.var2spl
            varname = nameof(operation(var))
            if varname in is_unknown
                dvar = nothing # compute derivative from ODE f
            else
                dvarname = Symbol(varname, Symbol("̇")) # add \dot # TODO: generalize to e.g. expand_derivatives(D(var))?
                dvar = only(@variables($dvarname(t)))
            end
            splval = spline(bg, var, dvar)
            update_vars = merge(update_vars, Dict(spl => splval))
            verbose && println("Splining $var with $(length(splval.t)) points")
        end
        update = ModelingToolkit.setsym_oop(ptprob0, collect(keys(update_vars)))
        newu0, newp = update(ptprob0, collect(values(update_vars)))
        ptprob0 = remake(ptprob0; tspan = ivspan, u0 = newu0, p = newp)

        kidx = ModelingToolkit.parameter_index(ptprob0, k)
        ptprob_tlv = TaskLocalValue{ODEProblem}(() -> remake(ptprob0; u0 = copy(ptprob0.u0) #= p is copied below =#)) # prevent conflicts where different tasks modify same problem: https://discourse.julialang.org/t/solving-ensembleproblem-efficiently-for-large-systems-memory-issues/116146/11 (alternatively copy just p and u0: https://github.com/SciML/ModelingToolkit.jl/issues/3056)
        nmode, nmodes = Atomic{Int}(0), length(ks)
        ptprobs = EnsembleProblem(; safetycopy = false, prob = ptprob0, prob_func = (_, i, _) -> begin
            ptprob = ptprob_tlv[]
            p = copy(ptprob0.p) # see https://github.com/SciML/ModelingToolkit.jl/issues/3346 and https://github.com/SciML/ModelingToolkit.jl/issues/3056
            setindex!(p, ks[i], kidx)
            return Setfield.@set ptprob.p = p
        end, output_func = (sol, i) -> begin
            if verbose
                atomic_add!(nmode, 1)
                progress = Int(round(nmode[] / nmodes * 100))
                print("\rSolved perturbation modes $(nmode[])/$(nmodes) = $(round(progress)) %")
            end
            return sol, false
        end)
        ensemblealg = thread ? EnsembleThreads() : EnsembleSerial()
        if !haskey(ptopts, :alg)
            ptopts = merge(ptopts, (alg = KenCarp4(),)) # default if unspecified
        end
        ptsols = solve(ptprobs; ensemblealg, trajectories = length(ks), verbose, ptopts..., kwargs...) # TODO: test GPU parallellization
        verbose && println() # end line in output_func
    else
        ptsols = nothing
    end

    return CosmologySolution{typeof(bg), typeof(ks), typeof(ptsols), typeof(h)}(prob, bg, ks, ptsols, h)
end
function solve(prob::CosmologyProblem; kwargs...)
    return solve(prob, []; kwargs...)
end
function solve(prob::CosmologyProblem, k::Number; kwargs...)
    return solve(prob, [k]; kwargs...)
end

"""
    issuccess(sol::CosmologySolution)

Returns whether the solution of a cosmological problem was successful (i.e. not failing due to instability or too many time steps).
"""
function issuccess(sol::CosmologySolution)
    return successful_retcode(sol.bg) && (isnothing(sol.pts) || all(successful_retcode(pt) for pt in sol.pts))
end

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

function (sol::CosmologySolution)(ts::AbstractArray, is::AbstractArray)
    tmin, tmax = extrema(sol.bg.t[[begin, end]])
    #minimum(ts) >= tmin || minimum(ts) ≈ tmin || throw("Requested time t = $(minimum(ts)) is before initial time $tmin")
    #maximum(ts) <= tmax || maximum(ts) ≈ tmax || throw("Requested time t = $(maximum(ts)) is after final time $tmax")
    return permutedims(sol.bg(ts, idxs=is)[:, :])
end
(sol::CosmologySolution)(ts::AbstractArray, i::Num) = sol(ts, [i])[:, 1]
(sol::CosmologySolution)(t::Number, is::AbstractArray) = sol([t], is)[1, :]
(sol::CosmologySolution)(t::Number, i::Num) = sol([t], [i])[1, 1]

# similar to https://github.com/SciML/SciMLBase.jl/blob/c568c0eb554ba78440a83792f058073c286a55d3/src/solutions/ode_solutions.jl#L277
function getfunc(sol::ODESolution, var; deriv = Val{0}, continuity = :left)
    ps = SciMLBase.parameter_values(sol)
    return t -> begin
        state = SciMLBase.ProblemState(; u = sol.interp(t, nothing, deriv, ps, continuity), p = ps, t)
        return getsym(sol, var)(state)
    end
end

function neighboring_modes_indices(sol::CosmologySolution, k)
    k = k_dimensionless.(k, sol.h)
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

function (sol::CosmologySolution)(out::AbstractArray, ks::AbstractArray, ts::AbstractArray, is::AbstractArray; smart = true)
    ks = k_dimensionless.(ks, sol.h)
    isempty(sol.ks) && throw(error("No perturbations solved for. Pass ks to solve()."))
    kmin, kmax = extrema(sol.ks)
    minimum(ks) >= kmin || throw("Requested wavenumber k = $(minimum(ks)) is below the minimum solved wavenumber $kmin")
    maximum(ks) <= kmax || throw("Requested wavenumber k = $(maximum(ks)) is above the maximum solved wavenumber $kmax")
    # disabled; sensitive to exact endpoint values
    #tmin, tmax = extrema(sol.th.t)
    #minimum(ts) >= tmin || throw("Requested time t = $(minimum(ts)) is below minimum solved time $tmin")
    #maximum(ts) <= tmax || throw("Requested time t = $(maximum(ts)) is above maximum solved time $tmin")

    # Pre-allocate intermediate and output arrays
    v = zeros(eltype(sol), (length(is), length(ts)))
    v1 = zeros(eltype(sol), (length(is), length(ts)))
    v2 = zeros(eltype(sol), (length(is), length(ts)))

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
            v1 .= sol.pts[i1](ts; idxs=is)[:, :] # TODO: make in-place (https://github.com/SciML/OrdinaryDiffEq.jl/issues/2562)
            i1_prev = i1
        end
        if i2 != i2_prev || !smart
            v2 .= sol.pts[i2](ts; idxs=is)[:, :] # TODO: getu or similar for speed? possible while preserving interpolation?
            i2_prev = i2
        end
        v .= v1
        if i1 != i2
            # interpolate between solutions
            k1 = sol.ks[i1]
            k2 = sol.ks[i2]
            w = (log(k) - log(k1)) / (log(k2) - log(k1)) # quantities vary smoothest (?) when interpolated in log(k) # TODO: cubic spline?
            @. v += (v2 - v1) * w # add to v1 from above
        end
        for ii in eachindex(is)
            out[ik, :, ii] .= v[ii, :]
        end
    end

    return out
end
function (sol::CosmologySolution)(ks::AbstractArray, ts::AbstractArray, is::AbstractArray; kwargs...)
    out = zeros(eltype(sol), length(ks), length(ts), length(is))
    return sol(out, ks, ts, is; kwargs...)
end
(sol::CosmologySolution)(k::Number, ts::AbstractArray, is::AbstractArray) = sol([k], ts, is)[1, :, :]
(sol::CosmologySolution)(ks::AbstractArray, t::Number, is::AbstractArray) = sol(ks, [t], is)[:, 1, :]
(sol::CosmologySolution)(ks::AbstractArray, ts::AbstractArray, i::Num) = sol(ks, ts, [i])[:, :, 1]
(sol::CosmologySolution)(k::Number, t::Number, is::AbstractArray) = sol([k], [t], is)[1, 1, :]
(sol::CosmologySolution)(k::Number, ts::AbstractArray, i::Num) = sol([k], ts, [i])[1, :, 1]
(sol::CosmologySolution)(ks::AbstractArray, t::Number, i::Num) = sol(ks, [t], [i])[:, 1, 1]
(sol::CosmologySolution)(k::Number, t::Number, i::Num) = sol([k], [t], [i])[1, 1, 1]

function (sol::CosmologySolution)(tmap::Pair, is)
    tvar, ts = tmap
    ts = timeseries(sol, tvar, ts)
    return sol(ts, is)
end

function (sol::CosmologySolution)(ks, tmap::Pair, is)
    tvar, ts = tmap
    ts = timeseries(sol, tvar, ts)
    return sol(ks, ts, is)
end

function timeseries(sol::CosmologySolution; kwargs...)
    ts = sol[t]
    return timeseries(ts; kwargs...)
end
function timeseries(sol::CosmologySolution, k; kwargs...)
    i1, i2 = neighboring_modes_indices(sol, k)
    t1s = sol[i1, t]
    t2s = sol[i2, t]
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
    ivspan = extrema(sol[t])
    prob = IntervalNonlinearProblem(f, ivspan, vals[1]; kwargs...)
    return map(val -> solve(remake(prob; p = val); alg).u, vals)
end
"""
    timeseries(sol::CosmologySolution, var, dvar, vals::AbstractArray)

Find the times when some variable `var` equals some values `vals` with a Hermite spline, also taking the derivative `dvar` of `var` into account.
"""
function timeseries(sol::CosmologySolution, var, dvar, vals::AbstractArray; kwargs...)
    ts = timeseries(sol; kwargs...)
    xs = sol(ts, var)
    ẋs = sol(ts, dvar) # dx/dt
    ṫs = 1 ./ ẋs # dt/dx
    spl = spline(ts, ṫs, xs)
    ts = spl(vals)
    return ts
end

"""
    time_today(sol::CosmologySolution; atoday = 1.0)

Return the time today, when the scale factor `a` equals `atoday`.
"""
function time_today(sol::CosmologySolution; atoday = 1.0)
    M = sol.prob.M
    afinal = sol[M.g.a][end]
    tfinal = sol[M.t][end]
    if sol[M.g.a][end] <= 1.0 && afinal ≈ atoday
        return tfinal # avoid root finding if signs are not different
    else
        return timeseries(sol, M.g.a, atoday)
    end
end

# TODO: more generic version that can do anything (e.g. S8)
# TODO: add dispatch with pars::AbstractArray that uses values in prob for initial guess
"""
    function shoot(
        prob::CosmologyProblem, pars::Dict, conditions::AbstractArray;
        alg = TrustRegion(), bg = true, pt = false, verbose = false, kwargs...
    )

Solve `prob` repeatedly to determine the values in `pars` (with initial guesses) so that the `conditions` are satisfied at the final time.
"""
function shoot(
    prob::CosmologyProblem, pars::Dict, conditions::AbstractArray;
    alg = TrustRegion(), bg = true, pt = false, verbose = false, kwargs...
)
    funcs = [ModelingToolkit.wrap(eq.lhs - eq.rhs) for eq in conditions] # expressions that should be 0 # TODO: shouldn't have to wrap

    function f(vals, _)
        newprob = remake(prob, Dict(keys(pars) .=> vals); bg, pt, shoot = false)
        sol = solve(newprob)
        return sol[funcs, :][:, end] # evaluate all expressions at final time (e.g. today)
    end

    nguess = collect(values(pars))
    nprob = NonlinearProblem(f, nguess)
    nsol = solve(nprob, alg; show_trace = Val(verbose), kwargs...)
    return Dict(keys(pars) .=> nsol.u)
end

"""
    parameters(prob::CosmologyProblem; bg = true, pt = true, spline = false)

Get all parameter values of the cosmological problem `prob`.
"""
function parameters(prob::CosmologyProblem; bg = true, pt = true, spline = false)
    bg = bg && !isnothing(prob.bg) ? parameters(prob.bg) : Dict()
    pt = pt && !isnothing(prob.pt) ? parameters(prob.pt) : Dict()
    pars = merge(bg, pt)
    !spline && delete!.(Ref(pars), values(prob.var2spl))
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
Base.broadcastable(sys::ODESystem) = Ref(sys)
Base.broadcastable(sol::CosmologySolution) = Ref(sol)
