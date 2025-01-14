import Base: nameof
import CommonSolve: solve
import SymbolicIndexingInterface: setp # all_variable_symbols, getname
import OhMyThreads: TaskLocalValue

function background(sys)
    sys = thermodynamics(sys)
    if :b in ModelingToolkit.get_name.(ModelingToolkit.get_systems(sys))
        sys = replace(sys, sys.b.rec => ODESystem([], t; name = :rec)) # remove recombination
    end
    return sys
end

function thermodynamics(sys)
    return transform((sys, _) -> taylor(sys, ϵ, [0]), sys)
end

function perturbations(sys; spline_thermo = true)
    if :b in ModelingToolkit.get_name.(ModelingToolkit.get_systems(sys)) && spline_thermo && !all([eq.rhs === 0 for eq in equations(sys.b.rec)])
        @named rec = thermodynamics_recombination_splined()
        sys = replace(sys, sys.b.rec => rec) # substitute in splined recombination
    end
    return transform((sys, _) -> taylor(sys, ϵ, 0:1), sys)
end

struct CosmologyModel
    sys::ODESystem

    bg::ODESystem
    th::ODESystem
    pt::ODESystem

    spline_thermo::Bool # whether to spline thermodynamics in perturbations
end

struct CosmologyProblem
    pars::Dict
    bg::Union{ODEProblem, Nothing}
    th::Union{ODEProblem, Nothing}
    pt::Union{ODEProblem, Nothing}
end

function CosmologyModel(sys::ODESystem; spline_thermo = true, debug = false)
    # TODO: move to CosmologyProblem; keep M as a simple ODESystem
    bg = structural_simplify(background(sys))
    th = structural_simplify(thermodynamics(sys))
    pt = structural_simplify(perturbations(sys; spline_thermo))

    if debug
        bg = debug_system(bg)
        th = debug_system(th)
        pt = debug_system(pt)
    end

    sys = complete(sys; flatten = false)
    return CosmologyModel(sys, bg, th, pt, spline_thermo)
end

# Forward property access to full system
Base.propertynames(M::CosmologyModel) = propertynames(getfield(M, :sys))
function Base.getproperty(M::CosmologyModel, prop::Symbol)
    if prop in propertynames(M)
        return Base.getproperty(getfield(M, :sys), prop)
    else
        return getfield(M, prop) # hidden access to other fields
    end
end

# Forward inspection functions to full system
nameof(M::CosmologyModel) = nameof(M.sys)
equations(M::CosmologyModel) = equations(M.sys)
observed(M::CosmologyModel) = observed(M.sys)
unknowns(M::CosmologyModel) = unknowns(M.sys)
parameters(M::CosmologyModel) = parameters(M.sys)
initialization_equations(M::CosmologyModel) = initialization_equations(M.sys)
defaults(M::CosmologyModel) = defaults(M.sys)
hierarchy(M::CosmologyModel; describe=true, kwargs...) = hierarchy(M.sys; describe, kwargs...)
Base.show(io::IO, mime::MIME"text/plain", M::CosmologyModel) = show(io, mime, M.sys) # chop off last excessive newline

struct CosmologySolution
    prob::CosmologyProblem
    bg::ODESolution
    th::ODESolution
    ks::AbstractArray
    pts::Union{EnsembleSolution, Nothing}
end

solvername(alg) = string(nameof(typeof(alg)))
solvername(alg::CompositeAlgorithm) = join(solvername.(alg.algs), "+")

function Base.show(io::IO, sol::CosmologySolution)
    print(io, "Cosmology solution for model ")
    printstyled(io, nameof(sol.prob.bg.f.sys), '\n'; bold = true)

    printstyled(io, "Parameters:\n"; bold = true)
    for (par, val) in sol.prob.pars
        print(io, "  $par = $val\n")
    end

    printstyled(io, "Stages:"; bold = true)
    print(io, "\n  1. background: solved with $(solvername(sol.bg.alg)), $(length(sol.bg)) points")
    if !isnothing(sol.th)
        print(io, "\n  2. thermodynamics: solved with $(solvername(sol.th.alg)), $(length(sol.th)) points")
    end
    if !isnothing(sol.pts)
        kmin, kmax = extrema(sol.ks)
        nmin, nmax = extrema(map(length, sol.pts))
        n = length(sol.pts)
        print(io, "\n  3. perturbations: solved with $(solvername(sol.pts[1].alg)), $nmin-$nmax points, x$n k ∈ [$kmin, $kmax] H₀/c (linear interpolation between log(k))")
    end
end

function CosmologyProblem(M::CosmologyModel, pars; aterm = 1.0, thermo = true, debug_initialization = false, guesses = Dict(), jac = false, sparse = false, kwargs...)
    # Split parameters into DifferentialEquations' "u0" and "p" convention
    params = merge(pars, Dict(M.k => NaN)) # k is unused, but must be set
    pars = intersect(keys(params), parameters(M)) # separate parameters from initial conditions
    vars = setdiff(keys(params), pars) # assume the rest are variables (do it without intersection to capture derivatives initial conditions)
    vars = Dict(var => params[var] for var in vars) # like u0
    pars = Dict(par => params[par] for par in pars) # like p

    # TODO: forwards only (makes more sense, can go beyond a=1, etc.)
    tspan = (1e-5, 4.0)
    # TODO: solve only one of bg,th (based on thermo bool)?

    bg = ODEProblem(M.bg, vars, tspan, pars; guesses, fully_determined = true, jac, sparse)
    th = ODEProblem(M.th, vars, tspan, pars; guesses, fully_determined = true, jac, sparse)
    pt = ODEProblem(M.pt, vars, tspan, pars; guesses, fully_determined = true, jac, sparse)

    delete!(params, k) # remove k dummy
    return CosmologyProblem(params, bg, th, pt)
end

# TODO: add generic function spline(sys::ODESystem, how_to_spline_different_vars) that splines the unknowns of a simplified ODESystem 
# TODO: use CommonSolve.step! to iterate background -> thermodynamics -> perturbations?
# TODO: solve thermodynamics only if parameters contain thermodynamics parameters?
# TODO: shoot to reach E = 1 today when integrating forwards
# TODO: want to use ODESolution's solver-specific interpolator instead of error-prone spline
"""
    solve(M::CosmologyModel, pars; aterm = 1.0, solver = Rodas4P(), reltol = 1e-10, kwargs...)

Solve `CosmologyModel` with parameters `pars` at the background level.
"""
function solve(prob::CosmologyProblem; aterm = 1.0, solver = Rodas4P(), reltol = 1e-10, debug_initialization = false, guesses = Dict(), jac = false, sparse = false, kwargs...)
    M = prob.bg.f.sys
    callback = callback_terminator(prob.bg, M.g.a, aterm)
    bg = solve(prob.bg, solver; callback, reltol, kwargs...)
    check_solution(bg.retcode)

    # Remove duplicate endpoints
    if bg.t[end] == bg.t[end-1]
        ilast = length(bg.t)
        deleteat!(bg.t, ilast)
        deleteat!(bg.u, ilast)
    end

    M = prob.th.f.sys
    th = solve(prob.th, solver; callback, reltol, kwargs...)
    check_solution(th.retcode)

    # Offset optical depth, so it's 0 today
    if have(M, :b)
        idx_τ = variable_index(prob.th, M.b.rec.τ)
        for i in 1:length(th.u)
            th.u[i][idx_τ] -= th.u[end][idx_τ]
        end
    end

    return CosmologySolution(prob, bg, th, [], nothing)
end

"""
    solve(M::CosmologyModel, pars, ks; aterm = 1.0, solver = KenCarp4(), reltol = 1e-8, verbose = false, thread = true, jac = false, sparse = false, kwargs...)

Solve `CosmologyModel` with parameters `pars` up to the perturbative level for wavenumbers `ks`.
"""
function solve(prob::CosmologyProblem, ks::AbstractArray; aterm = 1.0, solver = KenCarp4(), reltol = 1e-8, reltol_bg = 1e-10, verbose = false, thread = true, jac = false, sparse = false, kwargs...)
    ks = k_dimensionless.(ks, prob.pars[prob.bg.f.sys.g.h])
    !issorted(ks) && throw(error("ks = $ks are not sorted in ascending order"))
    
    if Threads.nthreads() == 1 && thread
        @warn "Multi-threading was requested, but disabled, since Julia is running with only 1 thread. Restart Julia with more threads (e.g. `julia --threads=auto`) to enable multi-threading, or pass thread = false to explicitly disable it."
        thread = false
    end

    sol = solve(prob; aterm, jac, sparse, reltol = reltol_bg, kwargs...)

    # If the thermodynamics solution should be splined,
    # solve it again (if needed) and update the spline parameters
    M = prob.pt.f.sys # TODO: can I exploit that the structure of the perturbation ODEs is ẏ = J * y with "constant" J?
    if true # || # TODO: M.spline_thermo
        ode_prob0 = prob.pt
        aini = sol[M.g.a][begin]
        τspline = spline(sol[M.b.rec.τ], sol[M.t]) # TODO: when solving thermo with low reltol: even though the solution is correct, just taking its points for splining can be insufficient. should increase number of points, so it won't mess up the perturbations
        τ̇spline = spline(sol[M.b.rec.τ̇], sol[M.t])
        cₛ²spline = spline(sol[M.b.rec.cₛ²], sol[M.t])
        splineset = ModelingToolkit.setsym_oop(ode_prob0, [M.g.a, M.b.rec.τspline, M.b.rec.τ̇spline, M.b.rec.cₛ²spline])
        newu0, newp = splineset(ode_prob0, [aini, τspline, τ̇spline, cₛ²spline])
        ode_prob0 = remake(ode_prob0, u0 = newu0, p = newp)
    end

    kset! = setp(prob.pt, k) # function that sets k on a problem
    ode_prob_tlv = TaskLocalValue{ODEProblem}(() -> deepcopy(ode_prob0)) # prevent conflicts where different tasks modify same problem: https://discourse.julialang.org/t/solving-ensembleproblem-efficiently-for-large-systems-memory-issues/116146/11 (alternatively copy just p and u0: https://github.com/SciML/ModelingToolkit.jl/issues/3056)
    ode_probs = EnsembleProblem(; safetycopy = false, prob = ode_prob0, prob_func = (_, i, _) -> begin
        ode_prob = ode_prob_tlv[]
        verbose && println("$i/$(length(ks)) k = $(ks[i]*k0) Mpc/h")
        kset!(ode_prob, ks[i])
        return ode_prob
    end)
    alg = thread ? EnsembleThreads() : EnsembleSerial()
    ode_sols = solve(ode_probs, solver, alg, trajectories = length(ks); reltol, kwargs...) # TODO: test GPU parallellization
    for i in 1:length(ode_sols)
        check_solution(ode_sols[i].retcode)
    end
    return CosmologySolution(prob, sol.bg, sol.th, ks, ode_sols)
end
function solve(M::CosmologyModel, pars, k::Number; kwargs...)
    return solve(M, pars, [k]; kwargs...)
end

function check_solution(code)
    code = Symbol(code)
    good = (:Success, :Terminated)
    code in good || @warn "Solver failed with status $code (expected $(join(String.(good), " or ")))"
end

# TODO: don't select time points as 2nd/3rd index, since these points will vary
const SymbolicIndex = Union{Num, AbstractArray{Num}}
function Base.getindex(sol::CosmologySolution, i::SymbolicIndex)
    if ModelingToolkit.isparameter(i) && i !== t # don't catch independent variable as parameter
        return sol.th.ps[i] # assume all parameters are in background/thermodynamics # TODO: index sol directly when this is fixed? https://github.com/SciML/ModelingToolkit.jl/issues/3267
    else
        return sol.th[i]
    end
end
function Base.getindex(sol::CosmologySolution, i::SymbolicIndex, j)
    return stack(sol.th[i, j])
end
Base.getindex(sol::CosmologySolution, i::Int, j::SymbolicIndex, k = :) = sol.pts[i][j, k]
Base.getindex(sol::CosmologySolution, i, j::SymbolicIndex, k = :) = [stack(sol[_i, j, k]) for _i in i]
Base.getindex(sol::CosmologySolution, i::Colon, j::SymbolicIndex, k = :) = sol[1:length(sol.pts), j, k]

function (sol::CosmologySolution)(ts::AbstractArray, is::AbstractArray)
    tmin, tmax = extrema(sol.th.t[[begin, end]])
    minimum(ts) >= tmin || minimum(ts) ≈ tmin || throw("Requested time t = $(minimum(ts)) is before initial time $tmin")
    maximum(ts) <= tmax || maximum(ts) ≈ tmax || throw("Requested time t = $(maximum(ts)) is after final time $tmax")
    return permutedims(sol.th(ts, idxs=is)[:, :])
end

function neighboring_modes_indices(sol::CosmologySolution, k)
    k = k_dimensionless.(k, sol.bg.ps[:h])
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

function (sol::CosmologySolution)(ks::AbstractArray, ts::AbstractArray, is::AbstractArray)
    ks = k_dimensionless.(ks, sol.bg.ps[:h])
    isempty(sol.ks) && throw(error("No perturbations solved for. Pass ks to solve()."))
    kmin, kmax = extrema(sol.ks)
    minimum(ks) >= kmin || throw("Requested wavenumber k = $(minimum(ks)) is below the minimum solved wavenumber $kmin")
    maximum(ks) <= kmax || throw("Requested wavenumber k = $(maximum(ks)) is above the maximum solved wavenumber $kmax")
    # disabled; sensitive to exact endpoint values
    #tmin, tmax = extrema(sol.th.t)
    #minimum(ts) >= tmin || throw("Requested time t = $(minimum(ts)) is below minimum solved time $tmin")
    #maximum(ts) <= tmax || throw("Requested time t = $(maximum(ts)) is above maximum solved time $tmin")

    # Pre-allocate intermediate and output arrays
    T = eltype(sol.pts[1])
    v1 = zeros(T, (length(is), length(ts)))
    v2 = zeros(T, (length(is), length(ts)))
    out = zeros(T, length(ks), length(ts), length(is))

    i1_prev, i2_prev = 0, 0 # cache previous looked up solution and reuse it, if possible
    for ik in eachindex(ks) # TODO: multithreading leads to trouble; what about tmap?
        k = ks[ik]
        # Find two wavenumbers to interpolate between
        i1, i2 = neighboring_modes_indices(sol, k)
        k1 = sol.ks[i1]
        k2 = sol.ks[i2]

        # Evaluate solutions for neighboring wavenumbers,
        # but reuse those from the previous iteration if we are still between the same neighboring wavenumbers
        if i1 != i1_prev
            v1 .= sol.pts[i1](ts; idxs=is)[:, :] # TODO: make in-place (https://github.com/SciML/OrdinaryDiffEq.jl/issues/2562)
            i1_prev = i1
        end
        if i2 != i2_prev
            v2 .= sol.pts[i2](ts; idxs=is)[:, :] # TODO: getu or similar for speed? possible while preserving interpolation?
            i2_prev = i2
        end
        v = v1
        if i1 != i2
            # interpolate between solutions
            w = (log(k) - log(k1)) / (log(k2) - log(k1)) # quantities vary smoothest (?) when interpolated in log(k) # TODO: cubic spline?
            v += (v2 - v1) * w # add to v1 from above
        end
        for ii in eachindex(is)
            out[ik, :, ii] = v[ii, :]
        end
    end

    return out
end

# Handle (ts, is) or (ks, ts, is) of arbitrary 0-dimensional and 1-dimensional combinations
function (sol::CosmologySolution)(args...)
    # Please read this function with a pirate's voice
    args_arr = [arg isa Number ? [arg] : arg for arg in args]
    args_outi = [arg isa Number ? 1 : Colon() for arg in args]
    out = sol(args_arr...) # convert to all-array call
    return out[args_outi...] # pick out dimensions for scalar ks/ts/is
end

function (sol::CosmologySolution)(tvar::Num, t, idxs)
    ts = timeseries.(sol, tvar, t)
    return sol(ts, idxs)
end

function (sol::CosmologySolution)(tvar::Num, k, t, idxs)
    ts = timeseries.(sol, tvar, t)
    return sol(k, ts, idxs)
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
        ts = exp.(extend_array(log.(ts), Nextra))
    end
    return ts
end
"""
    timeseries(sol::CosmologySolution, var, val::Number)

Find the time(s) when some variable equals some value with root finding.
"""
function timeseries(sol::CosmologySolution, var, val::Number)
    allequal(sign.(diff(sol[var]))) || error("$var is not monotonic")
    f(t) = sol(t, var) - val # var(t) == val when f(t) == 0
    tspan = extrema(sol[t])
    return find_zero(f, tspan)
end
"""
    timeseries(sol::CosmologySolution, var, vals::AbstractArray)

Find the times when some variable equals some values with splining.
"""
function timeseries(sol::CosmologySolution, var, vals::AbstractArray; kwargs...)
    ts = timeseries(sol; kwargs...)
    xs = sol(ts, var)
    spl = spline(ts, xs)
    ts = spl(vals)
    return ts
end

"""
    shoot(M::CosmologyModel, pars_fixed, pars_varying, conditions; solver = TrustRegion(), verbose = false, kwargs...)

Solve a cosmological model with fixed parameters, while varying some parameters (from their initial guesses) until the given `conditions` at the end time are satisfied.
"""
function shoot(M::CosmologyModel, pars_fixed, pars_varying, conditions; solver = TrustRegion(), verbose = false, kwargs...)
    funcs = [eq.lhs - eq.rhs for eq in conditions] .|> ModelingToolkit.wrap # expressions that should be 0 # TODO: shouldn't have to wrap

    function f(vals_varying, _)
        pars = merge(pars_fixed, Dict(keys(pars_varying) .=> vals_varying))
        prob = CosmologyProblem(M, pars) # TODO: remake
        sol = solve(prob; kwargs...)
        return sol[funcs, :][:, end] # evaluate all expressions at final time (e.g. today)
    end

    guess = collect(values(pars_varying))
    prob = NonlinearProblem(f, guess)
    sol = solve(prob, solver; show_trace = Val(verbose)) # TODO: speed up!
    check_solution(sol.retcode)
    return Dict(keys(pars_varying) .=> sol.u)
end

# Fix model/solution under broadcasted calls
Base.broadcastable(M::CosmologyModel) = Ref(M)
Base.broadcastable(sol::CosmologySolution) = Ref(sol)
