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

struct CosmologyProblem
    M::ODESystem

    bg::Union{ODEProblem, Nothing}
    th::Union{ODEProblem, Nothing}
    pt::Union{ODEProblem, Nothing}

    vars_fixed::Dict
    pars_fixed::Dict
    vars_shoot::Dict
    pars_shoot::Dict
    conditions::Vector{Num}

    spline_thermo::Bool
end

struct CosmologySolution
    prob::CosmologyProblem
    bg::Union{ODESolution, Nothing}
    th::Union{ODESolution, Nothing}
    ks::AbstractArray
    pts::Union{EnsembleSolution, Nothing}
end

algname(alg) = string(nameof(typeof(alg)))
algname(alg::CompositeAlgorithm) = join(algname.(alg.algs), "+")

function Base.show(io::IO, prob::CosmologyProblem; indent = "  ")
    print(io, "Cosmology problem for model ")
    printstyled(io, nameof(prob.bg.f.sys), '\n'; bold = true)

    printstyled(io, "Stages:"; bold = true)
    if !isnothing(prob.bg)
        print(io, '\n', indent, "Background")
    end
    if !isnothing(prob.th)
        print(io, '\n', indent, "Thermodynamics")
    end
    if !isnothing(prob.pt)
        print(io, '\n', indent, "Perturbations")
    end

    if !isempty(prob.vars_fixed) || !isempty(prob.pars_fixed)
        printstyled(io, "\nFixed (independent) parameters:"; bold = true)
        for (var, val) in merge(prob.vars_fixed, prob.pars_fixed)
            print(io, '\n', indent, "$var = $val")
        end
    end
    if !isempty(prob.vars_shoot) || !isempty(prob.pars_shoot)
        printstyled(io, "\nShooting (dependent) parameters:"; bold = true)
        for (var, val) in merge(prob.vars_shoot, prob.pars_shoot)
            print(io, '\n', indent, "$var [guess = $val]")
        end
    end
    if !isempty(prob.conditions)
        printstyled(io, "\nShooting conditions:"; bold=true)
        for condition in prob.conditions
            print(io, '\n', indent, condition, " = 0")
        end
    end
end

function Base.show(io::IO, sol::CosmologySolution; indent = "  ")
    print(io, "Cosmology solution for model ")
    printstyled(io, nameof(sol.prob.bg.f.sys), '\n'; bold = true)

    printstyled(io, "Stages:"; bold = true)
    if !isnothing(sol.bg)
        print(io, '\n', indent, "Background: solved with $(algname(sol.bg.alg)), $(length(sol.bg)) points")
    end
    if !isnothing(sol.th)
        print(io, '\n', indent, "Thermodynamics: solved with $(algname(sol.th.alg)), $(length(sol.th)) points")
    end
    if !isnothing(sol.pts)
        kmin, kmax = extrema(sol.ks)
        nmin, nmax = extrema(map(length, sol.pts))
        n = length(sol.pts)
        print(io, '\n', indent, "Perturbations: solved with $(algname(sol.pts[1].alg)), $nmin-$nmax points, x$n k ∈ [$kmin, $kmax] H₀/c (linear interpolation between log(k))")
    end

    if !isempty(sol.prob.vars_fixed) || !isempty(sol.prob.pars_fixed)
        printstyled(io, "\nFixed (independent) parameters:"; bold = true)
        for (var, val) in merge(sol.prob.vars_fixed, sol.prob.pars_fixed)
            print(io, '\n', indent, "$var = $val")
        end
    end
    if !isempty(sol.prob.vars_shoot) || !isempty(sol.prob.pars_shoot)
        printstyled(io, "\nShooting (dependent) parameters:"; bold = true)
        for var in keys(merge(sol.prob.vars_shoot, sol.prob.pars_shoot))
            var = ModelingToolkit.wrap(var)
            printstyled(io, '\n', indent, "$var = $(sol[var][begin])")
        end
    end
end

# Split parameters into DifferentialEquations' u0 and p convention
function split_u0_p(M::ODESystem, x::Dict)
    pars = intersect(keys(x), parameters(M)) # separate parameters from initial conditions
    vars = setdiff(keys(x), pars) # assume the rest are variables (do it without intersection to capture derivatives initial conditions)
    pars = Dict(par => x[par] for par in pars) # like p
    vars = Dict(var => x[var] for var in vars) # like u0
    return vars, pars
end

function CosmologyProblem(M::ODESystem, pars_fixed; spline_thermo = true, debug = false, aterm = 1.0, thermo = true, debug_initialization = false, guesses = Dict(), jac = false, sparse = false, shoot = Dict(), conditions = [], kwargs...)
    bg = structural_simplify(background(M))
    th = structural_simplify(thermodynamics(M))
    pt = structural_simplify(perturbations(M; spline_thermo))

    if debug
        bg = debug_system(bg)
        th = debug_system(th)
        pt = debug_system(pt)
    end

    vars_fixed, pars_fixed = split_u0_p(M, pars_fixed)
    vars_shoot, pars_shoot = split_u0_p(M, shoot)
    conditions = [eq.lhs - eq.rhs for eq in conditions] .|> ModelingToolkit.wrap # expressions that should be 0 # TODO: shouldn't have to wrap

    # TODO: solve only one of bg,th (based on thermo bool)?
    tspan = (1e-5, 4.0) # TODO: keyword argument
    vars = merge(vars_fixed, vars_shoot)
    pars = merge(pars_fixed, pars_shoot)
    parsk = merge(pars, Dict(M.k => NaN)) # k is unused, but must be set
    bg = ODEProblem(bg, vars, tspan, parsk; guesses, fully_determined = true, jac, sparse)
    th = ODEProblem(th, vars, tspan, parsk; guesses, fully_determined = true, jac, sparse)
    pt = ODEProblem(pt, vars, tspan, parsk; guesses, fully_determined = true, jac, sparse)

    return CosmologyProblem(M, bg, th, pt, vars_fixed, pars_fixed, vars_shoot, pars_shoot, conditions, spline_thermo)
end

# TODO: add generic function spline(sys::ODESystem, how_to_spline_different_vars) that splines the unknowns of a simplified ODESystem 
# TODO: solve thermodynamics only if parameters contain thermodynamics parameters?
# TODO: shoot to reach E = 1 today when integrating forwards
# TODO: want to use ODESolution's solver-specific interpolator instead of error-prone spline
"""
    solve(prob::CosmologyProblem; aterm = 1.0, solver = Rodas4P(), reltol = 1e-10, kwargs...)

Solve the `CosmologyProblem` at the background level.
"""
function solve(prob::CosmologyProblem, ks::AbstractArray = []; aterm = 1.0, bgopts = (alg = Rodas4P(), reltol = 1e-10,), thopts = (alg = Rodas4P(), reltol = 1e-10,), ptopts = (alg = KenCarp4(), reltol = 1e-8,), shootopts = (alg = NewtonRaphson(), reltol = 1e-3), guesses = Dict(), thread = true, jac = false, sparse = false, verbose = false)
    M = prob.bg.f.sys
    callback = callback_terminator(prob.bg, M.g.a, aterm)

    bg = solve(prob.bg; callback, bgopts...)

    # TODO: move shooting problem setup to CosmologyProblem
    vars_shoot = collect(keys(prob.vars_shoot))
    pars_shoot = collect(keys(prob.pars_shoot))
    all_shoot = [vars_shoot; pars_shoot] # defines ordering
    nvars = length(prob.vars_shoot)
    npars = length(prob.pars_shoot)

    if !isempty(all_shoot)
        function f(vals, _)
            # TODO: find more efficient way (https://docs.sciml.ai/ModelingToolkit/dev/examples/remake/)
            u0 = nvars == 0 ? missing : Dict(vars_shoot .=> vals[1:nvars]) # TODO: setsym instead?
            p = npars == 0 ? missing : Dict(ModelingToolkit.wrap.(pars_shoot) .=> vals[nvars+1:end]) # TODO: setsym instead?
            bgprob = remake(prob.bg; u0, p) # TODO: u0 from vars_shoot
            bg = solve(bgprob; callback, bgopts...)
            return bg[prob.conditions][end] # evaluate all expressions at final time (e.g. today)
        end
        shooting_guess = [[prob.vars_shoot[var] for var in vars_shoot]..., [prob.pars_shoot[par] for par in pars_shoot]...] # ... instead of ; for type-stability in case vars/pars is empty
        shooting_prob = NonlinearProblem(f, shooting_guess)
        shooting_sol = solve(shooting_prob; show_trace = Val(verbose), shootopts...) # TODO: speed up!
        check_solution(shooting_sol.retcode)
    else
        shooting_sol = []
    end

    check_solution(bg.retcode)

    # Remove duplicate endpoints due to callback termination
    if bg.t[end] == bg.t[end-1]
        ilast = length(bg.t)
        deleteat!(bg.t, ilast)
        deleteat!(bg.u, ilast)
    end

    M = prob.th.f.sys
    u0 = nvars == 0 ? missing : Dict(vars_shoot .=> shooting_sol[1:nvars]) # TODO: don't repeat, put in function # TODO: remake whole CosmologyProblem, and shoot it instead?
    p = npars == 0 ? missing : Dict(ModelingToolkit.wrap.(pars_shoot) .=> shooting_sol[nvars+1:end]) # TODO: setsym instead?
    thprob = remake(prob.th; tspan = extrema(bg.t), u0, p)
    th = solve(thprob; thopts...)
    check_solution(th.retcode)

    # Offset optical depth, so it's 0 today
    if have(M, :b)
        idx_τ = variable_index(prob.th, M.b.rec.τ)
        for i in 1:length(th.u)
            th.u[i][idx_τ] -= th.u[end][idx_τ]
        end
    end

    if isempty(ks)
        ptsols = nothing
    else
        ks = k_dimensionless.(ks, prob.pars_fixed[prob.bg.f.sys.g.h]) # TODO: always index _fixed?
        !issorted(ks) && throw(error("ks = $ks are not sorted in ascending order"))

        if Threads.nthreads() == 1 && thread
            @warn "Multi-threading was requested, but disabled, since Julia is running with only 1 thread. Restart Julia with more threads (e.g. `julia --threads=auto`) to enable multi-threading, or pass thread = false to explicitly disable it."
            thread = false
        end

        M = prob.pt.f.sys # TODO: can I exploit that the structure of the perturbation ODEs is ẏ = J * y with "constant" J?
        ptprob0 = prob.pt
        u0 = nvars == 0 ? missing : Dict(vars_shoot .=> shooting_sol[1:nvars]) # TODO: don't repeat, put in function # TODO: could do with setsym_oop below, but that seems to keep G.ϕ = 0.95 instead of using the shot value
        p = npars == 0 ? missing : Dict(ModelingToolkit.wrap.(pars_shoot) .=> shooting_sol[nvars+1:end]) # TODO: setsym instead?
        ptprob0 = remake(ptprob0; tspan = extrema(bg.t), u0, p)
    
        update_vars = Dict(M.g.a => bg[M.g.a][begin])
        if have(M, :b) && have(M.b, :rec) && prob.spline_thermo
            update_vars = merge(update_vars, Dict(
                M.b.rec.τspline => spline(th[M.b.rec.τ], th.t), # TODO: when solving thermo with low reltol: even though the solution is correct, just taking its points for splining can be insufficient. should increase number of points, so it won't mess up the perturbations
                M.b.rec.τ̇spline => spline(th[M.b.rec.τ̇], th.t),
                M.b.rec.cₛ²spline => spline(th[M.b.rec.cₛ²], th.t)
            ))
        end
        update = ModelingToolkit.setsym_oop(ptprob0, collect(keys(update_vars)))
        newu0, newp = update(ptprob0, collect(values(update_vars)))
        ptprob0 = remake(ptprob0, tspan = extrema(bg.t), u0 = newu0, p = newp)

        kset! = setp(ptprob0, k) # function that sets k on a problem
        ptprob_tlv = TaskLocalValue{ODEProblem}(() -> deepcopy(ptprob0)) # prevent conflicts where different tasks modify same problem: https://discourse.julialang.org/t/solving-ensembleproblem-efficiently-for-large-systems-memory-issues/116146/11 (alternatively copy just p and u0: https://github.com/SciML/ModelingToolkit.jl/issues/3056)
        ptprobs = EnsembleProblem(; safetycopy = false, prob = ptprob0, prob_func = (_, i, _) -> begin
            ptprob = ptprob_tlv[]
            verbose && println("$i/$(length(ks)) k = $(ks[i]*k0) Mpc/h")
            kset!(ptprob, ks[i])
            return ptprob
        end)
        ensemblealg = thread ? EnsembleThreads() : EnsembleSerial()
        ptsols = solve(ptprobs; ensemblealg, trajectories = length(ks), ptopts...) # TODO: test GPU parallellization
        for i in 1:length(ptsols)
            check_solution(ptsols[i].retcode)
        end
    end

    return CosmologySolution(prob, bg, th, ks, ptsols)
end

function solve(prob::CosmologyProblem, k::Number; kwargs...)
    return solve(prob, [k]; kwargs...)
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
    shoot(M::ODESystem, pars_fixed, pars_varying, conditions; solver = TrustRegion(), verbose = false, kwargs...)

Solve a cosmological model with fixed parameters, while varying some parameters (from their initial guesses) until the given `conditions` at the end time are satisfied.
"""
function shoot(M::ODESystem, pars_fixed, pars_varying, conditions; solver = TrustRegion(), verbose = false, kwargs...)
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
Base.broadcastable(sys::ODESystem) = Ref(sys) # TODO: add to ModelingToolkit?
Base.broadcastable(sol::CosmologySolution) = Ref(sol)
