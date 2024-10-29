import CommonSolve: solve
#import SymbolicIndexingInterface: all_variable_symbols, getname

function background(sys)
    sys = thermodynamics(sys)
    if :b in ModelingToolkit.get_name.(ModelingToolkit.get_systems(sys))
        sys = replace(sys, sys.b.rec => ODESystem([], t; name = :rec)) # remove recombination
    end
    return sys
end

function thermodynamics(sys)
    return transform((sys, _) -> extract_order(sys, [0]), sys)
end

function perturbations(sys; spline_thermo = true)
    if :b in ModelingToolkit.get_name.(ModelingToolkit.get_systems(sys)) && spline_thermo && !all([eq.rhs === 0 for eq in equations(sys.b.rec)])
        @named rec = thermodynamics_recombination_splined()
        sys = replace(sys, sys.b.rec => rec) # substitute in splined recombination
    end
    return transform((sys, _) -> extract_order(sys, [0, 1]), sys)
end

struct CosmologyModel
    sys::ODESystem

    bg::ODESystem
    th::ODESystem
    pt::ODESystem

    initE::Bool # need extra equation E = 1 to initialize when integrating background backwards?
    spline_thermo::Bool # whether to spline thermodynamics in perturbations
end

function CosmologyModel(sys::ODESystem; initE = true, spline_thermo = true, debug = false)
    if debug
        sys = debugize(sys) # TODO: make work with massive neutrinos
    end
    bg = structural_simplify(background(sys))
    th = structural_simplify(thermodynamics(sys))
    pt = structural_simplify(perturbations(sys; spline_thermo))
    sys = complete(sys; flatten = false)
    return CosmologyModel(sys, bg, th, pt, initE, spline_thermo)
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
equations(M::CosmologyModel) = equations(M.sys)
observed(M::CosmologyModel) = observed(M.sys)
unknowns(M::CosmologyModel) = unknowns(M.sys)
parameters(M::CosmologyModel) = parameters(M.sys)
initialization_equations(M::CosmologyModel) = initialization_equations(M.sys)
defaults(M::CosmologyModel) = defaults(M.sys)
hierarchy(M::CosmologyModel) = hierarchy(M.sys)
Base.show(io::IO, mime::MIME"text/plain", M::CosmologyModel) = show(io, mime, M.sys) # chop off last excessive newline

struct CosmologySolution
    bg::ODESolution
    th::ODESolution
    ks::AbstractArray
    pts::Union{EnsembleSolution, Nothing}
end

function Base.show(io::IO, sol::CosmologySolution)
    print(io, "Cosmology solution with stages")
    print(io, "\n  1. background: solved with $(nameof(typeof(sol.bg.alg))), $(length(sol.bg)) points")
    if !isnothing(sol.th)
        print(io, "\n  2. thermodynamics: solved with $(nameof(typeof(sol.th.alg))), $(length(sol.th)) points")
    end
    if !isnothing(sol.pts)
        solver = nameof(typeof(only(unique(map(pt -> pt.alg, sol.pts)))))
        kmin, kmax = extrema(map(pt -> pt.prob.ps[SymBoltz.k], sol.pts))
        nmin, nmax = extrema(map(pt -> length(pt), sol.pts))
        n = length(sol.pts)
        print(io, "\n  3. perturbations: solved with $solver, $nmin-$nmax points, x$(n) k ∈ [$kmin, $kmax] H₀/c (linear interpolation in-between)")
    end
end

# TODO: add generic function spline(sys::ODESystem, how_to_spline_different_vars) that splines the unknowns of a simplified ODESystem 
# TODO: use CommonSolve.step! to iterate background -> thermodynamics -> perturbations?
"""
    solve(M::CosmologyModel, pars; aini = 1e-7, solver = Rodas5P(), reltol = 1e-10, kwargs...)

Solve `CosmologyModel` with parameters `pars` at the background level.
"""
# TODO: solve thermodynamics only if parameters contain thermodynamics parameters?
# TODO: shoot to reach E = 1 today when integrating forwards
function solve(M::CosmologyModel, pars; aini = 1e-7, solver = Rodas5P(), reltol = 1e-15, backwards = true, thermo = true, debug_initialization = false, guesses = Dict(), kwargs...)
    # Split parameters into DifferentialEquations' "u0" and "p" convention # TODO: same in perturbations
    params = merge(pars, Dict(M.k => 0.0)) # k is unused, but must be set https://github.com/SciML/ModelingToolkit.jl/issues/3013 # TODO: remove
    pars = intersect(keys(params), parameters(M)) # separate parameters from initial conditions
    vars = setdiff(keys(params), pars) # assume the rest are variables (do it without intersection to capture derivatives initial conditions)
    vars = Dict(var => params[var] for var in vars) # like u0
    pars = Dict(par => params[par] for par in pars) # like p

    # First solve background forwards or backwards from today
    if backwards
        push!(vars, M.g.a => 1.0)
        M.initE && push!(vars, M.g.ℰ => 1.0)
        tspan = (0.0, -4.0) # integrate backwards
        aterm = aini # terminate at initial scale factor
    else
        push!(vars, M.g.a => aini)
        tspan = (0.0, +4.0) # integrate forwards
        aterm = 1.0 # terminate today
    end
    bg_prob = ODEProblem(M.bg, vars, tspan, pars; guesses, fully_determined = true)
    callback = callback_terminator(M.bg, M.g.a, aterm)
    debug_initialization = debug_initialization && !isnothing(bg_prob.f.initializeprob)
    if debug_initialization
        isys = bg_prob.f.initializeprob.f.sys
        println("Solving initialization equations")
        println(join(equations(isys), "\n"))
        println("for unknowns ", join(unknowns(isys), ", "), ":")
        solve(bg_prob.f.initializeprob; show_trace = Val(true))
    end
    bg_sol = solve(bg_prob, solver; callback, reltol, kwargs...)
    if debug_initialization
        vars = unknowns(isys) .=> bg_sol[unknowns(isys)][begin]
        println("Found variable values")
        println(join(vars, "\n"))
    end
    check_solution(bg_sol.retcode)

    # Then solve thermodynamics forwards till today
    if thermo
        iini = backwards ? length(bg_sol) : 1 # index corresponding to earliest time
        tini = 1 / bg_sol[M.g.ℰ][iini]
        Δt = abs(bg_sol[M.t][end] - tspan[1])
        tend = tini + Δt
        ics = unknowns(M.bg) .=> bg_sol[unknowns(M.bg)][iini]
        ics = filter(ic -> !contains(String(Symbol(ic.first)), "aˍt(t)"), ics) # remove ȧ initial condition
        th_prob = ODEProblem(M.th, ics, (tini, tend), pars; guesses, fully_determined = true)
        th_sol = solve(th_prob, solver; reltol, kwargs...)
        check_solution(th_sol.retcode)
    else
        th_sol = bg_sol
    end

    return CosmologySolution(bg_sol, th_sol, [], nothing)
end

function solve(M::CosmologyModel, pars, ks::AbstractArray; aini = 1e-7, solver = KenCarp4(), reltol = 1e-11, backwards = true, verbose = true, kwargs...)
    ks = k_dimensionless.(ks, pars[M.g.h])

    !issorted(ks) && throw(error("ks = $ks are not sorted in ascending order"))

    th_sol = solve(M, pars; aini, backwards, kwargs...)
    tini, tend = extrema(th_sol.th[t])
    if M.spline_thermo
        th_sol_spline = isempty(kwargs) ? th_sol : solve(M, pars; aini, backwards) # should solve again if given keyword arguments, like saveat
        pars = merge(pars, Dict(
            M.pt.b.rec.τspline => spline(th_sol_spline[M.b.rec.τ] .- th_sol_spline[M.b.rec.τ][end], th_sol_spline[M.t]), # TODO: more time points, spline log(t)?
            M.pt.b.rec.cs²spline => spline(th_sol_spline[M.b.rec.cs²], th_sol_spline[M.t])
        ))
    end

    ki = 1.0
    ics0 = unknowns(M.bg) .=> th_sol.bg[unknowns(M.bg)][backwards ? end : begin]
    ics0 = filter(ic -> !contains(String(Symbol(ic.first)), "aˍt"), ics0) # remove D(a)
    ics0 = Dict(ics0)
    push!(pars, k => ki)
    ode_prob0 = ODEProblem(M.pt, ics0, (tini, tend), pars; fully_determined = true) # TODO: why do I need this???
    #sol0 = solve(prob0, solver; reltol, kwargs...)
    ode_probs = EnsembleProblem(; safetycopy = false, prob = ode_prob0, prob_func = (ode_prob, i, _) -> begin
        verbose && println("$i/$(length(ks)) k = $(ks[i]*k0) Mpc/h")
        return ODEProblem(M.pt, ics0, (tini, tend), merge(pars, Dict(k => ks[i])), fully_determined = true) # TODO: use remake https://github.com/SciML/OrdinaryDiffEq.jl/pull/2228, https://github.com/SciML/ModelingToolkit.jl/issues/2799 etc. is fixed
        #= # TODO: this should work if I use defaults for perturbation ICs, but that doesnt work as it should because the initialization system becomes overdefined and 
        prob_new = remake(prob, u0 = [
            M.pt.th.bg.g.a => sol0[M.pt.th.bg.g.a][begin]
            M.pt.g1.Φ => sol0[M.pt.g1.Φ][begin]
            M.pt.cdm.θ => (ks[i]/ki)^2 * sol0[M.pt.cdm.θ][begin]
            M.pt.bar.θ => (ks[i]/ki)^2 * sol0[M.pt.bar.θ][begin]
            M.pt.ph.F[1] => (ks[i]/ki)^1 * sol0[M.pt.ph.F[1]][begin]
            M.pt.neu.F[1] => (ks[i]/ki)^1 * sol0[M.pt.neu.F[1]][begin]
            M.pt.neu.F[2] => (ks[i]/ki)^2 * sol0[M.pt.neu.F[2]][begin]
            collect(M.pt.mneu.ψ[:,1] .=> (ks[i]/ki)^1 * sol0[M.pt.mneu.ψ[:,1]][begin])...
            collect(M.pt.mneu.ψ[:,2] .=> (ks[i]/ki)^2 * sol0[M.pt.mneu.ψ[:,2]][begin])...
        ], tspan = (tini, 4.0), p = [k => ks[i]], use_defaults = true)
        return prob_new # BUG: prob_new's u0 does not match solution[begin]
        =#
    end)
    ode_sols = solve(ode_probs, solver, EnsembleThreads(), trajectories = length(ks); reltol, kwargs...) # TODO: test GPU parallellization
    for i in 1:length(ode_sols)
        check_solution(ode_sols[i].retcode)
    end
    return CosmologySolution(th_sol.bg, th_sol.th, ks, ode_sols)
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
function Base.getindex(sol::CosmologySolution, i::SymbolicIndex, j = :)
    if ModelingToolkit.isparameter(i) && i !== t && (j == :) # don't catch independent variable as parameter
        return sol.th.ps[i] # assume all parameters are in background/thermodynamics
    else
        return stack(sol.th[i, j])
    end
end
Base.getindex(sol::CosmologySolution, i::Int, j::SymbolicIndex, k = :) = sol.pts[i][j, k]
Base.getindex(sol::CosmologySolution, i, j::SymbolicIndex, k = :) = [stack(sol[_i, j, k]) for _i in i]
Base.getindex(sol::CosmologySolution, i::Colon, j::SymbolicIndex, k = :) = sol[1:length(sol.pts), j, k]

function (sol::CosmologySolution)(t, idxs)
    v = sol.th(t, idxs=idxs)
    if t isa AbstractArray
        v = v.u
        if idxs isa AbstractArray
            v = permutedims(stack(v))
        end
    end
    return v
end

function (sol::CosmologySolution)(k::Number, t, idxs)
    k = k_dimensionless.(k, sol.bg.ps[:h])

    isempty(sol.ks) && throw(error("no perturbations solved for; pass ks to solve()"))

    kmin, kmax = extrema(sol.ks)
    (kmin <= k <= kmax) || throw(error("k = $k is outside range ($kmin, $kmax)"))

    i2 = searchsortedfirst(sol.ks, k) # index above target k
    i1 = i2 - 1 # index below target k ()

    v2 = sol.pts[i2](t; idxs)
    if i1 == 0
        # k == kmin, no interpolation necessary
        v = v2
    else
        # k > kmin, linearly interpolate between k1 and k2
        v1 = sol.pts[i1](t; idxs)
        k1 = sol.ks[i1]
        k2 = sol.ks[i2]
        v = @. v1 + (v2 - v1) * (log(k) - log(k1)) / (log(k2) - log(k1)) # quantities vary smoother when interpolated in log(k) # TODO: use cubic spline?
    end

    if t isa AbstractArray
        v = v.u
        if idxs isa AbstractArray
            v = permutedims(stack(v))
        end
    end
    return v
end

function (sol::CosmologySolution)(k::AbstractArray, t, idxs)
    v = sol.(k, Ref(t), Ref(idxs))
    if t isa AbstractArray && idxs isa AbstractArray
        v = stack(v; dims=1)
    elseif t isa AbstractArray || idxs isa AbstractArray
        v = permutedims(stack(v))
    end
    return v
end

function (sol::CosmologySolution)(tvar::Num, t, idxs)
    ts = sol[SymBoltz.t]
    xs = sol(ts, tvar)
    ts = spline(ts, xs)(t)
    return sol(ts, idxs)
end

# TODO: change argument order to join with previous
function (sol::CosmologySolution)(tvar::Num, k, t, idxs)
    tmin, tmax = extrema(sol[SymBoltz.t])
    ts = exp.(range(log(tmin), log(tmax), length = 1000)) # TODO: select from solution
    xs = sol(ts, tvar)
    ts = CubicSpline(ts, xs; extrapolate=true)(t)
    return sol(k, ts, idxs)
end


function shoot(M::CosmologyModel, pars_fixed, pars_varying, conditions; solver = TrustRegion(), verbose = false, kwargs...)
    funcs = [eq.lhs - eq.rhs for eq in conditions] .|> ModelingToolkit.wrap # expressions that should be 0 # TODO: shouldn't have to wrap

    function f(vals_varying, _)
        pars = merge(pars_fixed, Dict(keys(pars_varying) .=> vals_varying))
        sol = solve(M, pars; kwargs...) # solve cosmology
        return sol[funcs][:, end] # evaluate all expressions at final time (e.g. today)
    end

    guess = collect(values(pars_varying))
    prob = NonlinearProblem(f, guess)
    sol = solve(prob, solver; show_trace = Val(verbose)) # TODO: speed up!
    check_solution(sol.retcode)
    return Dict(keys(pars_varying) .=> sol.u)
end
