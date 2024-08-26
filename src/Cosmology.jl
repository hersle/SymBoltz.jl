import CommonSolve: solve
#import SymbolicIndexingInterface: all_variable_symbols, getname

struct CosmologyModel
    sys::ODESystem

    bg::ODESystem
    th::ODESystem
    pt::ODESystem
end

function CosmologyModel(sys::ODESystem)
    bg = structural_simplify(background(sys))
    th = structural_simplify(thermodynamics(sys))
    pt = structural_simplify(perturbations(sys))
    return CosmologyModel(sys, bg, th, pt)
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
parameters(M::CosmologyModel) = parameters(M.sys)
initialization_equations(M::CosmologyModel) = initialization_equations(M.sys)
defaults(M::CosmologyModel) = defaults(M.sys)

Base.show(io::IO, M::CosmologyModel) = print(io, chop(sprint(print_tree, M.sys))) # chop off last excessive newline

struct CosmologySolution
    bg::ODESolution
    ks::Array{Float64}
    pts::Union{EnsembleSolution, Nothing}
end

function Base.show(io::IO, sol::CosmologySolution)
    print(io, "Cosmology solution with stages")
    print(io, "\n  1. background: solved with $(nameof(typeof(sol.bg.alg))), $(length(sol.bg)) points")
    if !isnothing(sol.pts)
        solver = nameof(typeof(only(unique(map(pt -> pt.alg, sol.pts)))))
        kmin, kmax = extrema(map(pt -> pt.prob.ps[SymBoltz.k], sol.pts))
        nmin, nmax = extrema(map(pt -> length(pt), sol.pts))
        n = length(sol.pts)
        print(io, "\n  2. perturbations: solved with $solver, $nmin-$nmax points, x$(n) k ∈ [$kmin, $kmax] H₀/c (linear interpolation in-between)")
    end
end

# TODO: don't select time points as 2nd/3rd index, since these points will vary
const SymbolicIndex = Union{Num, AbstractArray{Num}}
Base.getindex(sol::CosmologySolution, i::SymbolicIndex, j = :) = stack(sol.bg[i, j])
Base.getindex(sol::CosmologySolution, i::Int, j::SymbolicIndex, k = :) = sol.pts[i][j, k]
Base.getindex(sol::CosmologySolution, i, j::SymbolicIndex, k = :) = [stack(sol[_i, j, k]) for _i in i]
Base.getindex(sol::CosmologySolution, i::Colon, j::SymbolicIndex, k = :) = sol[1:length(sol.pts), j, k]

function (sol::CosmologySolution)(t, idxs)
    v = sol.bg(t, idxs=idxs)
    if t isa AbstractArray
        v = v.u
        if idxs isa AbstractArray
            v = permutedims(stack(v))
        end
    end
    return v
end

function (sol::CosmologySolution)(k::Number, t, idxs)
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
    tmin, tmax = extrema(sol[SymBoltz.t])
    ts = exp.(range(log(tmin), log(tmax), length = 1000))
    xs = sol(ts, tvar)
    ts = CubicSpline(ts, xs; extrapolate=true)(t)
    return sol(ts, idxs)
end

# TODO: change argument order to join with previous
function (sol::CosmologySolution)(tvar::Num, k, t, idxs)
    tmin, tmax = extrema(sol[SymBoltz.t])
    ts = exp.(range(log(tmin), log(tmax), length = 1000))
    xs = sol(ts, tvar)
    ts = CubicSpline(ts, xs; extrapolate=true)(t)
    return sol(k, ts, idxs)
end

# TODO: add generic function spline(sys::ODESystem, how_to_spline_different_vars) that splines the unknowns of a simplified ODESystem 
# TODO: use CommonSolve.step! to iterate background -> thermodynamics -> perturbations?
"""
    solve(prob::CosmologyModel, pars; tini = 1e-5, aend = 1e0, solver = Rodas5P(), reltol = 1e-13, kwargs...)

Solve `CosmologyModel` with parameters `pars` at the background level.
"""
function solve(prob::CosmologyModel, pars; tini = 1e-5, aend = 1e0, solver = Rodas5P(), reltol = 1e-13, kwargs...)
    ode_prob = ODEProblem(prob.bg, [], (tini, 4.0), pars; use_union = false)
    callback = callback_terminator(prob.bg, prob.bg.g.a, aend)
    ode_sol = solve(ode_prob, solver; callback, reltol, kwargs...)
    return CosmologySolution(ode_sol, [], nothing)
end

function solve(prob::CosmologyModel, pars, ks::AbstractArray; tini = 1e-5, aend = 1e0, solver = KenCarp47(), reltol = 1e-9, verbose = false, kwargs...)
    !issorted(ks) && throw(error("ks = $ks are not sorted in ascending order"))

    tend = 4.0
    if :b₊rec₊dτspline in Symbol.(parameters(prob.pt))
        bg_sol = solve(prob, pars; tini, aend) # TODO: forward kwargs...?
        tend = bg_sol.bg[t][end]
        ts = exp.(range(log(tini), log(tend), length=1024)) # TODO: select determine points adaptively from th_sol # TODO: CMB spectrum is sensitive to number of points here!
        pars = [pars;
            prob.pt.b.rec.dτspline => spline(bg_sol.bg(ts, idxs=log(-prob.bg.b.rec.dτ)).u, log.(ts)) # TODO: improve spline accuracy
            #prob.pt.th.rec.cs²spline => spline(log.(th_sol(ts, idxs=prob.th.rec.Tb).u), log.(ts)),
        ]
    end

    ki = 1.0
    ode_prob0 = ODEProblem(prob.pt, [], (tini, tend), [pars; k => ki]) # TODO: why do I need this???
    #sol0 = solve(prob0, solver; reltol, kwargs...)
    ode_probs = EnsembleProblem(; safetycopy = false, prob = ode_prob0, prob_func = (ode_prob, i, _) -> begin
        verbose && println("$i/$(length(ks)) k = $(ks[i]*k0) Mpc/h")
        return ODEProblem(prob.pt, [], (tini, tend), [pars; k => ks[i]]; use_union = false) # TODO: use remake https://github.com/SciML/OrdinaryDiffEq.jl/pull/2228, https://github.com/SciML/ModelingToolkit.jl/issues/2799 etc. is fixed
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
    ode_sols = solve(ode_probs, solver, EnsembleThreads(), trajectories = length(ks); reltol, progress=true, kwargs...) # TODO: test GPU parallellization
    return CosmologySolution(bg_sol.bg, ks, ode_sols)
end

function solve(M::CosmologyModel, pars, k::Number; kwargs...)
    return solve(M, pars, [k]; kwargs...)
end