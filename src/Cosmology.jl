import CommonSolve: solve
#import SymbolicIndexingInterface: all_variable_symbols, getname

# TODO: split further into BackgroundProblem and PerturbationProblem
struct CosmologyProblem
    bg::ODESystem
    pt::ODESystem
end

function CosmologyProblem(model::ODESystem)
    bg = structural_simplify(background(model))
    pt = structural_simplify(perturbations(model))
    return CosmologyProblem(bg, pt)
end

struct CosmologySolution
    bg::ODESolution
    pt::Union{ODESolution, Nothing} # TODO: multiple ks, source function-like splining?
end

function Base.show(io::IO, sol::CosmologySolution)
    print(io, "Cosmology solution with")
    print(io, "\n* background: solved with $(nameof(typeof(sol.bg.alg))), $(length(sol.bg)) points")
    !isnothing(sol.pt) && print(io, "\n* perturbation (k = $(sol.pt.prob.ps[k]) H₀/c): solved with $(nameof(typeof(sol.pt.alg))), $(length(sol.pt)) points")
end

#=
function subsol_with_variables(sol::CosmologySolution, vars)
    # TODO: what to do with ambiguous t? default to bg, but select pt if requesting only perturbed variables
    vars = getname.(vars)
    if !isnothing(sol.pt)
        bg_vars = getname.(all_variable_symbols(sol.bg)) # variables in background
        pt_vars = getname.(all_variable_symbols(sol.pt)) # variables in perturbations
        pt_vars = setdiff(pt_vars, bg_vars) # variables in perturbations, but not in background
        if !isdisjoint(vars, pt_vars)
            return sol.pt # one or more requested variables are perturbation variables
        end
    end
    return sol.bg
end

function subsol_with_variables(sol::CosmologySolution, var::Number)
    return subsol_with_variables(sol, [var])
end

Base.getindex(sol::CosmologySolution, i) = subsol_with_variables(sol, i)[i]
=#

function Base.getindex(sol::CosmologySolution, i)
    try # TODO: avoid try/catch
        return sol.bg[i]
    catch e
        e isa ArgumentError && return sol.pt[i]
        rethrow(e)
    end
end

function (sol::CosmologySolution)(t, idxs)
    try # TODO: avoid try/catch
        return sol.bg(t, idxs=idxs)
    catch e
        e isa ArgumentError && return sol.pt(t, idxs=idxs)
        rethrow(e)
    end
end

# TODO: add generic function spline(sys::ODESystem, how_to_spline_different_vars) that splines the unknowns of a simplified ODESystem 
# TODO: use CommonSolve.step! to iterate background -> thermodynamics -> perturbations?
function solve(prob::CosmologyProblem, pars; tini = 1e-5, aend = 1e0, solver = Rodas5P(), reltol = 1e-13, kwargs...)
    ode_prob = ODEProblem(prob.bg, [], (tini, 4.0), pars)
    callback = callback_terminator(prob.bg, prob.bg.g.a, aend)
    ode_sol = solve(ode_prob, solver; callback, reltol, kwargs...)
    return CosmologySolution(ode_sol, nothing)
end

function solve(prob::CosmologyProblem, pars, ks::AbstractArray; tini = 1e-5, aend = 1e0, solver = KenCarp47(), reltol = 1e-9, verbose = false, kwargs...)
    tend = 4.0
    if :b₊rec₊dτspline in Symbol.(parameters(prob.pt))
        bg_sol = solve(prob, pars; tini, aend) # TODO: forward kwargs...?
        tend = bg_sol[t][end]
        ts = exp.(range(log(tini), log(tend), length=1024)) # TODO: select determine points adaptively from th_sol # TODO: CMB spectrum is sensitive to number of points here!
        pars = [pars;
            prob.pt.b.rec.dτspline => spline(log.(.-bg_sol(ts, prob.bg.b.rec.dτ).u), log.(ts)) # TODO: improve spline accuracy
            #prob.pt.th.rec.cs²spline => spline(log.(th_sol(ts, idxs=prob.th.rec.Tb).u), log.(ts)),
        ]
    end

    ki = 1.0
    ode_prob0 = ODEProblem(prob.pt, [], (tini, tend), [pars; k => ki]) # TODO: why do I need this???
    #sol0 = solve(prob0, solver; reltol, kwargs...)
    ode_probs = EnsembleProblem(; safetycopy = false, prob = ode_prob0, prob_func = (ode_prob, i, _) -> begin
        verbose && println("$i/$(length(ks)) k = $(ks[i]*k0) Mpc/h")
        return ODEProblem(prob.pt, [], (tini, tend), [pars; k => ks[i]]) # TODO: use remake https://github.com/SciML/OrdinaryDiffEq.jl/pull/2228, https://github.com/SciML/ModelingToolkit.jl/issues/2799 etc. is fixed
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
    end, output_func = (sol, i) -> (CosmologySolution(bg_sol.bg, sol), false))
    return solve(ode_probs, solver, EnsembleThreads(), trajectories = length(ks); reltol, progress=true, kwargs...) # TODO: test GPU parallellization
end

function solve(M::CosmologyProblem, pars, k::Number; kwargs...)
    return solve(M, pars, [k]; kwargs...)[1]
end