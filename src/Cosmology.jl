import CommonSolve: solve

# TODO: split further into BackgroundProblem and PerturbationProblem
struct CosmologyProblem
    bg_sim::ODESystem # TODO: remove _sim suffix
    pt_sim::ODESystem
end

function CosmologyProblem(model::ODESystem)
    bg = background(model)
    bg_sim = structural_simplify(bg)

    pt = perturbations(model)
    pt_sim = structural_simplify(pt)

    return CosmologyProblem(bg_sim, pt_sim)
end

# TODO: add generic function spline(sys::ODESystem, how_to_spline_different_vars) that splines the unknowns of a simplified ODESystem 
# TODO: use CommonSolve.step! to iterate background -> thermodynamics -> perturbations?
function solve(prob::CosmologyProblem, pars; tini = 1e-5, aend = 1e0, solver = Rodas5P(), reltol = 1e-13, kwargs...)
    ode_prob = ODEProblem(prob.bg_sim, [], (tini, 4.0), pars)
    callback = callback_terminator(prob.bg_sim, prob.bg_sim.g.a, aend)
    return solve(ode_prob, solver; callback, reltol, kwargs...)
end

function solve(prob::CosmologyProblem, pars, ks::AbstractArray; tini = 1e-5, aend = 1e0, solver = KenCarp47(), reltol = 1e-9, verbose = false, kwargs...)
    tend = 4.0
    if :b₊rec₊dτspline in Symbol.(parameters(prob.pt_sim))
        bg_sol = solve(prob, pars; tini, aend) # TODO: forward kwargs...?
        tend = bg_sol[t][end]
        ts = exp.(range(log(tini), log(tend), length=1024)) # TODO: select determine points adaptively from th_sol # TODO: CMB spectrum is sensitive to number of points here!
        pars = [pars;
            prob.pt_sim.b.rec.dτspline => spline(log.(.-bg_sol(ts, idxs=prob.bg_sim.b.rec.dτ).u), log.(ts)) # TODO: improve spline accuracy
            #prob.pt_sim.th.rec.cs²spline => spline(log.(th_sol(ts, idxs=prob.th.rec.Tb).u), log.(ts)),
        ]
    end

    ki = 1.0
    ode_prob0 = ODEProblem(prob.pt_sim, [], (tini, tend), [pars; k => ki]) # TODO: why do I need this???
    #sol0 = solve(prob0, solver; reltol, kwargs...)
    ode_probs = EnsembleProblem(; safetycopy = false, prob = ode_prob0, prob_func = (ode_prob, i, _) -> begin
        verbose && println("$i/$(length(ks)) k = $(ks[i]*k0) Mpc/h")
        return ODEProblem(prob.pt_sim, [], (tini, tend), [pars; k => ks[i]]) # TODO: use remake https://github.com/SciML/OrdinaryDiffEq.jl/pull/2228, https://github.com/SciML/ModelingToolkit.jl/issues/2799 etc. is fixed
        #= # TODO: this should work if I use defaults for perturbation ICs, but that doesnt work as it should because the initialization system becomes overdefined and 
        prob_new = remake(prob, u0 = [
            M.pt_sim.th.bg.g.a => sol0[M.pt_sim.th.bg.g.a][begin]
            M.pt_sim.g1.Φ => sol0[M.pt_sim.g1.Φ][begin]
            M.pt_sim.cdm.θ => (ks[i]/ki)^2 * sol0[M.pt_sim.cdm.θ][begin]
            M.pt_sim.bar.θ => (ks[i]/ki)^2 * sol0[M.pt_sim.bar.θ][begin]
            M.pt_sim.ph.F[1] => (ks[i]/ki)^1 * sol0[M.pt_sim.ph.F[1]][begin]
            M.pt_sim.neu.F[1] => (ks[i]/ki)^1 * sol0[M.pt_sim.neu.F[1]][begin]
            M.pt_sim.neu.F[2] => (ks[i]/ki)^2 * sol0[M.pt_sim.neu.F[2]][begin]
            collect(M.pt_sim.mneu.ψ[:,1] .=> (ks[i]/ki)^1 * sol0[M.pt_sim.mneu.ψ[:,1]][begin])...
            collect(M.pt_sim.mneu.ψ[:,2] .=> (ks[i]/ki)^2 * sol0[M.pt_sim.mneu.ψ[:,2]][begin])...
        ], tspan = (tini, 4.0), p = [k => ks[i]], use_defaults = true)
        return prob_new # BUG: prob_new's u0 does not match solution[begin]
        =#
    end)
    return solve(ode_probs, solver, EnsembleThreads(), trajectories = length(ks); reltol, progress=true, kwargs...) # TODO: test GPU parallellization
end

function solve(M::CosmologyProblem, pars, k::Number; kwargs...)
    return solve(M, pars, [k]; kwargs...)[1]
end