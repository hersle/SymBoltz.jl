struct CosmologicalModel
    bg::ODESystem
    bg_sim::ODESystem
    pt::ODESystem
    pt_sim::ODESystem
end

function CosmologicalModel(sys::ODESystem)
    bg = background(sys)
    bg_sim = structural_simplify(bg)

    pt = perturbations(sys)
    pt_sim = structural_simplify(pt)

    return CosmologicalModel(bg, bg_sim, pt, pt_sim)
end

@kwdef struct CosmologicalParameters
    Ωγ0 = 5.5e-5
    Ωc0 = 0.267
    Ωb0 = 0.05
    h = 0.67
    Neff = 3.046
    As = 2.1e-9
    Yp = 0.245
end

# TODO: add generic function spline(sys::ODESystem, how_to_spline_different_vars) that splines the unknowns of a simplified ODESystem 
# TODO: just one solve(), but dispatch on whether or not k is specified to determine whether perturbations are solved (and recombination, in a different way)?

function solve_background(M::CosmologicalModel, par::CosmologicalParameters; tini = 1e-5, aend = 1e0, solver = Vern8(), reltol = 1e-15, kwargs...)
    prob = ODEProblem(M.bg_sim, [], (tini, 4.0), [
        M.bg_sim.r.Ω0 => par.Ωγ0,
        M.bg_sim.m.Ω0 => par.Ωc0,
        #M.bg_sim.bar.Ω0 => par.Ωb0,
        #M.bg_sim.neu.Neff => par.Neff,
        #M.th.bg.g.h => par.h
    ])
    callback = callback_terminator(M.bg_sim, M.bg_sim.g.a, aend)
    return solve(prob, solver; callback, reltol, kwargs...)
end

function solve_thermodynamics(M::CosmologicalModel, par::CosmologicalParameters; tini = 1e-5, aend = 1e0, solver = Rodas5P(), reltol = 1e-13, kwargs...) # need very small tolerance to get good csb²
    return solve_background(M, par; tini, aend, solver, reltol, kwargs...) # TODO: reenable true thermodynamics solution
    prob = ODEProblem(model.th_sim, [], (tini, 4.0), [
        model.th_sim.bg.ph.Ω0 => par.Ωγ0,
        model.th_sim.bg.cdm.Ω0 => par.Ωc0,
        model.th_sim.bg.bar.Ω0 => par.Ωb0,
        model.th_sim.bg.g.h => par.h,
        model.th_sim.bg.neu.Neff => par.Neff,
        model.th_sim.rec.Yp => par.Yp
    ])
    callback = callback_terminator(model.bg_sim, model.bg_sim.g.a, aend)
    return solve(prob, solver; callback, reltol, kwargs...)
end

function solve_perturbations(M::CosmologicalModel, ks::AbstractArray, par::CosmologicalParameters; tini = 1e-5, aend = 1e0, solver = KenCarp47(), reltol = 1e-9, verbose = false, kwargs...)
    pars = Pair{Any, Any}[ # TODO: avoid Any
        M.pt_sim.r.Ω0 => par.Ωγ0,
        M.pt_sim.m.Ω0 => par.Ωc0,
        #M.th.bg.bar.Ω0 => par.Ωb0,
        #M.th.bg.g.h => par.h,
        #M.th.bg.neu.Neff => par.Neff,
    ]

    if true # TODO: reenable
        th_sol = solve_thermodynamics(M, par; tini, aend) # TODO: forward kwargs...?
        tend = th_sol[t][end]
        #=
        ts = exp.(range(log(tini), log(tend), length=1024)) # TODO: select determine points adaptively from th_sol # TODO: CMB spectrum is sensitive to number of points here!
        push!(pars,
            M.pt_sim.th.rec.dτspline => spline(log.(.-th_sol(ts, idxs=M.th.rec.dτ).u), log.(ts)), # TODO: improve spline accuracy
            M.pt_sim.th.rec.cs²spline => spline(log.(th_sol(ts, idxs=M.th.rec.Tb).u), log.(ts)),
        )
        =#
    end

    ki = 1.0
    prob0 = ODEProblem(M.pt_sim, [], (tini, tend), [pars; k => ki]) # TODO: why do I need this???
    #sol0 = solve(prob0, solver; reltol, kwargs...)
    probs = EnsembleProblem(; safetycopy = false, prob = prob0, prob_func = (prob, i, _) -> begin
        verbose && println("$i/$(length(ks)) k = $(ks[i]*k0) Mpc/h")
        return ODEProblem(M.pt_sim, [], (tini, tend), [pars; k => ks[i]]) # TODO: use remake https://github.com/SciML/OrdinaryDiffEq.jl/pull/2228, https://github.com/SciML/ModelingToolkit.jl/issues/2799 etc. is fixed
        #= # TODO: this should work if I use defaults for perturbation ICs, but that doesnt work as it should because the initialization system becomes overdefined and 
        prob_new = remake(prob, u0 = [
            model.pt_sim.th.bg.g.a => sol0[model.pt_sim.th.bg.g.a][begin]
            model.pt_sim.g1.Φ => sol0[model.pt_sim.g1.Φ][begin]
            model.pt_sim.cdm.θ => (ks[i]/ki)^2 * sol0[model.pt_sim.cdm.θ][begin]
            model.pt_sim.bar.θ => (ks[i]/ki)^2 * sol0[model.pt_sim.bar.θ][begin]
            model.pt_sim.ph.F[1] => (ks[i]/ki)^1 * sol0[model.pt_sim.ph.F[1]][begin]
            model.pt_sim.neu.F[1] => (ks[i]/ki)^1 * sol0[model.pt_sim.neu.F[1]][begin]
            model.pt_sim.neu.F[2] => (ks[i]/ki)^2 * sol0[model.pt_sim.neu.F[2]][begin]
            collect(model.pt_sim.mneu.ψ[:,1] .=> (ks[i]/ki)^1 * sol0[model.pt_sim.mneu.ψ[:,1]][begin])...
            collect(model.pt_sim.mneu.ψ[:,2] .=> (ks[i]/ki)^2 * sol0[model.pt_sim.mneu.ψ[:,2]][begin])...
        ], tspan = (tini, 4.0), p = [k => ks[i]], use_defaults = true)
        return prob_new # BUG: prob_new's u0 does not match solution[begin]
        =#
    end)
    return solve(probs, solver, EnsembleThreads(), trajectories = length(ks); reltol, progress=true, kwargs...) # TODO: test GPU parallellization
end

function solve_perturbations(M::CosmologicalModel, ks::Number, par::CosmologicalParameters; kwargs...)
    return solve_perturbations(M, [ks], par; kwargs...)[1]
end