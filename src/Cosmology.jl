struct CosmologicalModel
    bg::ODESystem
    th::ODESystem
    pt::ODESystem

    bg_sim::ODESystem
    th_sim::ODESystem
    pt_sim::ODESystem

    spline_th::Bool # TODO: what about a parametric CosmologicalModel?
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

function ΛCDM(; lmax=6, spline_th = true)
    @named bg = background_ΛCDM()
    @named th = thermodynamics_ΛCDM(bg; spline=true) # just for constructing perturbation system
    @named pt = perturbations_ΛCDM(th, lmax; spline_th) # TODO: spline full background OR select quantities
    @named th = thermodynamics_ΛCDM(bg; spline=false)
    bg_sim = structural_simplify(bg)
    th_sim = structural_simplify(th)
    pt_sim = structural_simplify(pt)
    return CosmologicalModel(bg, th, pt, bg_sim, th_sim, pt_sim, spline_th)
end

# TODO: add aini, aend

function solve_background(model::CosmologicalModel, par::CosmologicalParameters; tini = 1e-5, aend = 1e0, solver = Vern8(), reltol = 1e-15, kwargs...)
    prob = ODEProblem(model.bg_sim, [], (tini, 4.0), [
        model.bg_sim.ph.Ω0 => par.Ωγ0,
        model.bg_sim.cdm.Ω0 => par.Ωc0,
        model.bg_sim.bar.Ω0 => par.Ωb0,
        model.bg_sim.neu.Neff => par.Neff,
        model.th.bg.g.h => par.h
    ])
    callback = callback_terminator(model.bg_sim, model.bg_sim.g.a, aend)
    return solve(prob, solver; callback, reltol, kwargs...)
end

function solve_thermodynamics(model::CosmologicalModel, par::CosmologicalParameters; tini = 1e-5, aend = 1e0, solver = Rodas5P(), reltol = 1e-13, kwargs...) # need very small tolerance to get good csb²
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

function solve_perturbations(model::CosmologicalModel, ks::AbstractArray, par::CosmologicalParameters; tini = 1e-5, aend = 1e0, solver = KenCarp47(), reltol = 1e-9, verbose = false, kwargs...)
    pars = Pair{Any, Any}[ # TODO: avoid Any
        model.th.bg.ph.Ω0 => par.Ωγ0,
        model.th.bg.cdm.Ω0 => par.Ωc0,
        model.th.bg.bar.Ω0 => par.Ωb0,
        model.th.bg.g.h => par.h,
        model.th.bg.neu.Neff => par.Neff,
    ]

    if model.spline_th
        th_sol = solve_thermodynamics(model, par; tini, aend) # TODO: forward kwargs...?
        tend = th_sol[t][end]
        ts = exp.(range(log(tini), log(tend), length=1024)) # TODO: select determine points adaptively from th_sol # TODO: CMB spectrum is sensitive to number of points here!
        push!(pars,
            model.pt_sim.th.rec.dτspline => spline(log.(.-th_sol(ts, idxs=model.th.rec.dτ).u), log.(ts)), # TODO: improve spline accuracy
            #model.pt_sim.th.rec.Tbspline => spline(log.(th_sol[model.th.rec.cs²]), log.(ts)),
            model.pt_sim.th.rec.cs²spline => spline(log.(th_sol(ts, idxs=model.th.rec.Tb).u), log.(ts)),
        )
    end

    ki = 1.0
    prob0 = ODEProblem(model.pt_sim, [], (tini, tend), [pars; k => ki]) # TODO: why do I need this???
    #sol0 = solve(prob0, solver; reltol, kwargs...)
    probs = EnsembleProblem(; safetycopy = false, prob = prob0, prob_func = (prob, i, _) -> begin
        verbose && println("$i/$(length(ks)) k = $(ks[i]*k0) Mpc/h")
        return ODEProblem(model.pt_sim, [], (tini, tend), [pars; k => ks[i]]) # TODO: use remake https://github.com/SciML/OrdinaryDiffEq.jl/pull/2228, https://github.com/SciML/ModelingToolkit.jl/issues/2799 etc. is fixed
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

function solve_perturbations(model::CosmologicalModel, ks::Number, par::CosmologicalParameters; kwargs...)
    return solve_perturbations(model, [ks], par; kwargs...)[1]
end