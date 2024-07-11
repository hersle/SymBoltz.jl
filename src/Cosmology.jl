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

function solve_background(M::CosmologicalModel, par::CosmologicalParameters; tini = 1e-5, aend = 1e0, solver = Rodas5P(), reltol = 1e-13, kwargs...)
    prob = ODEProblem(M.bg_sim, [], (tini, 4.0), [
        M.bg_sim.γ.Ω0 => par.Ωγ0
        M.bg_sim.c.Ω0 => par.Ωc0
        M.bg_sim.b.Ω0 => par.Ωb0
        #M.bg_sim.neu.Neff => par.Neff
        M.bg_sim.g.h => par.h
        M.bg_sim.b.rec.Yp => par.Yp
    ])
    callback = callback_terminator(M.bg_sim, M.bg_sim.g.a, aend)
    return solve(prob, solver; callback, reltol, kwargs...)
end

function solve_perturbations(M::CosmologicalModel, ks::AbstractArray, par::CosmologicalParameters; tini = 1e-5, aend = 1e0, solver = KenCarp47(), reltol = 1e-9, verbose = false, kwargs...)
    pars = Pair{Any, Any}[ # TODO: avoid Any
        M.pt_sim.γ.Ω0 => par.Ωγ0
        M.pt_sim.c.Ω0 => par.Ωc0
        M.pt_sim.b.Ω0 => par.Ωb0
        M.pt_sim.g.h => par.h
        #M.th.bg.neu.Neff => par.Neff,
    ]

    tend = 4.0
    if :b₊rec₊dτspline in Symbol.(parameters(M.pt_sim))
        th_sol = solve_background(M, par; tini, aend) # TODO: forward kwargs...?
        tend = th_sol[t][end]
        ts = exp.(range(log(tini), log(tend), length=1024)) # TODO: select determine points adaptively from th_sol # TODO: CMB spectrum is sensitive to number of points here!
        push!(pars,
            M.pt_sim.b.rec.dτspline => spline(log.(.-th_sol(ts, idxs=M.pt_sim.b.rec.dτ).u), log.(ts)), # TODO: improve spline accuracy
            #M.pt_sim.th.rec.cs²spline => spline(log.(th_sol(ts, idxs=M.th.rec.Tb).u), log.(ts)),
        )
    else
        push!(pars, M.pt_sim.b.rec.Yp => par.Yp)
    end

    ki = 1.0
    prob0 = ODEProblem(M.pt_sim, [], (tini, tend), [pars; k => ki]) # TODO: why do I need this???
    #sol0 = solve(prob0, solver; reltol, kwargs...)
    probs = EnsembleProblem(; safetycopy = false, prob = prob0, prob_func = (prob, i, _) -> begin
        verbose && println("$i/$(length(ks)) k = $(ks[i]*k0) Mpc/h")
        return ODEProblem(M.pt_sim, [], (tini, tend), [pars; k => ks[i]]) # TODO: use remake https://github.com/SciML/OrdinaryDiffEq.jl/pull/2228, https://github.com/SciML/ModelingToolkit.jl/issues/2799 etc. is fixed
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
    return solve(probs, solver, EnsembleThreads(), trajectories = length(ks); reltol, progress=true, kwargs...) # TODO: test GPU parallellization
end

function solve_perturbations(M::CosmologicalModel, ks::Number, par::CosmologicalParameters; kwargs...)
    return solve_perturbations(M, [ks], par; kwargs...)[1]
end