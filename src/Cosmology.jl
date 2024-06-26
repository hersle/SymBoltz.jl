struct CosmologicalModel
    bg::ODESystem
    th::ODESystem
    pt::ODESystem
    bg_sim::ODESystem
    th_sim::ODESystem
    pt_sim::ODESystem
end

@kwdef struct CosmologicalParameters
    Ωγ0 = 5.5e-5
    Ων0 = 3.046 * 7/8 * (4/11)^(4/3) * 5.5e-5 # TODO: handle more elegantly with Neff/Tν0
    Ωc0 = 0.267
    Ωb0 = 0.05
    h = 0.67
    As = 2.1e-9
    Yp = 0.245
end

function ΛCDM(; lmax=6)
    @named bg = background_ΛCDM()
    @named th = thermodynamics_ΛCDM(bg)
    @named pt = perturbations_ΛCDM(th, lmax)

    bg_sim = structural_simplify(bg)
    th_sim = structural_simplify(th)
    pt_sim = structural_simplify(pt)

    return CosmologicalModel(bg, th, pt, bg_sim, th_sim, pt_sim)
end

# TODO: add aini, aend

function solve_background(model::CosmologicalModel, par::CosmologicalParameters; aend = 1e0, solver = Vern8(), reltol = 1e-8, kwargs...)
    prob = ODEProblem(model.bg_sim, [], (1e-5, 4.0), [
        model.bg_sim.ph.Ω0 => par.Ωγ0,
        model.bg_sim.neu.Ω0 => par.Ων0,
        model.bg_sim.cdm.Ω0 => par.Ωc0,
        model.bg_sim.bar.Ω0 => par.Ωb0,
        model.bg_sim.g.h => NaN
    ])
    aindex = variable_index(model.bg_sim, model.bg_sim.g.a)
    callback = ContinuousCallback((u, _, _) -> (a = u[aindex]; a - aend), terminate!) # stop when a == 1 today
    return solve(prob, solver; callback, reltol, kwargs...)
end

function solve_thermodynamics(model::CosmologicalModel, par::CosmologicalParameters; solver = Rodas5P(), reltol = 1e-6, kwargs...)
    prob = ODEProblem(model.th_sim, [], (1e-5, 4.0), [
        model.bg.ph.Ω0 => par.Ωγ0,
        model.bg.neu.Ω0 => par.Ων0,
        model.bg.cdm.Ω0 => par.Ωc0,
        model.bg.bar.Ω0 => par.Ωb0,
        model.bg.g.h => par.h,
        model.th_sim.Yp => par.Yp
    ])
    return solve(prob, solver; reltol, kwargs...)
end

function solve_perturbations(model::CosmologicalModel, ks::AbstractArray, par::CosmologicalParameters; solver = Rodas5P(), reltol = 1e-6, verbose = false, kwargs...)
    probs = EnsembleProblem(; safetycopy = false, prob = nothing, prob_func = (_, i, _) -> begin
        verbose && println("$i/$(length(ks)) k = $(ks[i]*k0) Mpc/h")
        return ODEProblem(model.pt_sim, [], (1e-5, 4.0), [
            model.th.bg.ph.Ω0 => par.Ωγ0,
            model.th.bg.neu.Ω0 => par.Ων0,
            model.th.bg.cdm.Ω0 => par.Ωc0,
            model.th.bg.bar.Ω0 => par.Ωb0,
            model.th.bg.g.h => par.h,
            model.th.Yp => par.Yp,
            k => ks[i]
        ])
    end)
    return solve(probs, solver, EnsembleSerial(), trajectories = length(ks); reltol, kwargs...) # TODO: test GPU parallellization
end

function solve_perturbations(model::CosmologicalModel, ks::Number, par::CosmologicalParameters; kwargs...)
    return solve_perturbations(model, [ks], par; kwargs...)[1]
end