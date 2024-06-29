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

function solve_background(model::CosmologicalModel, par::CosmologicalParameters; aend = 1e0, solver = Vern8(), reltol = 1e-8, kwargs...)
    prob = ODEProblem(model.bg_sim, [], (1e-5, 4.0), [
        model.bg_sim.ph.Ω0 => par.Ωγ0,
        model.bg_sim.cdm.Ω0 => par.Ωc0,
        model.bg_sim.bar.Ω0 => par.Ωb0,
        model.th.bg.g.h => par.h
    ])
    callback = callback_terminator(model.bg_sim, model.bg_sim.g.a, aend)
    return solve(prob, solver; callback, reltol, kwargs...)
end

function solve_thermodynamics(model::CosmologicalModel, par::CosmologicalParameters; aend = NaN, solver = Rodas5P(), reltol = 1e-13, kwargs...) # need very small tolerance to get good csb²
    prob = ODEProblem(model.th_sim, [], (1e-5, 4.0), [
        model.th_sim.bg.ph.Ω0 => par.Ωγ0,
        model.th_sim.bg.cdm.Ω0 => par.Ωc0,
        model.th_sim.bg.bar.Ω0 => par.Ωb0,
        model.th_sim.bg.g.h => par.h,
        model.th_sim.rec.Yp => par.Yp
    ])
    callback = nothing # callback_terminator(model.bg_sim, model.bg_sim.g.a, aend)
    return solve(prob, solver; callback, reltol, kwargs...)
end

function solve_perturbations(model::CosmologicalModel, ks::AbstractArray, par::CosmologicalParameters; solver = KenCarp47(), reltol = 1e-6, verbose = false, kwargs...)
    pars = Pair{Any, Any}[ # TODO: avoid Any
        model.th.bg.ph.Ω0 => par.Ωγ0,
        model.th.bg.cdm.Ω0 => par.Ωc0,
        model.th.bg.bar.Ω0 => par.Ωb0,
        model.th.bg.g.h => par.h,
    ]

    if model.spline_th
        ts = exp.(range(log(1e-5 + 1e-10), log(4.0 - 1e-10), length=1024)) # TODO: select determine points adaptively from th_sol # TODO: CMB spectrum is sensitive to number of points here!
        th_sol = solve_thermodynamics(model, par; saveat = ts) # TODO: forward kwargs...?
        push!(pars,
            model.pt_sim.th.rec.dτspline => CubicSpline(log.(.-th_sol[model.th.rec.dτ]), log.(ts); extrapolate=true), # TODO: improve spline accuracy
            model.pt_sim.th.rec.Tbspline => CubicSpline(log.(th_sol[model.th.rec.cs²]), log.(ts); extrapolate=true),
            model.pt_sim.th.rec.cs²spline => CubicSpline(log.(th_sol[model.th.rec.Tb]), log.(ts); extrapolate=true),
        )
    end

    prob_dummy = ODEProblem(model.pt_sim, [], (1e-5, 4.0), [pars; k => 1.0]) # TODO: why do I need this???
    probs = EnsembleProblem(; safetycopy = false, prob = prob_dummy, prob_func = (_, i, _) -> begin
        verbose && println("$i/$(length(ks)) k = $(ks[i]*k0) Mpc/h")
        return ODEProblem(model.pt_sim, [], (1e-5, 4.0), [pars; k => ks[i]]) # TODO: use remake https://github.com/SciML/OrdinaryDiffEq.jl/pull/2228, https://github.com/SciML/ModelingToolkit.jl/issues/2799 etc. is fixed
    end)
    return solve(probs, solver, EnsembleThreads(), trajectories = length(ks); reltol, progress=true, kwargs...) # TODO: test GPU parallellization
end

function solve_perturbations(model::CosmologicalModel, ks::Number, par::CosmologicalParameters; kwargs...)
    return solve_perturbations(model, [ks], par; kwargs...)[1]
end