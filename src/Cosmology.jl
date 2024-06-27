struct CosmologicalModel
    bg::ODESystem
    pt::ODESystem

    bg_sim::ODESystem
    pt_sim::ODESystem

    spline_th::Bool # TODO: what about a parametric CosmologicalModel?
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

function ΛCDM(; lmax=6, spline_th = true)
    @named bg = background_ΛCDM(; thermo=false) # without thermodynamics; just for use in constructing perturbations
    @named pt = perturbations_ΛCDM(bg, lmax; spline_th) # TODO: spline full background OR select quantities
    @named bg = background_ΛCDM(; thermo=true)
    bg_sim = structural_simplify(bg)
    pt_sim = structural_simplify(pt)
    return CosmologicalModel(bg, pt, bg_sim, pt_sim, spline_th)
end

# TODO: add aini, aend

function solve_background(model::CosmologicalModel, par::CosmologicalParameters; aend = NaN, solver = Rodas5P(), reltol = 1e-6, kwargs...)
    prob = ODEProblem(model.bg_sim, [], (1e-5, 4.0), [
        model.bg_sim.ph.Ω0 => par.Ωγ0,
        model.bg_sim.neu.Ω0 => par.Ων0,
        model.bg_sim.cdm.Ω0 => par.Ωc0,
        model.bg_sim.bar.Ω0 => par.Ωb0,
        model.bg_sim.g.h => par.h,
        model.bg_sim.th.Yp => par.Yp
    ])
    callback = callback_terminator(model.bg_sim, model.bg_sim.g.a, aend)
    return solve(prob, solver; callback, reltol, kwargs...)
end

function solve_perturbations(model::CosmologicalModel, ks::AbstractArray, par::CosmologicalParameters; solver = KenCarp47(), reltol = 1e-6, verbose = false, kwargs...)
    pars = Pair{Any, Any}[ # TODO: avoid Any
        model.bg.ph.Ω0 => par.Ωγ0,
        model.bg.neu.Ω0 => par.Ων0,
        model.bg.cdm.Ω0 => par.Ωc0,
        model.bg.bar.Ω0 => par.Ωb0,
        model.bg.g.h => par.h,
        model.bg.th.Yp => par.Yp
    ]

    if model.spline_th
        ts = exp.(range(log(1e-5 + 1e-10), log(4.0 - 1e-10), length=1024)) # TODO: select determine points adaptively from th_sol # TODO: CMB spectrum is sensitive to number of points here!
        bg_sol = solve_background(model, par; saveat = ts) # TODO: forward kwargs...?
        dτs = bg_sol[model.bg.th.dτ] # TODO: interpolate directly from ODESolution?
        dτspline = CubicSpline(log.(.-dτs), log.(ts); extrapolate=true) # spline this logarithmically for accuracy during integration
        push!(pars, model.pt_sim.th.dτspline => dτspline)
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