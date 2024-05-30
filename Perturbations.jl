struct PerturbationsSystem
    sys::ODESystem
    ssys::ODESystem
    prob::ODEProblem
    bg::BackgroundSystem
    th::ThermodynamicsSystem
end

# TODO: merge with background metric!
function perturbations_metric(; name)
    @variables Φ(η) Ψ(η)
    @parameters k # perturbation wavenumber
    k, Φ, Ψ = GlobalScope.([k, Φ, Ψ])
    return ODESystem(Equation[], η, [Φ, Ψ], [k]; name)
end

function perturbations_radiation(g, interact=false; name)
    @variables Θ0(η) Θ1(η) δ(η)
    interaction = interact ? only(@variables interaction(η)) : 0
    return ODESystem([
        Dη(Θ0) + k*Θ1 ~ -Dη(g.Φ) # Dodelson (5.67) or (8.10)
        Dη(Θ1) - k/3*Θ0 ~ k/3*g.Ψ + interaction # Dodelson (5.67) or (8.11)
        δ ~ 4*Θ0
    ], η; name)
end

function perturbations_photon_hierarchy(g, lmax=6, polarization=true; name)
    @variables Θ0(η) Θ(η)[1:lmax] δ(η) dτ(η) ub(η)
    eqs = [
        Dη(Θ0) + g.k*Θ[1] ~ -Dη(g.Φ)
        Dη(Θ[1]) - g.k/3*(Θ0-2*Θ[2]) ~ g.k/3*g.Ψ - dτ/3 * (ub - 3*Θ[1])
        [Dη(Θ[l]) ~ g.k/(2*l+1) * (l*Θ[l-1] - (l+1)*Θ[l+1]) + dτ * (Θ[l] - Θ[2]/10*δkron(l,2)) for l in 2:lmax-1]...
        Dη(Θ[lmax]) ~ g.k*Θ[lmax-1] - (lmax+1) * Θ[lmax] / η + dτ * (Θ[lmax] - Θ[2]/10*δkron(lmax,2))
        δ ~ 4*Θ0
    ]
    if polarization
        @variables ΘP0(η) ΘP(η)[1:lmax] Π(η)
        push!(eqs,
            Dη(ΘP0) + g.k*ΘP[1] ~ dτ * (ΘP0 - Π/2),
            Dη(ΘP[1]) - g.k/(2*1+1) * (1*ΘP0 - (1+1)*ΘP[1+1]) ~ dτ * (ΘP[1] - Π/10*δkron(1,2)),
            [Dη(ΘP[l]) - g.k/(2*l+1) * (l*ΘP[l-1] - (l+1)*ΘP[l+1]) ~ dτ * (ΘP[l] - Π/10*δkron(l,2)) for l in 2:lmax-1]...,
            Dη(ΘP[lmax]) ~ g.k*ΘP[lmax-1] - (lmax+1) * ΘP[lmax] / η + dτ * ΘP[lmax],
            Π ~ Θ[2] + ΘP[2] + ΘP0
        )
    end
    return ODESystem(eqs, η; name)
end

function perturbations_matter(gbg, gpt, interact=false; name)
    @variables δ(η) u(η) Δ(η) interaction(η)
    eqs = [
        Dη(δ) + gpt.k*u ~ -3*Dη(gpt.Φ) # Dodelson (5.69) or (8.12) with i*uc -> uc
        Dη(u) + u*gbg.ℰ ~ gpt.k*gpt.Ψ + interaction # Dodelson (5.70) or (8.13) with i*uc -> uc (opposite sign convention from Hans' website) # TODO: treat interaction explicitly
        Δ ~ δ + 3 * gbg.ℰ/gpt.k * u # TODO: correct??? gauge-invariant density perturbation (https://arxiv.org/pdf/1307.1459#equation.2.12)
    ]
    if !interact
        push!(eqs, interaction ~ 0)
    end
    return ODESystem(eqs, η, [δ, u, Δ, interaction], []; name)
end

function perturbations_gravity(gbg, gpt; name)
    @variables δρ(η) Π(η)
    return ODESystem([
        Dη(gpt.Φ) ~ (4π*gbg.a^2*δρ - gpt.k^2*gpt.Φ + 3*gbg.ℰ^2*gpt.Ψ) / (3*gbg.ℰ) # Dodelson (6.41) # TODO: write in more natural form?
        gpt.k^2 * (gpt.Ψ + gpt.Φ) ~ Π # anisotropic stress
    ], η; name)
end

function perturbations_ΛCDM(th::ThermodynamicsSystem, lmax::Int; name)
    bg = th.bg
    @named gpt = Symboltz.perturbations_metric()
    @named ph = Symboltz.perturbations_photon_hierarchy(gpt, lmax, true)
    @named cdm = Symboltz.perturbations_matter(bg.sys.g, gpt, false)
    @named bar = Symboltz.perturbations_matter(bg.sys.g, gpt, true)
    @named gravpt = Symboltz.perturbations_gravity(bg.sys.g, gpt)
    return Symboltz.PerturbationsSystem(bg, th, gpt, gravpt, ph, cdm, bar; name)
end

# TODO: take list of species, each of which "exposes" contributions to δρ and Π
@register_symbolic dτfunc(η, spl) # TODO: improve somehow
dτfunc(η, spl) = -exp(spl(log(η))) # TODO: type-stable? @code_warntype dτfunc(1e0) gives warnings, but code seems fast?
function PerturbationsSystem(bg::BackgroundSystem, th::ThermodynamicsSystem, g::ODESystem, grav::ODESystem, ph::ODESystem, cdm::ODESystem, bar::ODESystem; name)
    @parameters fb dτspline # TODO: get rid of
    @variables ρc(η) ρb(η) δρr(η) δρc(η) δρb(η) R(η) dτ(η) Δm(η) # TODO: get rid of
    lmax = lastindex(ph.Θ)
    defaults = [
        g.Φ => 2/3 # TODO: why 2/3? # arbitrary normalization (from primordial curvature power spectrum?)
        ph.Θ0 => 1/2 * g.Φ # Dodelson (7.89)
        ph.Θ[1] => -1/6 * g.k/bg.sys.g.ℰ * g.Φ # Dodelson (7.95)
        ph.Θ[2] => -8/15#=-20/45=#*g.k/dτ * ph.Θ[1] # TODO: change with/without polarization; another argument for merging photons+polarization
        [ph.Θ[l] => -l/(2*l+1) * g.k/dτ * ph.Θ[l-1] for l in 3:lmax]...
        ph.ΘP0 => 5/4 * ph.Θ[2]
        ph.ΘP[1] => -1/4 * g.k/dτ * ph.Θ[2]
        ph.ΘP[2] => 1/4 * ph.Θ[2]
        [ph.ΘP[l] => -l/(2*l+1) * g.k/dτ * ph.ΘP[l-1] for l in 3:lmax]...
        cdm.δ => 3 * ph.Θ0 # Dodelson (7.94)
        bar.δ => 3 * ph.Θ0
        cdm.u => 3 * ph.Θ[1] # Dodelson (7.95)
        bar.u => 3 * ph.Θ[1]
    ]
    connections = ODESystem([
        ρb ~ fb * bg.sys.mat.ρ
        ρc ~ bg.sys.mat.ρ - ρb
        δρr ~ ph.δ * bg.sys.rad.ρ
        δρc ~ cdm.δ * ρc
        δρb ~ bar.δ * ρb
        grav.δρ ~ δρr + δρc + δρb # total energy density perturbation

        Δm ~ (ρc * cdm.Δ + ρb * bar.Δ) / bg.sys.mat.ρ # total gauge-invariant matter overdensity
    
        # baryon-photon interactions: Compton (Thomson) scattering # TODO: define connector type?
        R ~ 3/4 * ρb / bg.sys.rad.ρ # Dodelson (5.74)
        bar.interaction ~ +dτ/R * (bar.u - 3*ph.Θ[1])
        ph.ub ~ bar.u
        ph.dτ ~ dτ
        dτ ~ dτfunc(η, dτspline) # TODO: spline over η
    
        # gravity shear stress
        grav.Π ~ -32π*bg.sys.g.a^2 * bg.sys.rad.ρ*ph.Θ[2] # TODO: add neutrinos
    ], η, [Δm, dτ], [fb, dτspline]; name, defaults)
    sys = compose(connections, g, grav, ph, bar, cdm, bg.sys) # TODO: add background stuff?
    ssys = structural_simplify(sys)
    prob = ODEProblem(ssys, unknowns(ssys) .=> NaN, (0.0, 4.0), parameters(ssys) .=> NaN; jac=true) 
    return PerturbationsSystem(sys, ssys, prob, bg, th)
end

function solve(pt::PerturbationsSystem, ks::AbstractArray, Ωr0, Ωm0, Ωb0, h, Yp; aini = 1e-8, reltol=1e-8)
    th = pt.th
    bg = th.bg
    fb = Ωb0 / Ωm0; @assert fb <= 1 # TODO: avoid duplication thermo logic
    th_sol = solve(th, Ωr0, Ωm0, Ωb0, h, Yp) # update spline for dτ (e.g. to propagate derivative information through recombination, if called with dual numbers) TODO: use th_sol(a; idxs=th.dτ) directly in a type-stable way?
    ηs = th_sol[η]
    ηini, ηtoday = ηs[begin], ηs[end]
    ρr, ρm, ρΛ, a, ℰ = th_sol(ηini; idxs=[bg.sys.rad.ρ, bg.sys.mat.ρ, bg.sys.de.ρ, bg.sys.g.a, bg.sys.g.ℰ]) # TODO: avoid duplicate logic
    τs = th_sol[th.sys.τ] .- th_sol[th.sys.τ][end]
    τspline = CubicSpline(τs, ηs)
    dτs = DataInterpolations.derivative.(Ref(τspline), ηs)
    dτspline = CubicSpline(log.(.-dτs), log.(ηs); extrapolate=true) # TODO: extrapolate valid?
    dτ = dτs[begin]

    p = Dict(pt.ssys.fb => fb, pt.ssys.gpt.k => 0.0, bg.sys.g.H0 => NaN, pt.ssys.dτspline => NaN)
    u0 = Dict(bg.sys.rad.ρ => ρr, bg.sys.mat.ρ => ρm, bg.sys.de.ρ => ρΛ, bg.sys.g.a => a, bg.sys.g.ℰ => ℰ, pt.ssys.dτ => dτ) # merge(ModelingToolkit.defaults(pt.ssys), Dict(pt.ssys.dτ => dτ, bg.sys.g.ℰ => ℰ, bg.sys.rad.ρ => ρr, bg.sys.mat.ρ => ρm, bg.sys.de.ρ => ρΛ, bg.sys.g.a => a))    
    prob = ODEProblem(pt.ssys, u0, (ηini, ηtoday), p; jac = true)
    probs = EnsembleProblem(; prob, prob_func = (prob, i, _) -> begin
        k = ks[i]
        println("$i/$(length(ks)) k = $(k*k0) Mpc/h")
        return remake(prob; p = [pt.ssys.gpt.k => k, pt.ssys.dτspline => dτspline])
    end)

    return solve(probs, KenCarp4(), EnsembleThreads(), trajectories = length(ks); reltol) # KenCarp4 and Kvaerno5 seem to work well # TODO: test GPU parallellization
end