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

function perturbations_photon_hierarchy(g, lmax=6, interact=false; name)
    @variables Θ0(η) Θ(η)[1:lmax] δ(η) interactions(η)[1:lmax]
    eqs = [
        Dη(Θ0) + g.k*Θ[1] ~ -Dη(g.Φ)
        Dη(Θ[1]) - g.k/3*(Θ0-2*Θ[2]) ~ g.k/3*g.Ψ + interactions[1]
        [Dη(Θ[l]) ~ g.k/(2*l+1) * (l*Θ[l-1] - (l+1)*Θ[l+1]) + interactions[l] for l in 2:lmax-1]...
        Dη(Θ[lmax]) ~ g.k*Θ[lmax-1] - (lmax+1) * Θ[lmax] / η + interactions[lmax]
        δ ~ 4*Θ0
    ]
    if !interact
        push!(eqs, collect(interactions .~ 0)...)
    end
    return ODESystem(eqs, η; name)
end

# TODO: merge with photon_hierarchy and have polarization flag
function perturbations_polarization_hierarchy(g, lmax=6; name)
    @variables Θ0(η) Θ(η)[1:lmax] dτ(η) Π(η) # TODO: index Θ[l=0] when fixed: https://github.com/SciML/ModelingToolkit.jl/pull/2671
    eqs = [
        Dη(Θ0) + g.k*Θ[1] ~ dτ * (Θ0 - Π/2)
        Dη(Θ[1]) - g.k/(2*1+1) * (1*Θ0 - (1+1)*Θ[1+1]) ~ dτ * (Θ[1] - Π/10*δkron(1,2))
        [Dη(Θ[l]) - g.k/(2*l+1) * (l*Θ[l-1] - (l+1)*Θ[l+1]) ~ dτ * (Θ[l] - Π/10*δkron(l,2)) for l in 2:lmax-1]...
        Dη(Θ[lmax]) ~ g.k*Θ[lmax-1] - (lmax+1) * Θ[lmax] / η + dτ * Θ[lmax]
    ]
    return ODESystem(eqs, η; name)
end

function perturbations_matter(gbg, gpt, interact=false; name)
    @variables δ(η) u(η) Δ(η) interaction(η)
    eqs = [
        Dη(δ) + gpt.k*u ~ -3*Dη(gpt.Φ) # Dodelson (5.69) or (8.12) with i*uc -> uc
        Dη(u) + u*gbg.ℰ ~ gpt.k*gpt.Ψ + interaction # Dodelson (5.70) or (8.13) with i*uc -> uc (opposite sign convention from Hans' website)
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
    @named pol = Symboltz.perturbations_polarization_hierarchy(gpt, lmax)
    @named cdm = Symboltz.perturbations_matter(bg.sys.g, gpt, false)
    @named bar = Symboltz.perturbations_matter(bg.sys.g, gpt, true)
    @named gravpt = Symboltz.perturbations_gravity(bg.sys.g, gpt)
    return Symboltz.PerturbationsSystem(bg, th, gpt, gravpt, ph, pol, cdm, bar; name)
end

# TODO: take list of species, each of which "exposes" contributions to δρ and Π
@register_symbolic dτfunc(η, spl) # TODO: improve somehow
dτfunc(η, spl) = -exp(spl(log(η))) # TODO: type-stable? @code_warntype dτfunc(1e0) gives warnings, but code seems fast?
function PerturbationsSystem(bg::BackgroundSystem, th::ThermodynamicsSystem, g::ODESystem, grav::ODESystem, ph::ODESystem, pol::ODESystem, cdm::ODESystem, bar::ODESystem; name)
    @parameters fb dτspline # TODO: get rid of
    @variables ρc(η) ρb(η) δρr(η) δρc(η) δρb(η) R(η) dτ(η) Δm(η) # TODO: get rid of
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
        bar.interaction     ~ +dτ/R * (bar.u - 3*ph.Θ[1])
        ph.interactions[1] ~ -dτ/3 * (bar.u - 3*ph.Θ[1])
        [ph.interactions[l] ~ dτ * (ph.Θ[l] - ph.Θ[2]/10*δkron(l,2)) for l in 2:lastindex(ph.interactions)]...
        dτ ~ dτfunc(η, dτspline) # TODO: spline over η
    
        # polarization
        pol.dτ ~ dτ
        pol.Π ~ ph.Θ[2] + pol.Θ[2] + pol.Θ0
    
        # gravity shear stress
        grav.Π ~ -32π*bg.sys.g.a^2 * bg.sys.rad.ρ*ph.Θ[2] # TODO: add neutrinos
    ], η, [Δm], [fb, dτspline]; name)
    sys = compose(connections, g, grav, ph, pol, bar, cdm, bg.sys) # TODO: add background stuff?
    ssys = structural_simplify(sys; simplify=true, allow_symbolic=true)
    prob = ODEProblem(ssys, unknowns(ssys) .=> NaN, (0.0, 4.0), parameters(ssys) .=> NaN; jac=true) 
    return PerturbationsSystem(sys, ssys, prob, bg, th)
end

function solve(pt::PerturbationsSystem, kvals::AbstractArray, Ωr0, Ωm0, Ωb0, h, Yp; aini = 1e-8, reltol=1e-8)
    th = pt.th
    bg = th.bg
    fb = Ωb0 / Ωm0; @assert fb <= 1 # TODO: avoid duplication thermo logic
    th_sol = solve(th, Ωr0, Ωm0, Ωb0, h, Yp) # update spline for dτ (e.g. to propagate derivative information through recombination, if called with dual numbers) TODO: use th_sol(a; idxs=th.dτ) directly in a type-stable way?
    ηs = th_sol[η]
    ηini, ηtoday = ηs[begin], ηs[end]
    ρrini, ρmini, ρΛini, Eini = th_sol(ηini; idxs=[bg.sys.rad.ρ, bg.sys.mat.ρ, bg.sys.de.ρ, bg.sys.g.E]) # TODO: avoid duplicate logic
    τs = th_sol[th.sys.τ] .- th_sol[th.sys.τ][end]
    τspline = CubicSpline(τs, ηs)
    dτs = DataInterpolations.derivative.(Ref(τspline), ηs)
    dτspline = CubicSpline(log.(.-dτs), log.(ηs); extrapolate=true) # TODO: extrapolate valid?
    dτini = dτfunc(ηini, dτspline)

    lmax = lastindex(pt.sys.ph.Θ)
    T = eltype([Ωr0, Ωm0, Ωb0, h, Yp])

    function prob_func(_, i, _)
        kval = kvals[i]
        println("$i/$(length(kvals)) k = $(kval*k0) Mpc/h")
    
        Θrini = Vector{T}(undef, lmax) # TODO: avoid Any, use eltype logic? # TODO: allocate outside, but beware of race condition?
        ΘPini = Vector{T}(undef, lmax) # TODO: allow lrmax ≠ lPmax # TODO: avoid Any, use eltype logic?

        Φini = 2/3 # TODO: why 2/3? # arbitrary normalization (from primordial curvature power spectrum?)
        Θrini0 = Φini/2 # Dodelson (7.89)
        Θrini[1] = -kval*Φini/(6*aini*Eini) # Dodelson (7.95)
        Θrini[2] = -8/15#=-20/45=#*kval/dτini * Θrini[1] # TODO: change with/without polarization; another argument for merging photons+polarization
        for l in 3:lmax
            Θrini[l] = -l/(2*l+1) * kval/dτini * Θrini[l-1]
        end
        ΘPini0 = 5/4 * Θrini[2]
        ΘPini[1] = -kval/(4*dτini) * Θrini[2]
        ΘPini[2] = 1/4 * Θrini[2]
        for l in 3:lmax
            ΘPini[l] = -l/(2*l+1) * kval/dτini * ΘPini[l-1]
        end
        δcini = δbini = 3*Θrini0 # Dodelson (7.94)
        ucini = ubini = 3*Θrini[1] # Dodelson (7.95)
        return remake(pt.prob; tspan = (ηini, ηtoday), u0 = Dict(pt.ssys.gpt.Φ => Φini, pt.ssys.ph.Θ0 => Θrini0, [pt.ssys.ph.Θ[l] => Θrini[l] for l in 1:lmax]..., pt.ssys.pol.Θ0 => ΘPini0, [pt.ssys.pol.Θ[l] => ΘPini[l] for l in 1:lmax]..., pt.ssys.bar.δ => δbini, pt.ssys.bar.u => ubini, pt.ssys.cdm.δ => δcini, pt.ssys.cdm.u => ucini, bg.sys.rad.ρ => ρrini, bg.sys.mat.ρ => ρmini, bg.sys.de.ρ => ρΛini, bg.sys.g.a => aini), p = [pt.ssys.fb => fb, pt.ssys.gpt.k => kval, pt.ssys.dτspline => dτspline])
    end

    probs = EnsembleProblem(prob = nothing, prob_func = prob_func)
    return solve(probs, KenCarp4(), EnsembleThreads(), trajectories = length(kvals); reltol) # KenCarp4 and Kvaerno5 seem to work well # TODO: test GPU parallellization
end