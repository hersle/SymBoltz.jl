struct PerturbationsSystem
    sys::ODESystem
    ssys::ODESystem
    bg::BackgroundSystem
    th::ThermodynamicsSystem
end

# TODO: merge with background metric!
function perturbations_metric(; name)
    @variables Φ(η) Ψ(η)
    @parameters k # perturbation wavenumber
    k, Φ, Ψ = GlobalScope.([k, Φ, Ψ])
    defaults = [Φ => 2/3] # TODO: why 2/3? # arbitrary normalization (from primordial curvature power spectrum?)
    return ODESystem(Equation[], η, [Φ, Ψ], [k]; defaults, name)
end

function perturbations_species(gbg, gpt, w, cs² = w, w′ = 0, σ = 0; uinteract = false, name)
    @assert w′ == 0 && σ == 0 # TODO: relax (need to include in ICs)
    @variables δ(η) u(η) Δ(η) uinteraction(η)
    eqs = [
        Dη(δ) ~ -(1+w)*(gpt.k*u+3*Dη(gpt.Φ)) - 3*gbg.ℰ*(cs²-w)*δ # Bertschinger & Ma (30) with Φ -> -Φ; or Baumann (4.4.173) with Φ -> -Φ
        Dη(u) ~ -gbg.ℰ*(1-3*w)*u - w′/(1+w)*u + cs²/(1+w)*gpt.k*δ + gpt.k*gpt.Ψ - gpt.k*σ + uinteraction # Bertschinger & Ma (30) with θ = kv
        Δ ~ δ + 3 * (1+w) * gbg.ℰ/gpt.k * u # Baumann (4.2.144) with v -> -u
    ]
    defaults = [
        δ => 3/2 * (1+w) * gpt.Φ # adiabatic: δᵢ/(1+wᵢ) == δⱼ/(1+wⱼ) (https://cmb.wintherscoming.no/theory_initial.php#adiabatic)
        u => -1/2 * gpt.k / gbg.ℰ * gpt.Φ # TODO: include σ ≠ 0 # solve u′ + ℋ(1-3w)u = w/(1+w)*kδ + kΨ with Ψ=const, IC for δ, Φ=-Ψ, ℋ=H₀√(Ωᵣ₀)/a after converting ′ -> d/da by gathering terms with u′ and u in one derivative using the trick to multiply by exp(X(a)) such that X′(a) will "match" the terms in front of u
    ]
    if !uinteract
        push!(eqs, uinteraction ~ 0)
    end
    return ODESystem(eqs, η, [δ, u, Δ, uinteraction], []; defaults, name)
end

perturbations_matter(gbg, gpt; kwargs...) = perturbations_species(gbg, gpt, 0; kwargs...)
perturbations_radiation(gbg, gpt; kwargs...) = perturbations_species(gbg, gpt, 0; kwargs...)
perturbations_cosmological_constant(gbg, gpt; kwargs...) = perturbations_species(gbg, gpt, 0; kwargs...) # TODO: ill-defined?

function perturbations_photon_hierarchy(gbg, gpt, lmax=6, polarization=true; name)
    @variables Θ0(η) Θ(η)[1:lmax] δ(η) dτ(η) ub(η) Π(η)
    eqs = [
        Dη(Θ0) + gpt.k*Θ[1] ~ -Dη(gpt.Φ)
        Dη(Θ[1]) - gpt.k/3*(Θ0-2*Θ[2]) ~ gpt.k/3*gpt.Ψ - dτ/3 * (ub - 3*Θ[1])
        [Dη(Θ[l]) ~ gpt.k/(2*l+1) * (l*Θ[l-1] - (l+1)*Θ[l+1]) + dτ * (Θ[l] - Π/10*δkron(l,2)) for l in 2:lmax-1]... # TODO: Π in last term here?
        Dη(Θ[lmax]) ~ gpt.k*Θ[lmax-1] - (lmax+1) * Θ[lmax] / η + dτ * (Θ[lmax] - Π/10*δkron(lmax,2))
        δ ~ 4*Θ0
    ]
    defaults = Dict(
        Θ0 => 1/2 * gpt.Φ, # Dodelson (7.89)
        Θ[1] => -1/6 * gpt.k/gbg.ℰ * gpt.Φ, # Dodelson (7.95)
        # Θ[2] is set differently depending on polarization below
        [Θ[l] => -l/(2*l+1) * gpt.k/dτ * Θ[l-1] for l in 3:lmax]...
    )
    if polarization
        @variables ΘP0(η) ΘP(η)[1:lmax]
        push!(eqs,
            Dη(ΘP0) + gpt.k*ΘP[1] ~ dτ * (ΘP0 - Π/2),
            Dη(ΘP[1]) - gpt.k/(2*1+1) * (1*ΘP0 - (1+1)*ΘP[1+1]) ~ dτ * (ΘP[1] - Π/10*δkron(1,2)),
            [Dη(ΘP[l]) - gpt.k/(2*l+1) * (l*ΘP[l-1] - (l+1)*ΘP[l+1]) ~ dτ * (ΘP[l] - Π/10*δkron(l,2)) for l in 2:lmax-1]...,
            Dη(ΘP[lmax]) ~ gpt.k*ΘP[lmax-1] - (lmax+1) * ΘP[lmax] / η + dτ * ΘP[lmax],
            Π ~ Θ[2] + ΘP[2] + ΘP0
        )
        push!(defaults,
            Θ[2] => -8/15 * gpt.k/dτ * Θ[1], # with polarization, overwrite from above
            ΘP0 => 5/4 * Θ[2],
            ΘP[1] => -1/4 * gpt.k/dτ * Θ[2],
            ΘP[2] => 1/4 * Θ[2],
            [ΘP[l] => -l/(2*l+1) * gpt.k/dτ * ΘP[l-1] for l in 3:lmax]...
        )
    else
        push!(eqs, Π ~ Θ[2])
        push!(defaults, Θ[2] => -20/45 * gpt.k/dτ * Θ[1]) # without polarization
    end
    return ODESystem(eqs, η; defaults, name)
end

function perturbations_gravity(gbg, gpt; name)
    @variables δρ(η) Π(η)
    return ODESystem([
        Dη(gpt.Φ) ~ 4π*gbg.a^2*δρ/(3*gbg.ℰ) - gpt.k^2*gpt.Φ/(3*gbg.ℰ) + gbg.ℰ*gpt.Ψ  # (4π*gbg.a^2*δρ - gpt.k^2*gpt.Φ + 3*gbg.ℰ^2*gpt.Ψ) / (3*gbg.ℰ) # Dodelson (6.41) # TODO: write in more natural form?
        gpt.k^2 * (gpt.Ψ + gpt.Φ) ~ Π # anisotropic stress
    ], η; name)
end

function perturbations_ΛCDM(th::ThermodynamicsSystem, lmax::Int; name)
    bg = th.bg
    @named gpt = Symboltz.perturbations_metric()
    @named ph = Symboltz.perturbations_photon_hierarchy(bg.sys.g, gpt, lmax, true)
    @named cdm = Symboltz.perturbations_matter(bg.sys.g, gpt; uinteract=false)
    @named bar = Symboltz.perturbations_matter(bg.sys.g, gpt; uinteract=true)
    @named gravpt = Symboltz.perturbations_gravity(bg.sys.g, gpt)
    return Symboltz.PerturbationsSystem(bg, th, gpt, gravpt, ph, cdm, bar; name)
end

# TODO: take list of species, each of which "exposes" contributions to δρ and Π
function PerturbationsSystem(bg::BackgroundSystem, th::ThermodynamicsSystem, g::ODESystem, grav::ODESystem, ph::ODESystem, cdm::ODESystem, bar::ODESystem; name)
    @parameters fb τspline::CubicSpline dτspline::CubicSpline ddτspline::CubicSpline csb²spline::CubicSpline # TODO: get rid of
    @variables ρc(η) ρb(η) δρr(η) δρc(η) δρb(η) R(η) csb²(η) τ(η) dτ(η) ddτ(η) v(η) dv(η) Δm(η) SSW(η) SD(η) SISW(η) S(η) # TODO: get rid of
    # TODO: do various IC types (adiabatic, isocurvature, ...) from here?
    connections = ODESystem([
        # gravity density and shear stress
        ρb ~ fb * bg.sys.mat.ρ
        ρc ~ bg.sys.mat.ρ - ρb
        δρr ~ ph.δ * bg.sys.rad.ρ
        δρc ~ cdm.δ * ρc
        δρb ~ bar.δ * ρb
        grav.δρ ~ δρr + δρc + δρb # total energy density perturbation
        grav.Π ~ -32π*bg.sys.g.a^2 * bg.sys.rad.ρ*ph.Θ[2] # TODO: add neutrinos
    
        # baryon-photon interactions: Compton (Thomson) scattering # TODO: define connector type?
        R ~ 3/4 * ρb / bg.sys.rad.ρ # Dodelson (5.74)
        bar.uinteraction ~ #=g.k*csb²*bar.δ +=# dτ/R * (bar.u - 3*ph.Θ[1]) # TODO: enable csb² when it seems stable...
        ph.ub ~ bar.u
        ph.dτ ~ dτ

        # gauge-independent matter overdensity for matter power spectrum
        Δm ~ (ρc * cdm.Δ + ρb * bar.Δ) / bg.sys.mat.ρ

        # source functions for CMB power spectrum
        τ ~ spleval(η, τspline)
        dτ ~ -exp(spleval(log(η), dτspline))
        ddτ ~ spleval(η, ddτspline)
        csb² ~ spleval(log(η), csb²spline)
        v ~ -dτ * exp(-τ) # visibility function
        dv ~ (dτ^2 - ddτ) * exp(-τ)
        SSW ~ v * (ph.Θ0 + g.Ψ + ph.Π/4) # SW
        SD ~ (dv * bar.u + v * Dη(bar.u)) / g.k # Doppler
        SISW ~ exp(-τ) * (Dη(g.Ψ) - Dη(g.Φ)) # ISW
        S ~ SSW + SD + SISW # TODO: add polarization
    ], η, [Δm, csb², τ, dτ, ddτ, v, dv, SSW, SD, SISW, S], [fb, csb²spline, τspline, dτspline, ddτspline]; name)
    sys = compose(connections, g, grav, ph, bar, cdm, bg.sys) # TODO: add background stuff?
    ssys = structural_simplify(sys)
    return PerturbationsSystem(sys, ssys, bg, th)
end

# TODO: get rid of allocations!!! use SVector somehow? see https://docs.sciml.ai/ModelingToolkit/dev/basics/FAQ/#Change-the-unknown-variable-vector-type, also follow 
function solve(pt::PerturbationsSystem, ks::AbstractArray, Ωr0, Ωm0, Ωb0, h, Yp; solver = KenCarp4(), aini = 1e-8, reltol = 1e-10, verbose = false, kwargs...)
    th = pt.th
    bg = th.bg
    fb = Ωb0 / Ωm0; @assert fb <= 1 # TODO: avoid duplication thermo logic
    th_sol = solve(th, Ωr0, Ωm0, Ωb0, h, Yp) # update spline for dτ (e.g. to propagate derivative information through recombination, if called with dual numbers) TODO: use th_sol(a; idxs=th.dτ) directly in a type-stable way?
    ηs = th_sol[η]
    ηini, ηtoday = ηs[begin], ηs[end]
    ρr, ρm, ρΛ, a, ℰ = th_sol(ηini; idxs=[bg.sys.rad.ρ, bg.sys.mat.ρ, bg.sys.de.ρ, bg.sys.g.a, bg.sys.g.ℰ]) # TODO: avoid duplicate logic
    τs = th_sol[th.sys.τ] .- th_sol[th.sys.τ][end]
    τspline = CubicSpline(τs, ηs; extrapolate=true)
    dτs = DataInterpolations.derivative.(Ref(τspline), ηs)
    dτspline = CubicSpline(log.(.-dτs), log.(ηs); extrapolate=true) # spline this logarithmically for accurayc during integration # TODO: extrapolate valid?
    dτ = dτs[begin]
    ddτs = DataInterpolations.derivative.(Ref(τspline), ηs, 2) # spline normally (just observed anyway)
    ddτspline = CubicSpline(ddτs, ηs; extrapolate=true)

    csb²s = th_sol[th.sys.cs²]
    csb²spline = CubicSpline(csb²s, log.(ηs); extrapolate=true)

    p = [pt.ssys.fb => fb, pt.ssys.gpt.k => 0.0, bg.sys.g.H0 => NaN, pt.ssys.τspline => τspline, pt.ssys.dτspline => dτspline, pt.ssys.ddτspline => ddτspline, pt.ssys.csb²spline => csb²spline]
    u0 = [bg.sys.rad.ρ => ρr, bg.sys.mat.ρ => ρm, bg.sys.de.ρ => ρΛ, bg.sys.g.a => a, bg.sys.g.ℰ => ℰ, pt.ssys.dτ => dτ]
    prob = ODEProblem(pt.ssys, u0, (ηini, ηtoday), p; jac = true, sparse = false) # TODO: sparse fails with dual numbers
    probs = EnsembleProblem(; safetycopy = false, prob, prob_func = (prob, i, _) -> begin
        k = ks[i]
        if verbose
            println("$i/$(length(ks)) k = $(k*k0) Mpc/h")
        end
        return remake(prob; p = [pt.ssys.gpt.k => k])
    end)

    return solve(probs, solver, EnsembleThreads(), trajectories = length(ks); reltol, kwargs...) # KenCarp4 and Kvaerno5 seem to work well # TODO: test GPU parallellization
end

# proxy function for evaluating a spline
@register_symbolic spleval(x, spline::CubicSpline)
spleval(x, spline) = spline(x)