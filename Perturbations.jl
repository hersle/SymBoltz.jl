struct PerturbationsSystem
    sys::ODESystem
    ssys::ODESystem
    bg::BackgroundSystem
    th::ThermodynamicsSystem
end

# TODO: merge with background metric!
function perturbations_metric(; kwargs...)
    Φ, Ψ = GlobalScope.(@variables Φ(η) Ψ(η))
    k, = GlobalScope.(@parameters k) # perturbation wavenumber
    return ODESystem(Equation[], η, [Φ, Ψ], [k]; defaults = [Φ => 2/3], kwargs...)
end

function perturbations_species(g0, g1, w, cs² = w, w′ = 0, σ = 0; uinteract = false, kwargs...)
    @assert w′ == 0 && σ == 0 # TODO: relax (need to include in ICs)
    @variables δ(η) u(η) Δ(η) uinteraction(η)
    eqs = [
        Dη(δ) ~ -(1+w)*(g1.k*u+3*Dη(g1.Φ)) - 3*g0.ℰ*(cs²-w)*δ # Bertschinger & Ma (30) with Φ -> -Φ; or Baumann (4.4.173) with Φ -> -Φ
        Dη(u) ~ -g0.ℰ*(1-3*w)*u - w′/(1+w)*u + cs²/(1+w)*g1.k*δ + g1.k*g1.Ψ - g1.k*σ + uinteraction # Bertschinger & Ma (30) with θ = kv
        Δ ~ δ + 3 * (1+w) * g0.ℰ/g1.k * u # Baumann (4.2.144) with v -> -u
    ]
    defaults = [
        δ => 3/2 * (1+w) * g1.Φ # adiabatic: δᵢ/(1+wᵢ) == δⱼ/(1+wⱼ) (https://cmb.wintherscoming.no/theory_initial.php#adiabatic)
        u => -1/2 * g1.k / g0.ℰ * g1.Φ # TODO: include σ ≠ 0 # solve u′ + ℋ(1-3w)u = w/(1+w)*kδ + kΨ with Ψ=const, IC for δ, Φ=-Ψ, ℋ=H₀√(Ωᵣ₀)/a after converting ′ -> d/da by gathering terms with u′ and u in one derivative using the trick to multiply by exp(X(a)) such that X′(a) will "match" the terms in front of u
    ]
    !uinteract && push!(eqs, uinteraction ~ 0)
    return ODESystem(eqs, η; defaults, kwargs...)
end

perturbations_matter(g0, g1; kwargs...) = perturbations_species(g0, g1, 0; kwargs...)
perturbations_radiation(g0, g1; kwargs...) = perturbations_species(g0, g1, 1//3; kwargs...)
perturbations_cosmological_constant(g0, g1; kwargs...) = perturbations_species(g0, g1, -1; kwargs...) # TODO: ill-defined?

function perturbations_photon_hierarchy(g0, g1, lmax=6, polarization=true; kwargs...)
    @variables Θ0(η) Θ(η)[1:lmax] δ(η) dτ(η) ub(η) Π(η) ΘP0(η) ΘP(η)[1:lmax]
    eqs = [
        Dη(Θ0) + g1.k*Θ[1] ~ -Dη(g1.Φ)
        Dη(Θ[1]) - g1.k/3*(Θ0-2*Θ[2]) ~ g1.k/3*g1.Ψ - dτ/3 * (ub - 3*Θ[1])
        [Dη(Θ[l]) ~ g1.k/(2*l+1) * (l*Θ[l-1] - (l+1)*Θ[l+1]) + dτ * (Θ[l] - Π/10*δkron(l,2)) for l in 2:lmax-1]... # TODO: Π in last term here?
        Dη(Θ[lmax]) ~ g1.k*Θ[lmax-1] - (lmax+1) * Θ[lmax] / η + dτ * (Θ[lmax] - Π/10*δkron(lmax,2))
        δ ~ 4*Θ0
        Π ~ Θ[2] + ΘP[2] + ΘP0
        (polarization ? [
            Dη(ΘP0) + g1.k*ΘP[1] ~ dτ * (ΘP0 - Π/2)
            Dη(ΘP[1]) - g1.k/(2*1+1) * (1*ΘP0 - (1+1)*ΘP[1+1]) ~ dτ * (ΘP[1] - Π/10*δkron(1,2))
            [Dη(ΘP[l]) - g1.k/(2*l+1) * (l*ΘP[l-1] - (l+1)*ΘP[l+1]) ~ dτ * (ΘP[l] - Π/10*δkron(l,2)) for l in 2:lmax-1]...
            Dη(ΘP[lmax]) ~ g1.k*ΘP[lmax-1] - (lmax+1) * ΘP[lmax] / η + dτ * ΘP[lmax]
        ] : [
            ΘP0 ~ 0, collect(ΘP .~ 0)... # pin to zero
        ])...
    ]
    defaults = [
        Θ0 => 1/2 * g1.Φ, # Dodelson (7.89)
        Θ[1] => -1/6 * g1.k/g0.ℰ * g1.Φ, # Dodelson (7.95)
        Θ[2] => (polarization ? -8/15 : -20/45) * g1.k/dτ * Θ[1], # # depends on whether polarization is included
        [Θ[l] => -l/(2*l+1) * g1.k/dτ * Θ[l-1] for l in 3:lmax]...,
        ΘP0 => 5/4 * Θ[2],
        ΘP[1] => -1/4 * g1.k/dτ * Θ[2],
        ΘP[2] => 1/4 * Θ[2],
        [ΘP[l] => -l/(2*l+1) * g1.k/dτ * ΘP[l-1] for l in 3:lmax]...
    ]
    return ODESystem(eqs, η; defaults, kwargs...)
end

function perturbations_gravity(g0, g1; kwargs...)
    @variables δρ(η) Π(η)
    return ODESystem([
        Dη(g1.Φ) ~ 4π*g0.a^2*δρ/(3*g0.ℰ) - g1.k^2*g1.Φ/(3*g0.ℰ) + g0.ℰ*g1.Ψ
        g1.k^2 * (g1.Ψ + g1.Φ) ~ Π # anisotropic stress
    ], η; kwargs...)
end

function perturbations_ΛCDM(th::ThermodynamicsSystem, lmax::Int; kwargs...)
    bg = th.bg
    @named g1 = Symboltz.perturbations_metric()
    @named ph = Symboltz.perturbations_photon_hierarchy(bg.sys.g, g1, lmax, true)
    @named cdm = Symboltz.perturbations_matter(bg.sys.g, g1; uinteract=false)
    @named bar = Symboltz.perturbations_matter(bg.sys.g, g1; uinteract=true)
    @named gravpt = Symboltz.perturbations_gravity(bg.sys.g, g1)
    return Symboltz.PerturbationsSystem(bg, th, g1, gravpt, ph, cdm, bar; kwargs...)
end

# TODO: take list of species, each of which "exposes" contributions to δρ and Π
# TODO: support nicer and more general spline input interface
function PerturbationsSystem(bg::BackgroundSystem, th::ThermodynamicsSystem, g::ODESystem, grav::ODESystem, ph::ODESystem, cdm::ODESystem, bar::ODESystem; kwargs...)
    @parameters τspline::CubicSpline dτspline::CubicSpline ddτspline::CubicSpline csb²spline::CubicSpline # TODO: get rid of
    @variables ρc(η) ρb(η) δρr(η) δρc(η) δρb(η) R(η) csb²(η) τ(η) dτ(η) ddτ(η) v(η) dv(η) Δm(η) SSW(η) SD(η) SISW(η) S(η) # TODO: get rid of
    # TODO: do various IC types (adiabatic, isocurvature, ...) from here?
    connections = ODESystem([
        # gravity density and shear stress
        ρb ~ bg.sys.bar.ρ
        ρc ~ bg.sys.cdm.ρ
        δρr ~ ph.δ * bg.sys.rad.ρ
        δρc ~ cdm.δ * ρc
        δρb ~ bar.δ * ρb
        grav.δρ ~ δρr + δρc + δρb # total energy density perturbation
        grav.Π ~ -32π*bg.sys.g.a^2 * bg.sys.rad.ρ*ph.Θ[2] # TODO: add neutrinos
    
        # baryon-photon interactions: Compton (Thomson) scattering # TODO: define connector type?
        R ~ 3/4 * ρb / bg.sys.rad.ρ # Dodelson (5.74)
        bar.uinteraction ~ #=g.k*csb²*bar.δ +=# dτ/R * (bar.u - 3*ph.Θ[1]) # TODO: enable csb² when it seems stable... # TODO: define some common interaction type, e.g. momentum transfer
        ph.ub ~ bar.u
        ph.dτ ~ dτ

        # gauge-independent matter overdensity for matter power spectrum
        Δm ~ (ρc * cdm.Δ + ρb * bar.Δ) / (bg.sys.cdm.ρ + bg.sys.bar.ρ)

        # source functions for CMB power spectrum
        τ ~ spleval(η, τspline) # optical depth
        dτ ~ -exp(spleval(log(η), dτspline))
        ddτ ~ spleval(η, ddτspline)
        csb² ~ spleval(log(η), csb²spline)
        v ~ -dτ * exp(-τ) # visibility function # TODO: take derivatives without additional splines
        dv ~ (dτ^2 - ddτ) * exp(-τ)
        SSW ~ v * (ph.Θ0 + g.Ψ + ph.Π/4) # Sachs-Wolfe
        SD ~ (dv * bar.u + v * Dη(bar.u)) / g.k # Doppler
        SISW ~ exp(-τ) * (Dη(g.Ψ) - Dη(g.Φ)) # integrated Sachs-Wolfe
        S ~ SSW + SD + SISW # CMB source function TODO: add polarization
    ], η; kwargs...)
    sys = compose(connections, g, grav, ph, bar, cdm, bg.sys) # TODO: add background stuff?
    ssys = structural_simplify(sys)
    return PerturbationsSystem(sys, ssys, bg, th)
end

# TODO: get rid of allocations!!! use SVector somehow? see https://docs.sciml.ai/ModelingToolkit/dev/basics/FAQ/#Change-the-unknown-variable-vector-type, also follow 
function solve(pt::PerturbationsSystem, ks::AbstractArray, Ωr0, Ωc0, Ωb0, h, Yp; solver = KenCarp4(), aini = 1e-8, reltol = 1e-10, verbose = false, kwargs...)
    th = pt.th
    bg = th.bg
    th_sol = solve(th, Ωr0, Ωc0, Ωb0, h, Yp) # update spline for dτ (e.g. to propagate derivative information through recombination, if called with dual numbers)
    ΩΛ0 = th_sol.ps[bg.sys.de.Ω0]
    ηs = th_sol[η]
    ηini, ηtoday = ηs[begin], ηs[end]
    ηs = exp.(range(log(ηini), log(ηtoday), length=1024)) # TODO: select determine points adaptively from th_sol # TODO: use saveat in th sol # TODO: CMB spectrum is sensitive to number of points here!
    a, ℰ = th_sol(ηini; idxs=[bg.sys.g.a, bg.sys.g.ℰ])
    τs = th_sol(ηs, idxs=th.sys.τ).u .- th_sol(ηtoday, idxs=th.sys.τ)
    τspline = CubicSpline(τs, ηs; extrapolate=true)
    dτs = DataInterpolations.derivative.(Ref(τspline), ηs)
    dτspline = CubicSpline(log.(.-dτs), log.(ηs); extrapolate=true) # spline this logarithmically for accurayc during integration # TODO: extrapolate valid?
    ddτs = DataInterpolations.derivative.(Ref(τspline), ηs, 2) # spline normally (just observed anyway)
    ddτspline = CubicSpline(ddτs, ηs; extrapolate=true)
    csb²s = th_sol(ηs, idxs=th.sys.cs²).u
    csb²spline = CubicSpline(csb²s, log.(ηs); extrapolate=true)

    p = [bg.sys.rad.Ω0 => Ωr0, bg.sys.cdm.Ω0 => Ωc0, bg.sys.bar.Ω0 => Ωb0, bg.sys.de.Ω0 => ΩΛ0, pt.ssys.g1.k => 0.0, bg.sys.g.H0 => NaN, pt.ssys.τspline => τspline, pt.ssys.dτspline => dτspline, pt.ssys.ddτspline => ddτspline, pt.ssys.csb²spline => csb²spline] # TODO: copy/merge background parameters
    u0 = [bg.sys.g.a => a, bg.sys.g.ℰ => ℰ, pt.ssys.dτ => dτs[begin]]
    prob = ODEProblem(pt.ssys, u0, (ηini, ηtoday), p; jac = true, sparse = false) # TODO: sparse fails with dual numbers # TODO: move into PerturbationsSystem again?
    probs = EnsembleProblem(; safetycopy = false, prob, prob_func = (prob, i, _) -> begin
        k = ks[i]
        verbose && println("$i/$(length(ks)) k = $(k*k0) Mpc/h")
        return remake(prob; p = [pt.ssys.g1.k => k]) # TODO: need use_defaults?
    end)

    return solve(probs, solver, EnsembleThreads(), trajectories = length(ks); reltol, kwargs...) # KenCarp4 and Kvaerno5 seem to work well # TODO: test GPU parallellization
end

# proxy function for evaluating a spline
@register_symbolic spleval(x, spline::CubicSpline)
spleval(x, spline) = spline(x)