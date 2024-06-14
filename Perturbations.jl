struct PerturbationsSystem
    sys::ODESystem
    ssys::ODESystem
    bg::BackgroundSystem
    th::ThermodynamicsSystem
end

# TODO: merge with background metric!
function perturbations_metric(; kwargs...)
    Φ, Ψ = GlobalScope.(@variables Φ(η) Ψ(η))
    return ODESystem(Equation[], η, [Φ, Ψ], [k]; kwargs...)
end

function perturbations_species(g0, g1, w, cs² = w, w′ = 0, σ = 0; uinteract = false, kwargs...)
    @assert w′ == 0 && σ == 0 # TODO: relax (need to include in ICs)
    @variables δ(η) u(η) Δ(η) uinteraction(η)
    eqs = [
        Dη(δ) ~ -(1+w)*(k*u+3*Dη(g1.Φ)) - 3*g0.ℰ*(cs²-w)*δ # Bertschinger & Ma (30) with Φ -> -Φ; or Baumann (4.4.173) with Φ -> -Φ
        Dη(u) ~ -g0.ℰ*(1-3*w)*u - w′/(1+w)*u + cs²/(1+w)*k*δ + k*g1.Ψ - k*σ + uinteraction # Bertschinger & Ma (30) with θ = kv
        Δ ~ δ + 3 * (1+w) * g0.ℰ/k * u # Baumann (4.2.144) with v -> -u
    ]
    initialization_eqs = [
        δ ~ -3/2 * (1+w) * g1.Ψ # adiabatic: δᵢ/(1+wᵢ) == δⱼ/(1+wⱼ) (https://cmb.wintherscoming.no/theory_initial.php#adiabatic)
        u ~ -1/2 * k*η * g1.Ψ # TODO: include σ ≠ 0 # solve u′ + ℋ(1-3w)u = w/(1+w)*kδ + kΨ with Ψ=const, IC for δ, Φ=-Ψ, ℋ=H₀√(Ωᵣ₀)/a after converting ′ -> d/da by gathering terms with u′ and u in one derivative using the trick to multiply by exp(X(a)) such that X′(a) will "match" the terms in front of u
    ]
    !uinteract && push!(eqs, uinteraction ~ 0)
    return ODESystem(eqs, η; initialization_eqs, kwargs...)
end

perturbations_matter(g0, g1; kwargs...) = perturbations_species(g0, g1, 0; kwargs...)
perturbations_radiation(g0, g1; kwargs...) = perturbations_species(g0, g1, 1//3; kwargs...)
perturbations_cosmological_constant(g0, g1; kwargs...) = perturbations_species(g0, g1, -1; kwargs...) # TODO: ill-defined?

function perturbations_photon_hierarchy(g0, g1, lmax=6, polarization=true; kwargs...)
    @variables Θ0(η) Θ(η)[1:lmax] δ(η) dτ(η) ub(η) Π(η) ΘP0(η) ΘP(η)[1:lmax]
    eqs = [
        Dη(Θ0) + k*Θ[1] ~ -Dη(g1.Φ)
        Dη(Θ[1]) - k/3*(Θ0-2*Θ[2]) ~ k/3*g1.Ψ - dτ/3 * (ub - 3*Θ[1])
        [Dη(Θ[l]) ~ k/(2*l+1) * (l*Θ[l-1] - (l+1)*Θ[l+1]) + dτ * (Θ[l] - Π/10*δkron(l,2)) for l in 2:lmax-1]... # TODO: Π in last term here?
        Dη(Θ[lmax]) ~ k*Θ[lmax-1] - (lmax+1) * Θ[lmax] / η + dτ * (Θ[lmax] - Π/10*δkron(lmax,2))
        δ ~ 4*Θ0
        Π ~ Θ[2] + ΘP[2] + ΘP0
        (polarization ? [
            Dη(ΘP0) + k*ΘP[1] ~ dτ * (ΘP0 - Π/2)
            Dη(ΘP[1]) - k/(2*1+1) * (1*ΘP0 - (1+1)*ΘP[1+1]) ~ dτ * (ΘP[1] - Π/10*δkron(1,2))
            [Dη(ΘP[l]) - k/(2*l+1) * (l*ΘP[l-1] - (l+1)*ΘP[l+1]) ~ dτ * (ΘP[l] - Π/10*δkron(l,2)) for l in 2:lmax-1]...
            Dη(ΘP[lmax]) ~ k*ΘP[lmax-1] - (lmax+1) * ΘP[lmax] / η + dτ * ΘP[lmax]
        ] : [
            ΘP0 ~ 0, collect(ΘP .~ 0)... # pin to zero
        ])...
    ]
    initialization_eqs = [
        Θ0 ~ -1/2 * g1.Ψ, # Dodelson (7.89)
        Θ[1] ~ 1/6 * k*η * g1.Ψ, # Dodelson (7.95)
        Θ[2] ~ (polarization ? -8/15 : -20/45) * k/dτ * Θ[1], # depends on whether polarization is included # TODO: move to initialization_eqs?
        [Θ[l] ~ -l/(2*l+1) * k/dτ * Θ[l-1] for l in 3:lmax]...,
        ΘP0 ~ 5/4 * Θ[2],
        ΘP[1] ~ -1/4 * k/dτ * Θ[2],
        ΘP[2] ~ 1/4 * Θ[2],
        [ΘP[l] ~ -l/(2*l+1) * k/dτ * ΘP[l-1] for l in 3:lmax]...
    ]
    return ODESystem(eqs, η; initialization_eqs, kwargs...)
end

function perturbations_massless_neutrino_hierarchy(g0, g1, neu0, ph0, lmax=6; kwargs...)
    @variables Θ0(η) Θ(η)[1:lmax] δ(η)
    eqs = [
        Dη(Θ0) ~ -k*Θ[1] - Dη(g1.Φ)
        Dη(Θ[1]) ~ k/3 * (Θ0 - 2*Θ[2] + g1.Ψ)
        [Dη(Θ[l]) ~ k/(2*l+1) * (l*Θ[l-1] - (l+1)*Θ[l+1]) for l in 2:lmax-1]...
        Dη(Θ[lmax]) ~ k*Θ[lmax-1] - (lmax+1)/η*Θ[lmax]
        δ ~ 4*Θ0
    ]
    initialization_eqs = [
        Θ0 ~ -1/2 * g1.Ψ,
        Θ[1] ~ 1/6 * k*η * g1.Ψ,
        Θ[2] ~ 1/30 * (k*η)^2 * g1.Ψ, # Dodelson (7.122) and (7.123), # (k*g0.a)^2 / (80π*ρr0) * g1.Ψ, # 2/15 * (k*η)^2 * g1.Ψ, # -k^2*g0.a^2 / (32π * (15/4*ρr0 + ρν0)), # TODO: how to set ICs consistently with Ψ, Π and Θν2?
        [Θ[l] ~ 1/(2*l+1) * k*η * Θ[l-1] for l in 3:lmax]...
    ]
    return ODESystem(eqs, η; initialization_eqs, kwargs...)
end

function perturbations_gravity(g0, g1; kwargs...)
    @variables δρ(η) Π(η)
    return ODESystem([
        Dη(g1.Φ) ~ 4π*g0.a^2*δρ/(3*g0.ℰ) - k^2*g1.Φ/(3*g0.ℰ) + g0.ℰ*g1.Ψ
        k^2 * (g1.Ψ + g1.Φ) ~ Π # anisotropic stress
    ], η; kwargs...)
end

function perturbations_ΛCDM(th::ThermodynamicsSystem, lmax::Int; kwargs...)
    bg = th.bg
    @named g1 = perturbations_metric()
    @named ph = perturbations_photon_hierarchy(bg.sys.g, g1, lmax, true)
    @named neu = perturbations_massless_neutrino_hierarchy(bg.sys.g, g1, bg.sys.neu, bg.sys.ph, lmax)
    @named cdm = perturbations_matter(bg.sys.g, g1; uinteract=false)
    @named bar = perturbations_matter(bg.sys.g, g1; uinteract=true)
    @named gravpt = perturbations_gravity(bg.sys.g, g1)
    fν = bg.sys.neu.Ω0 / (bg.sys.ph.Ω0 + bg.sys.neu.Ω0) # TODO: make proper parameter
    defaults = [
        g1.Ψ => -1 / (3/2 + 2*fν/5), # Φ found from solving initialization system
        #g1.Φ => (1 + 2/5*fν) / (3/2 + 2*fν/5), # Ψ found from solving initialization system
    ]
    return PerturbationsSystem(bg, th, g1, gravpt, ph, neu, cdm, bar; defaults, guesses = [g1.Ψ => 1.0, neu.Θ[1] => 1e-5, ph.Θ[6] => 0.0], kwargs...)
end

# TODO: take list of species, each of which "exposes" contributions to δρ and Π
# TODO: support nicer and more general spline input interface
function PerturbationsSystem(bg::BackgroundSystem, th::ThermodynamicsSystem, g::ODESystem, grav::ODESystem, ph::ODESystem, neu::ODESystem, cdm::ODESystem, bar::ODESystem; kwargs...)
    @parameters τspline::CubicSpline dτspline::CubicSpline ddτspline::CubicSpline csb²spline::CubicSpline # TODO: get rid of
    @variables δργ(η) δρν(η) δρc(η) δρb(η) R(η) csb²(η) τ(η) dτ(η) ddτ(η) v(η) dv(η) Δm(η) SSW(η) SD(η) SISW(η) S(η) # TODO: get rid of
    # TODO: do various IC types (adiabatic, isocurvature, ...) from here?
    connections = ODESystem([
        # gravity density and shear stress
        δργ ~ ph.δ * bg.sys.ph.ρ
        δρν ~ neu.δ * bg.sys.neu.ρ
        δρc ~ cdm.δ * bg.sys.cdm.ρ
        δρb ~ bar.δ * bg.sys.bar.ρ
        grav.δρ ~ δργ + δρν + δρc + δρb # total energy density perturbation
        grav.Π ~ -32π*bg.sys.g.a^2 * (bg.sys.ph.ρ*ph.Θ[2] + bg.sys.neu.ρ*neu.Θ[2])
    
        # baryon-photon interactions: Compton (Thomson) scattering # TODO: define connector type?
        R ~ 3/4 * bg.sys.bar.ρ / bg.sys.ph.ρ # Dodelson (5.74)
        bar.uinteraction ~ #=g.k*csb²*bar.δ +=# dτ/R * (bar.u - 3*ph.Θ[1]) # TODO: enable csb² when it seems stable... # TODO: define some common interaction type, e.g. momentum transfer
        ph.ub ~ bar.u
        ph.dτ ~ dτ

        # gauge-independent matter overdensity for matter power spectrum
        Δm ~ (bg.sys.cdm.ρ * cdm.Δ + bg.sys.bar.ρ * bar.Δ) / (bg.sys.cdm.ρ + bg.sys.bar.ρ)

        # source functions for CMB power spectrum
        τ ~ spleval(η, τspline) # optical depth
        dτ ~ -exp(spleval(log(η), dτspline))
        ddτ ~ spleval(η, ddτspline)
        csb² ~ spleval(log(η), csb²spline)
        v ~ -dτ * exp(-τ) # visibility function # TODO: take derivatives without additional splines
        dv ~ (dτ^2 - ddτ) * exp(-τ)

        # not working after added initialization system TODO: restore
        SSW ~ v * (ph.Θ0 + g.Ψ + ph.Π/4) # Sachs-Wolfe
        SD ~ (dv * bar.u + v * Dη(bar.u)) / g.k # Doppler
        SISW ~ exp(-τ) * (Dη(g.Ψ) - Dη(g.Φ)) # integrated Sachs-Wolfe
        S ~ SSW + SD + SISW # CMB source function TODO: add polarization
    ], η; kwargs...)
    sys = compose(connections, g, grav, ph, neu, bar, cdm, bg.sys) # TODO: add background stuff?
    ssys = structural_simplify(sys)
    return PerturbationsSystem(sys, ssys, bg, th)
end

# TODO: get rid of allocations!!! use SVector somehow? see https://docs.sciml.ai/ModelingToolkit/dev/basics/FAQ/#Change-the-unknown-variable-vector-type, also follow 
function solve(pt::PerturbationsSystem, ks::AbstractArray, Ωγ0, Ων0, Ωc0, Ωb0, h, Yp; solver = KenCarp4(), aini = 1e-8, reltol = 1e-10, verbose = false, kwargs...)
    th = pt.th
    bg = th.bg
    th_sol = solve(th, Ωγ0, Ων0, Ωc0, Ωb0, h, Yp) # update spline for dτ (e.g. to propagate derivative information through recombination, if called with dual numbers)
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

    u0 = [bg.sys.g.a => a] # TODO: does this cause overinitialization?
    p = [bg.sys.ph.Ω0 => Ωγ0, bg.sys.neu.Ω0 => Ων0, bg.sys.cdm.Ω0 => Ωc0, bg.sys.bar.Ω0 => Ωb0, bg.sys.de.Ω0 => ΩΛ0, k => 1.0, bg.sys.g.H0 => NaN, pt.ssys.τspline => τspline, pt.ssys.dτspline => dτspline, pt.ssys.ddτspline => ddτspline, pt.ssys.csb²spline => csb²spline] # TODO: copy/merge background parameters
    prob_uninit = ODEProblem(pt.ssys, u0, (ηini, ηtoday), p; jac = true) # TODO: sparse fails with dual numbers # TODO: cache in PerturbationsSystem again?
    iprob_uninit = ModelingToolkit.InitializationProblem(pt.ssys, ηini, u0, p; warn_initialize_determined = false) # TODO: cache in PerturbationsSystem again?
    
    vars = unknowns(pt.ssys)
    # TODO: improve performance!
    probs = EnsembleProblem(; safetycopy = false, prob = prob_uninit, prob_func = (prob_uninit, i, _) -> begin
        verbose && println("$i/$(length(ks)) k = $(ks[i]*k0) Mpc/h")
        iprob = remake(iprob_uninit, p = [k => ks[i]])
        isol = solve(iprob; verbose) # TODO: reduce time spent here!
        return remake(prob_uninit, u0 = isol[vars], p = [k => ks[i]]) # TODO: avoid vars indexing?
        #return ODEProblem(pt.ssys, u0, (ηini, ηtoday), merge(Dict(p), Dict(k => ks[i]))) # works, but even slower
    end)

    return solve(probs, solver, EnsembleThreads(), trajectories = length(ks); reltol, kwargs...) # KenCarp4 and Kvaerno5 seem to work well # TODO: test GPU parallellization
end

# proxy function for evaluating a spline
@register_symbolic spleval(x, spline::CubicSpline)
spleval(x, spline) = spline(x)