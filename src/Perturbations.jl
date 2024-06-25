struct PerturbationsSystem
    sys::ODESystem
    ssys::ODESystem
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
    bg = th.sys.bg
    @named g1 = perturbations_metric()
    @named ph = perturbations_photon_hierarchy(bg.g, g1, lmax, true)
    @named neu = perturbations_massless_neutrino_hierarchy(bg.g, g1, bg.neu, bg.ph, lmax)
    @named cdm = perturbations_matter(bg.g, g1; uinteract=false)
    @named bar = perturbations_matter(bg.g, g1; uinteract=true)
    @named gravpt = perturbations_gravity(bg.g, g1)
    fν = bg.neu.Ω0 / (bg.ph.Ω0 + bg.neu.Ω0) # TODO: make proper parameter
    defaults = [
        g1.Ψ => -1 / (3/2 + 2*fν/5), # Φ found from solving initialization system
        #g1.Φ => (1 + 2/5*fν) / (3/2 + 2*fν/5), # Ψ found from solving initialization system
    ]
    return PerturbationsSystem(th, g1, gravpt, ph, neu, cdm, bar; defaults, guesses = [g1.Ψ => 1.0; collect(ph.Θ .=> 0.0); collect(neu.Θ .=> 0.0)], kwargs...)
end

# TODO: take list of species, each of which "exposes" contributions to δρ and Π
# TODO: support nicer and more general spline input interface
function PerturbationsSystem(th::ThermodynamicsSystem, g::ODESystem, grav::ODESystem, ph::ODESystem, neu::ODESystem, cdm::ODESystem, bar::ODESystem; kwargs...)
    @variables δργ(η) δρν(η) δρc(η) δρb(η) R(η) Δm(η) # TODO: get rid of
    # TODO: do various IC types (adiabatic, isocurvature, ...) from here?
    bg = th.sys.bg
    connections = ODESystem([
        # gravity density and shear stress
        δργ ~ ph.δ * bg.ph.ρ
        δρν ~ neu.δ * bg.neu.ρ
        δρc ~ cdm.δ * bg.cdm.ρ
        δρb ~ bar.δ * bg.bar.ρ
        grav.δρ ~ δργ + δρν + δρc + δρb # total energy density perturbation
        grav.Π ~ -32π*bg.g.a^2 * (bg.ph.ρ*ph.Θ[2] + bg.neu.ρ*neu.Θ[2])
    
        # baryon-photon interactions: Compton (Thomson) scattering # TODO: define connector type?
        R ~ 3/4 * bg.bar.ρ / bg.ph.ρ # Dodelson (5.74)
        bar.uinteraction ~ #=g.k*csb²*bar.δ +=# th.sys.dτ/R * (bar.u - 3*ph.Θ[1]) # TODO: enable csb² when it seems stable... # TODO: define some common interaction type, e.g. momentum transfer
        ph.ub ~ bar.u
        ph.dτ ~ th.sys.dτ

        # gauge-independent matter overdensity for matter power spectrum
        Δm ~ (bg.cdm.ρ * cdm.Δ + bg.bar.ρ * bar.Δ) / (bg.cdm.ρ + bg.bar.ρ)
    ], η; kwargs...)
    sys = compose(connections, g, grav, ph, neu, bar, cdm, th.sys)
    ssys = structural_simplify(sys)
    return PerturbationsSystem(sys, ssys, th)
end

# TODO: get rid of allocations!!! use SVector somehow? see https://docs.sciml.ai/ModelingToolkit/dev/basics/FAQ/#Change-the-unknown-variable-vector-type, also follow 
function solve(pt::PerturbationsSystem, ks::AbstractArray, Ωγ0, Ων0, Ωc0, Ωb0, h, Yp; solver = Rodas5P(), aini = 1e-7, reltol = 1e-7, verbose = false, kwargs...)
    th = pt.th
    bg = th.sys.bg # th.bg
    th_sol = solve(th, Ωγ0, Ων0, Ωc0, Ωb0, h, Yp; aini) # update spline for dτ (e.g. to propagate derivative information through recombination, if called with dual numbers)
    ηs = th_sol[η]
    ηini, ηtoday = ηs[begin], ηs[end]
    a = th_sol(ηini; idxs=bg.g.a)

    u0 = [bg.g.a => a] # TODO: does this cause overinitialization?
    p = [bg.ph.Ω0 => Ωγ0, bg.neu.Ω0 => Ων0, bg.cdm.Ω0 => Ωc0, bg.bar.Ω0 => Ωb0, bg.g.H0 => H100 * h, th.sys.Yp => Yp] # TODO: copy/merge background parameters
    # TODO: improve performance!
    probs = EnsembleProblem(; safetycopy = false, prob = nothing#=prob_uninit=#, prob_func = (prob_uninit, i, _) -> begin
        verbose && println("$i/$(length(ks)) k = $(ks[i]*k0) Mpc/h")
        return ODEProblem(pt.ssys, u0, (ηini, ηtoday), merge(Dict(p), Dict(k => ks[i]))) # works, but even slower
    end)

    return solve(probs, solver, EnsembleThreads(), trajectories = length(ks); reltol, kwargs...) # KenCarp4 and Kvaerno5 seem to work well # TODO: test GPU parallellization
end