# TODO: merge with background metric!
function perturbations_metric(; kwargs...)
    Φ, Ψ = GlobalScope.(@variables Φ(t) Ψ(t))
    return ODESystem(Equation[], t, [Φ, Ψ], [k]; kwargs...)
end

function perturbations_species(g0, g1, w, cs² = w, ẇ = 0, σ = 0; θinteract = false, kwargs...)
    @assert ẇ == 0 && σ == 0 # TODO: relax (need to include in ICs)
    @variables δ(t) θ(t) Δ(t) θinteraction(t)
    eqs = [
        D(δ) ~ -(1+w)*(θ-3*D(g1.Φ)) - 3*g0.ℰ*(cs²-w)*δ # Bertschinger & Ma (30) with Φ -> -Φ; or Baumann (4.4.173) with Φ -> -Φ
        D(θ) ~ -g0.ℰ*(1-3*w)*θ - ẇ/(1+w)*θ + cs²/(1+w)*k^2*δ - k^2*σ + k^2*g1.Ψ + θinteraction # Bertschinger & Ma (30) with θ = kv
        # TODO: Δ ~ δ + 3(1+w) * g0.ℰ/θ # Baumann (4.2.144) with v -> -u
    ]
    initialization_eqs = [
        δ ~ -3/2 * (1+w) * g1.Ψ # adiabatic: δᵢ/(1+wᵢ) == δⱼ/(1+wⱼ) (https://cmb.wintherscoming.no/theory_initial.php#adiabatic)
        θ ~ 1/2 * (k^2*t) * g1.Ψ # # TODO: fix with θ: -1/2 * k*t * g1.Ψ # TODO: include σ ≠ 0 # solve u′ + ℋ(1-3w)u = w/(1+w)*kδ + kΨ with Ψ=const, IC for δ, Φ=-Ψ, ℋ=H₀√(Ωᵣ₀)/a after converting ′ -> d/da by gathering terms with u′ and u in one derivative using the trick to multiply by exp(X(a)) such that X′(a) will "match" the terms in front of u
    ]
    !θinteract && push!(eqs, θinteraction ~ 0)
    return ODESystem(eqs, t; initialization_eqs, kwargs...)
end

perturbations_matter(g0, g1; kwargs...) = perturbations_species(g0, g1, 0; kwargs...)
perturbations_radiation(g0, g1; kwargs...) = perturbations_species(g0, g1, 1//3; kwargs...)
perturbations_cosmological_constant(g0, g1; kwargs...) = perturbations_species(g0, g1, -1; kwargs...) # TODO: ill-defined?

function perturbations_photon_hierarchy(g0, g1, lmax=6, polarization=true; kwargs...) # TODO: enable polarization
    @variables F0(t) F(t)[1:lmax] δ(t) θ(t) σ(t) τ̇(t) θb(t) Π(t) G0(t) G(t)[1:lmax]
    eqs = [
        D(F0) ~ -k*F[1] + 4*D(g1.Φ)
        D(F[1]) ~ k/3*(F0-2*F[2]+4*g1.Ψ) - 4/3 * τ̇/k * (θb - θ) # TODO: ±τ̇ here and elsewhere ? # TODO: factor k in last term?
        D(F[2]) ~ 2/5*k*F[1] - 3/5*k*F[3] + 9/10*τ̇*F[2] - 1/10*τ̇*(G0+G[2])
        [D(F[l]) ~ k/(2*l+1) * (l*F[l-1] - (l+1)*F[l+1]) + τ̇*F[l] for l in 3:lmax-1]... # TODO: Π in last term here?
        D(F[lmax]) ~ k*F[lmax-1] - (lmax+1) / t * F[lmax] + τ̇ * F[lmax] # TODO: assumes lmax ≥ ???
        δ ~ F0
        θ ~ 3/4*k*F[1]
        σ ~ F[2]/2
        Π ~ F[2] + G0 + G[2]
        (polarization ? [
            D(G0) ~ k * (-G[1]) - τ̇ * (-G0 + Π/2)
            D(G[1]) ~ k/3 * (G0 - 2*G[2]) - τ̇ * (-G[1])
            [D(G[l]) ~ k/(2*l+1) * (l*G[l-1] - (l+1)*G[l+1]) - τ̇ * (-G[l] + Π/10*δkron(l,2)) for l in 2:lmax-1]... # TODO: collect all equations here once G[0] works
            D(G[lmax]) ~ k*G[lmax-1] - (lmax+1) / t * G[lmax] + τ̇ * G[lmax]
        ] : [
            G0 ~ 0, collect(G .~ 0)... # pin to zero
        ])...
    ]
    initialization_eqs = [
        δ ~ -2 * g1.Ψ # Dodelson (7.89)
        θ ~ 1/2 * (k^2*t) * g1.Ψ # Dodelson (7.95)
        F[2] ~ 0 # (polarization ? -8/15 : -20/45) * k/dτ * Θ[1], # depends on whether polarization is included # TODO: move to initialization_eqs?
        [F[l] ~ 0 #=-l/(2*l+1) * k/dτ * Θ[l-1]=# for l in 3:lmax]...
        G0 ~ 0 #5/4 * Θ[2],
        G[1] ~ 0 #-1/4 * k/dτ * Θ[2],
        G[2] ~ 0 #1/4 * Θ[2],
        [G[l] ~ 0 #=-l/(2*l+1) * k/dτ * ΘP[l-1]=# for l in 3:lmax]...
    ]
    return ODESystem(eqs, t; initialization_eqs, kwargs...)
end

function perturbations_massless_neutrino_hierarchy(g0, g1, neu0, ph0, lmax=6; kwargs...)
    @variables F0(t) F(t)[1:lmax] δ(t) θ(t) σ(t)
    eqs = [
        D(F0) ~ -k*F[1] + 4*D(g1.Φ)
        D(F[1]) ~ k/3*(F0-2*F[2]+4*g1.Ψ)
        [D(F[l]) ~ k/(2*l+1) * (l*F[l-1] - (l+1)*F[l+1]) for l in 2:lmax-1]...
        D(F[lmax]) ~ k*F[lmax-1] - (lmax+1) / t * F[lmax] # this is the same cutoff as for photons; gives best agreement with CLASS
        #F[lmax] ~ (2*lmax-1) / (k*t) * F[lmax-1] - F[lmax-2] # Bertschinger & Ma say they use this cutoff, which is different from the photon cutoff; it gives results somewhat different from CLASS
        δ ~ F0
        θ ~ 3/4*k*F[1]
        σ ~ F[2]/2
    ]
    initialization_eqs = [
        δ ~ -2 * g1.Ψ # adiabatic: δᵢ/(1+wᵢ) == δⱼ/(1+wⱼ) (https://cmb.wintherscoming.no/theory_initial.php#adiabatic)
        θ ~ 1/2 * (k^2*t) * g1.Ψ
        σ ~ 1/15 * (k*t)^2 * g1.Ψ # TODO: how to set ICs consistently with Ψ, Π and Θν2?
        [F[l] ~ 0 #=1/(2*l+1) * k*t * Θ[l-1]=# for l in 3:lmax]...
    ]
    return ODESystem(eqs, t; initialization_eqs, kwargs...)
end

# TODO: don't duplicate things in background neutrinos
# TODO: use vector equations?
# TODO: check with MB
function perturbations_massive_neutrino_hierarchy(g0, g1; nx=5, lmax=4, kwargs...)
    x, W = gauss(f0, nx, 0.0, 1000.0) # reduced momentum bins x = q*c / (kB*T0) # these points give accurate integral for Iρmν in the background, at least # TODO: ok for perturbations? # TODO: also include common x^2 factor in weighting?
    vars = @variables T(t) y(t) ρ(t) δρ(t) δ(t) P(t) σ(t) u(t) ψ0(t)[1:nx] ψ(t)[1:nx, 1:lmax] dlnf0_dlnx(t)[1:nx] ϵ(t)[1:nx]
    eqs = [
        ρ ~ 1/g0.a^4 * ∫(collect(@. x^2 * ϵ), W) # analytical solution with initial y≈0 is 7/120*π^4 / a^4 # TODO: don't duplicate background
        δρ ~ 1/g0.a^4 * ∫(collect(@. x^2 * ϵ * ψ0), W) # analytical solution with initial y≈0 and ψ0 below is -7/60*π^4 * Ψ/a^4
        δ ~ δρ / ρ

        P ~ 1/3 / g0.a^4 * ∫(collect(@. x^4 / ϵ), W)
        σ ~ 8π/3 / g0.a^4 * ∫(collect(@. x^4 / ϵ * ψ[:,2]), W) / (ρ + P)

        u ~ 4π / g0.a^4 * ∫(collect(@. x^3 * ψ[:, 1]), W) / (ρ + P)
    ]
    initialization_eqs = []
    for i in 1:nx
        push!(eqs, [ # TODO: write shorter with vector equations and collect
            ϵ[i] ~ √(x[i]^2 + y^2) # TODO: use z to match reduced x and y?
            dlnf0_dlnx[i] ~ -x[i] / (1 + exp(-x[i]))
            D(ψ0[i]) ~ -k * x[i]/ϵ[i] * ψ[i,1] - D(g1.Φ) * dlnf0_dlnx[i]
            D(ψ[i,1]) ~ +k/3 * x[i]/ϵ[i] * (ψ0[i] - 2ψ[i,2]) - k/3 * ϵ[i]/x[i] * g1.Ψ * dlnf0_dlnx[i]
            [D(ψ[i,l]) ~ k/(2l+1) * x[i]/ϵ[i] * (l*ψ[i,l-1] - (l+1)*ψ[i,l+1]) for l in 2:lmax-1]...
            ψ[i,lmax] ~ (2lmax-1) * ϵ[i]/x[i] * ψ[i,lmax-1] / (k*t) - ψ[i,lmax-2]
        ]...)
        push!(initialization_eqs, [
            ψ0[i] ~ +1/2 * g1.Ψ * dlnf0_dlnx[i]
            ψ[i,1] ~ 0 #-1/6 * k*t * g1.Ψ * ϵ/x * dlnf0_dlnx
            [ψ[i,l] ~ 0 for l in 2:lmax-1] # TODO: proper ICs    
        ]...)
    end
    return ODESystem(eqs, t, vars, []; initialization_eqs, kwargs...)
end

function perturbations_gravity(g0, g1; kwargs...)
    @variables δρ(t) Π(t)
    return ODESystem([
        D(g1.Φ) ~ -4π/3*g0.a^2/g0.ℰ*δρ - k^2/(3*g0.ℰ)*g1.Φ - g0.ℰ*g1.Ψ
        k^2 * (g1.Φ - g1.Ψ) ~ 12π * g0.a^2 * Π # TODO: Π = (ρ+P)*σ # anisotropic stress
    ], t; kwargs...)
end

# TODO: take list of species, each of which "exposes" contributions to δρ and Π
# TODO: support nicer and more general spline input interface
function perturbations_ΛCDM(th::ODESystem, lmax::Int; spline_th=false, kwargs...)
    bg = th.bg
    @named g1 = perturbations_metric()
    @named ph = perturbations_photon_hierarchy(bg.g, g1, lmax, true)
    @named neu = perturbations_massless_neutrino_hierarchy(bg.g, g1, bg.neu, bg.ph, lmax)
    #@named mneu = perturbations_massive_neutrino_hierarchy(bg.g, g1)
    @named cdm = perturbations_matter(bg.g, g1; θinteract=false)
    @named bar = perturbations_matter(bg.g, g1; θinteract=true)
    @named grav = perturbations_gravity(bg.g, g1)

    # TODO: do various IC types (adiabatic, isocurvature, ...) from here?
    pars = convert(Vector{Any}, @parameters fν C)
    vars = @variables δργ(t) δρν(t) δρmν(t) δρc(t) δρb(t) R(t) Δm(t) dτ(t) πν(t)
    defaults = [
        C => +1/2 # TODO: +1/2?
        fν => bg.neu.Ω0 / (bg.neu.Ω0 + bg.ph.Ω0)
        g1.Ψ => 20C / (15 + 4fν) # Φ found from solving initialization system # TODO: is this correct when having both massless and massive neutrinos?
        #g1.Φ => (1 + 2/5*fν) / (3/2 + 2*fν/5) # Ψ found from solving initialization system
    ]
    guesses = [
        g1.Ψ => 1.0;
        collect(ph.F .=> 0.0);
        collect(neu.F .=> 0.0)
    ]
    eqs = [
        # gravity density and shear stress
        grav.δρ ~ ph.δ*bg.ph.ρ + cdm.δ*bg.cdm.ρ + bar.δ*bg.bar.ρ + neu.δ*bg.neu.ρ # + mneu.δ*bg.mneu.ρ # total energy density perturbation
        grav.Π ~ (bg.ph.ρ+bg.ph.P)*ph.σ + (bg.neu.ρ+bg.neu.P)*neu.σ # + (bg.mneu.ρ+bg.mneu.P)*mneu.σ
    
        # baryon-photon interactions: Compton (Thomson) scattering # TODO: define connector type?
        R ~ 4*bg.ph.ρ / (3*bg.bar.ρ)
        bar.θinteraction ~ #=g.k^2*csb²*bar.δ +=# -th.rec.dτ * R * (ph.θ - bar.θ) # TODO: enable csb² when it seems stable... # TODO: define some common interaction type, e.g. momentum transfer
        ph.θb ~ bar.θ
        ph.τ̇ ~ th.rec.dτ

        # gauge-independent matter overdensity for matter power spectrum
        #Δm ~ (bg.cdm.ρ * cdm.Δ + bg.bar.ρ * bar.Δ) / (bg.cdm.ρ + bg.bar.ρ)

        # TODO: combine bg+pt systems
        #mneu.T ~ bg.mneu.T
        #mneu.y ~ bg.mneu.y
    ]
    
    connections = ODESystem(eqs, t, vars, pars; defaults, guesses, kwargs...)
    return compose(connections, g1, grav, ph, neu, #=mneu,=# bar, cdm, th)
end