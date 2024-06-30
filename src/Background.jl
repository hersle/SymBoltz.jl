function background_metric(; kwargs...)
    a, ℰ, E, H, ℋ = GlobalScope.(@variables a(t) ℰ(t) E(t) H(t) ℋ(t)) # TODO: more natural way to connect them?
    H0, h = GlobalScope.(@parameters H0 h)
    return ODESystem([
        ℰ ~ D(a) / a # ℰ = ℋ/ℋ0
        E ~ ℰ / a # E = H/H0
        ℋ ~ ℰ * H0
        H ~ E * H0
    ], t, [a, ℰ, E, H, ℋ], [H0, h]; defaults = [H0 => H100 * h], kwargs...)
end

function background_gravity_GR(g; kwargs...)
    @variables ρ(t) ρcrit(t)
    return ODESystem([
        D(g.a) ~ √(8π/3 * ρ) * g.a^2 # Friedmann equation
        ρcrit ~ 3/8π * g.E^2 # critical density (H² = 8πG/3 * ρcrit)
    ], t; kwargs...)
end

function background_species(g, w; kwargs...)
    @parameters Ω0 ρ0 = 3/8π*Ω0
    @variables ρ(t) P(t)
    return ODESystem([
        P ~ w * ρ # equation of state
        ρ ~ ρ0 * g.a^(-3*(1+w)) # alternative derivative: D(ρ) ~ -3 * g.ℰ * (ρ + P) 
    ], t; kwargs...)
end
background_radiation(g; kwargs...) = background_species(g, 1//3; kwargs...)
background_matter(g; kwargs...) = background_species(g, 0; kwargs...)
background_cosmological_constant(g; kwargs...) = background_species(g, -1; kwargs...)

function background_photons(g; kwargs...)
    @parameters T0
    @variables T(t)    
    return extend(background_radiation(g; kwargs...), ODESystem([
        T ~ T0 / g.a # alternative derivative: D(Tγ) ~ -1*Tγ * g.ℰ
    ], t; name=:ph))
end

function background_massive_neutrinos(g; nx=20, xmax=50.0, kwargs...)
    pars = @parameters Ω0_massless Ω0 ρ0 m ∑m T0 y0
    vars = @variables ρ(t) T(t) y(t) P(t) w(t)

    # TODO: also allow adaptive quadgk
    x, W = gauss(nx, 0.0, xmax)
    f0(x) = 1 / (exp(x) + 1) # TODO: exp(x) or exp(ϵ)
    ϵ(x, y) = √(x^2 + y^2)
    Iρ(y) = ∫((@. x^2 * ϵ(x, y) * f0(x)), W) # Iρ(0) = 7π^4/120
    IP(y) = ∫((@. x^4 / ϵ(x, y) * f0(x)), W) # IP(0) = Iρ(0)

    eqs = [
        T ~ T0 / g.a
        y ~ m*c^2 / (kB*T)
        ρ ~ ρ0/g.a^4 * Iρ(y) / Iρ(y0) # have ρ = Cρ * Iρ(y) / a⁴, so Cρ = ρ0 * 1⁴ / Iρ(y0)
        P ~ 1/3 * ρ0/g.a^4 * IP(y) / Iρ(y0) # have P = CP * IP(y) / a⁴, and in the early universe Iρ(y→0) → IP(y→0) and P/ρ = CP * IP(y) / (Cρ * Iρ(y)) → CP/Cρ → 1/3, so CP = Cρ/3
        w ~ P / ρ
    ]
    defaults = [
        Ω0 => Ω0_massless * Iρ(y0)/Iρ(0) * y0 # ≈ Ω0_massless * (3ζ(3)/2)/(7π^4/120)*y0 for y0 → ∞
        ρ0 => 3/8π * Ω0
        ∑m => 0.06 * eV/c^2 # total mass
        m => ∑m / 3 # single mass
        y0 => m*c^2 / (kB*T0)
    ]
    return ODESystem(eqs, t, vars, pars; defaults, kwargs...)
end

function background_ΛCDM(; kwargs...)
    @named g = background_metric()
    @named grav = background_gravity_GR(g)
    @named ph = background_photons(g)
    @named neu = background_radiation(g)
    @named mneu = background_massive_neutrinos(g)
    @named cdm = background_matter(g)
    @named bar = background_matter(g)
    @named de = background_cosmological_constant(g)
    species = [ph, neu, mneu, cdm, bar, de]
    initialization_eqs = [g.a ~ (2t)^(1/2) * (ph.Ω0 + neu.Ω0 + mneu.Ω0_massless)^(1/4)] # analytical radiation-dominated solution # TODO: add effect from massive neutrinos
    defaults = [
        species[end].Ω0 => 1 - sum(s.Ω0 for s in species[begin:end-1]) # TODO: solve nonlinear system
        ph.T0 => (ph.ρ0 * 15/π^2 * g.H0^2/G * ħ^3*c^5)^(1/4) / kB # TODO: move to photon system
        neu.Ω0 => (3.046/3) * (4/11)^(4/3) * ph.Ω0 # TODO: make Neff parameter
        mneu.T0 => (3.046/3)^(1/4) * (4/11)^(1/3) * ph.T0 # same as for massless neutrinos # TODO: make Neff parameter etc.
        mneu.Ω0_massless => neu.Ω0 # Ω0 for corresponding massless neutrinos
    ]
    eqs = [grav.ρ ~ sum(s.ρ for s in species)]
    comps = [g; grav; species]
    bg = ODESystem(eqs, t; initialization_eqs, defaults, kwargs...)
    return compose(bg, comps...)
end