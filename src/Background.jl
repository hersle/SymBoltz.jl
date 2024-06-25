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

function background_ΛCDM(; kwargs...)
    @named g = background_metric()
    @named grav = background_gravity_GR(g)
    @named ph = background_radiation(g)
    @named neu = background_radiation(g)
    @named cdm = background_matter(g)
    @named bar = background_matter(g)
    @named de = background_cosmological_constant(g)
    species = [ph, neu, cdm, bar, de]
    initialization_eqs = [
        g.a ~ √(ph.Ω0 + neu.Ω0) * t # analytical radiation-dominated solution
    ]
    defaults = [
        species[end].Ω0 => 1 - sum(s.Ω0 for s in species[begin:end-1])
    ]
    connections = ODESystem([
        grav.ρ ~ sum(s.ρ for s in species);
    ], t; initialization_eqs, defaults, kwargs...)
    return compose(connections, [g; grav; species]...)
end