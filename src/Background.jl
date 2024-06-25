function background_metric(; kwargs...)
    a, ℰ, E = GlobalScope.(@variables a(η) ℰ(η) E(η)) # TODO: more natural way to connect them?
    return ODESystem([
        ℰ ~ Dη(a) / a # ℰ = ℋ/ℋ0
        E ~ ℰ / a # E = H/H0
    ], η; kwargs...)
end

function background_gravity_GR(g; kwargs...)
    @variables ρ(η) ρcrit(η)
    return ODESystem([
        Dη(g.a) ~ √(8π/3 * ρ) * g.a^2 # Friedmann equation
        ρcrit ~ 3/8π * g.E^2 # critical density (H² = 8πG/3 * ρcrit)
    ], η; kwargs...)
end

function background_species(g, w; kwargs...)
    @parameters Ω0 ρ0 = 3/8π*Ω0
    @variables ρ(η) P(η)
    return ODESystem([
        P ~ w * ρ # equation of state
        ρ ~ ρ0 * g.a^(-3*(1+w)) # alternative derivative: Dη(ρ) ~ -3 * g.ℰ * (ρ + P) 
    ], η; kwargs...)
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
        g.a ~ √(ph.Ω0 + neu.Ω0) * η # analytical radiation-dominated solution
    ]
    defaults = [
        species[end].Ω0 => 1 - sum(s.Ω0 for s in species[begin:end-1])
    ]
    connections = ODESystem([
        grav.ρ ~ sum(s.ρ for s in species);
    ], η; initialization_eqs, defaults, kwargs...)
    return compose(connections, [g; grav; species]...)
end