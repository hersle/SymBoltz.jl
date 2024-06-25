struct BackgroundSystem
    # TODO: add fields like metric, gravity, species, ...?
    sys::ODESystem # unsimplified
    ssys::ODESystem # simplified
    prob::ODEProblem
end

function background_metric(; kwargs...)
    a, ℰ, E, H, ℋ = GlobalScope.(@variables a(η) ℰ(η) E(η) H(η) ℋ(η)) # TODO: more natural way to connect them?
    H0, = GlobalScope.(@parameters H0)
    return ODESystem([
        ℰ ~ Dη(a) / a # ℰ = ℋ/ℋ0
        E ~ ℰ / a # E = H/H0
        ℋ ~ ℰ * H0
        H ~ E * H0
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
    defaults = [species[end].Ω0 => 1 - sum(s.Ω0 for s in species[begin:end-1])]
    return BackgroundSystem(g, grav, species; defaults, kwargs...)
end

function BackgroundSystem(g::ODESystem, grav::ODESystem, species::AbstractArray{ODESystem}; jac=true, kwargs...)
    components = [g; grav; species]
    connections = ODESystem([
        grav.ρ ~ sum(s.ρ for s in species);
    ], η; kwargs...)
    sys = compose(connections, components...)
    ssys = structural_simplify(sys) # simplified system
    prob = ODEProblem(ssys, unknowns(ssys) .=> NaN, (NaN, NaN), parameters(ssys) .=> NaN; jac)
    return BackgroundSystem(sys, ssys, prob)
end

function solve(bg::BackgroundSystem, Ωγ0, Ων0, Ωc0, Ωb0; aini=1e-8, aend=1.0, solver=Vern8(), reltol=1e-8, kwargs...)
    # TODO: handle with MTK initialization when this is fixed? https://github.com/SciML/ModelingToolkit.jl/pull/2686
    # TODO: take symbolic IC map
    Ωr0 = Ωγ0 + Ων0 # TODO: gather m = c + b and r = γ + ν in the ODESystem
    ηini = aini / √(Ωr0) # analytical radiation-dominated solution # TODO: use init system
    
    prob = remake(
        bg.prob;
        tspan = (ηini, 4.0),
        u0 = [bg.ssys.g.a => aini],
        p = [bg.ssys.ph.Ω0 => Ωγ0, bg.ssys.neu.Ω0 => Ων0, bg.ssys.cdm.Ω0 => Ωc0, bg.ssys.bar.Ω0 => Ωb0]
    )

    # integrate until a == aend # TODO: just use η interval instead
    aindex = variable_index(bg.ssys, bg.ssys.g.a)
    callback = ContinuousCallback((u, _, _) -> (a = u[aindex]; a - aend), terminate!)

    return solve(prob, solver; callback, reltol, kwargs...)
end