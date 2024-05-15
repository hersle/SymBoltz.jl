struct BackgroundSystem
    # TODO: add fields like metric, gravity, species, ...?
    sys::ODESystem # unsimplified
    ssys::ODESystem # simplified
    prob::ODEProblem
end

function background_metric(; name)
    @variables a(η) ℰ(η) E(η) H(η) ℋ(η)
    @parameters H0
    a, ℰ, E, H, H0 = GlobalScope.((a, ℰ, E, H, H0)) # spacetime underlies all components # TODO: more natural way to connect them?
    eqs = [
        ℰ ~ Dη(a) / a # ℰ = ℋ/ℋ0
        E ~ ℰ / a # E = H/H0
        ℋ ~ ℰ * H0
        H ~ E * H0
    ]
    return ODESystem(eqs, η, [a, ℰ, E, ℋ, H], [H0]; name)
end

function background_gravity_GR(g; name)
    @variables ρ(η) ρcrit(η)
    eqs = [
        Dη(g.a) ~ √(8π/3 * ρ) * g.a^2 # Friedmann equation
        ρcrit ~ 3/8π * g.E^2 # critical density (H² = 8πG/3 * ρcrit)
    ]
    return ODESystem(eqs, η, [ρ, ρcrit], []; name)
end

function background_species(g, eos; name)
    @variables ρ(η) P(η)
    eqs = [
        P ~ eos(ρ) # equation of state
        Dη(ρ) ~ -3 * g.ℰ * (ρ + P) # alternative analytical solution: ρ ~ ρ0 / a^(3*(1+w))
    ]
    return ODESystem(eqs, η; name)
end
background_radiation(g; kwargs...) = background_species(g, ρ -> ρ/3; kwargs...)
background_matter(g; kwargs...) = background_species(g, ρ -> 0; kwargs...)
background_cosmological_constant(g; kwargs...) = background_species(g, ρ -> -ρ; kwargs...)

function background_ΛCDM(; name)
    @named g = background_metric()
    @named grav = background_gravity_GR(g)
    @named rad = background_radiation(g)
    @named mat = background_matter(g)
    @named de = background_cosmological_constant(g)
    return BackgroundSystem(g, grav, [rad, mat, de]; name)
end

function BackgroundSystem(g::ODESystem, grav::ODESystem, species::AbstractArray{ODESystem}; name, jac=true)
    components = [g; grav; species]
    connections = ODESystem([
        grav.ρ ~ sum(s.ρ for s in species);
    ], η; name)
    sys = compose(connections, components...)
    ssys = structural_simplify(sys) # simplified system
    prob = ODEProblem(ssys, unknowns(ssys) .=> NaN, (NaN, NaN), parameters(ssys) .=> NaN; jac)
    return BackgroundSystem(sys, ssys, prob)
end

function solve(bg::BackgroundSystem, Ωr0, Ωm0; aini=1e-8, aend=1.0, solver=KenCarp4(), reltol=1e-8, kwargs...)
    # TODO: handle with MTK initialization when this is fixed? https://github.com/SciML/ModelingToolkit.jl/pull/2686
    # TODO: take symbolic IC map
    ΩΛ0 = 1 - Ωr0 - Ωm0
    ηini = aini / √(Ωr0) # analytical radiation-dominated solution
    ρrini = 3/8π * Ωr0 / aini^4
    ρmini = 3/8π * Ωm0 / aini^3
    ρΛini = 3/8π * ΩΛ0
    
    prob = remake(bg.prob; tspan=(ηini, 3.1), u0 = [bg.ssys.g.a => aini, bg.ssys.rad.ρ => ρrini, bg.ssys.mat.ρ => ρmini, bg.ssys.de.ρ => ρΛini])

    # TODO: stop when a == aend
    aindex = variable_index(bg.ssys, bg.ssys.g.a)
    callback = ContinuousCallback((u, _, _) -> (a = u[aindex]; a - aend), terminate!)

    return solve(prob, solver; #=callback,=# reltol, kwargs...)
end