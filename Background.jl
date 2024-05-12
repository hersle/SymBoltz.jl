struct BackgroundSystem
    sys::ODESystem
    ssys::ODESystem
    prob::ODEProblem
end

function background_metric(; name)
    @variables a(η) ℰ(η) E(η)
    a, ℰ, E = GlobalScope.((a, ℰ, E)) # spacetime underlies all components # TODO: more natural way to connect them?
    eqs = [
        ℰ ~ ∂η(a) / a # in units of H0
        E ~ ℰ / a # in units of H0
    ]
    return ODESystem(eqs, η, [a, ℰ, E], []; name)
end

function background_gravity_GR(g; name)
    @variables ρ(η)
    eqs = [∂η(g.a) ~ √(8π/3*ρ) * g.a^2] # TODO: 8π/3 factor? # TODO: write E?
    return ODESystem(eqs, η; name)
end

# TODO: have species "class", sum automatically sum species.ρ in BackgroundSystem
function background_species_constant_eos(g, w; name)
    @variables ρ(η) P(η) Ω(η) ρcrit(η)
    eqs = [
        P ~ w*ρ
        ∂η(ρ) ~ -3 * g.ℰ * (ρ + P) # alternative analytical solution: ρ ~ ρ0 / a^(3*(1+w))
        Ω ~ ρ / ρcrit
    ]
    return ODESystem(eqs, η; name)
end
background_radiation(g; kwargs...) = background_species_constant_eos(g, 1//3; kwargs...)
background_matter(g; kwargs...) = background_species_constant_eos(g, 0; kwargs...)
background_cosmological_constant(g; kwargs...) = background_species_constant_eos(g, -1; kwargs...)

function BackgroundSystem(g::ODESystem, grav::ODESystem, species::AbstractArray{ODESystem}; name)
    components = [g; grav; species]
    connections = ODESystem([
        grav.ρ ~ sum(s.ρ for s in species);
        [s.ρcrit ~ grav.ρ for s in species] # TODO: 8π/(3E^2) or similar?
    ], η; name)
    sys = compose(connections, components...)
    ssys = structural_simplify(sys) # simplified system
    prob = ODEProblem(ssys, unknowns(ssys) .=> NaN, (0.0, 4.0), parameters(ssys) .=> NaN; jac=true)
    return BackgroundSystem(sys, ssys, prob)
end

# TODO: take symbolic IC map
function solve(bg::BackgroundSystem, Ωr0, Ωm0; aini=1e-8, aend=1.0, solver=KenCarp4(), reltol=1e-8, kwargs...)
    # TODO: handle with MTK initialization when this is fixed? https://github.com/SciML/ModelingToolkit.jl/pull/2686
    ΩΛ0 = 1 - Ωr0 - Ωm0
    ηini = aini / √(Ωr0) # analytical radiation-dominated solution
    ρrini = 3/(8π) * Ωr0 / aini^4
    ρmini = 3/(8π) * Ωm0 / aini^3
    ρΛini = 3/(8π) * ΩΛ0
    
    prob = remake(bg.prob; tspan=(ηini, 4.0), u0 = [bg.ssys.g.a => aini, bg.ssys.rad.ρ => ρrini, bg.ssys.mat.ρ => ρmini, bg.ssys.de.ρ => ρΛini])

    # stop when a == aend
    aindex = variable_index(bg.ssys, bg.ssys.g.a)
    callback = ContinuousCallback((u, _, _) -> (a = u[aindex]; a - aend), terminate!)

    return solve(prob, solver; callback, reltol, kwargs...) # using KenCarp4 here leads to difference in AD vs FD
end