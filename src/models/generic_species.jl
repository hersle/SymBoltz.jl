"""
    species_constant_eos(g, w, ẇ = 0, σ = 0; analytical = true, adiabatic = false, interact = false, continuity_pressure = true, name = :s, kwargs...)

Create a symbolic component for a particle species with equation of state ``w = P/ρ`` in the spacetime with the metric `g`,
with the generic species evolution equations in [Ma & Bertschinger](https://arxiv.org/abs/astro-ph/9506072).

# Keyword arguments

- `analytical`: whether the background continuity equation is integrated analytically (requires constant ``w``).
- `adiabatic`: whether to set ``cₛ² = cₐ²`` and ``δP = cₐ² δρ``; otherwise the caller must specify ``cₛ²`` and ``δP``.
- `interact`: whether the momentum transfer ``f`` from other species is left unspecified so it can be set elsewhere; otherwise it is set to ``f = 0``.
- `continuity_pressure`: whether the perturbed continuity equation includes the pressure term.
"""
function species_constant_eos(g, w, ẇ = 0, σ = 0; analytical = true, adiabatic = false, interact = false, continuity_pressure = true, name = :s, kwargs...)
    _w, _σ = w, σ # w and σ are redefined as symbolic variables below
    @assert ẇ == 0 && _σ == 0 # TODO: relax (need to include in ICs)
    if analytical
        pars = @parameters Ω₀, [description = "Reduced background density today"]
    else
        pars = @parameters ρᵢ, [description = "Initial background density"]
    end
    vars = @variables begin
        w(τ), [description = "Equation of state"]
        ρ(τ), [description = "Background density"]
        P(τ), [description = "Background pressure"]
        Ω(τ), [description = "Reduced background density"]
        cₛ²(τ), [description = "Speed of sound squared"]
        cₐ²(τ), [description = "Adiabatic speed of sound squared"]
        δP(τ, k), [description = "Pressure perturbation"]
        δ(τ, k), [description = "Overdensity (gauge-dependent)"]
        Δ(τ, k), [description = "Overdensity (gauge-independent)"]
        θ(τ, k), [description = "Velocity divergence"]
        f(τ, k), [description = "Momentum transfer from other species"]
        σ(τ, k), [description = "Shear stress"]
        u(τ, k), [description = "Velocity"]
        u̇(τ, k), [description = "Velocity derivative"]
    end
    n = 3 * (1 + _w)
    n = Symbolics.symbolic_type(_w) == Symbolics.NotSymbolic() && isinteger(n) ? Int(n) : n
    if analytical
        eqs = [
            Ω ~ Ω₀ / g.a^n
            ρ ~ 3/(8*Num(π)) * Ω
        ]
        initial_conditions = []
    else
        eqs = [
            D(ρ) ~ -3 * g.ℋ * (ρ + P)
            Ω ~ 8π/3 * ρ
        ]
        initial_conditions = [ρ => ρᵢ]
    end
    append!(eqs, [
        w ~ _w
        P ~ w * ρ

        cₐ² ~ w - (iszero(ẇ) ? 0 : ẇ/(3*g.ℋ*(1+w))) # Ṗ/ρ̇ (branch avoids 0/0 for w = -1)
        D(δ) ~ -(1+w)*(θ-3*D(g.Φ)) - (continuity_pressure ? 3*g.ℋ*(δP/ρ-w*δ) : 0) # Bertschinger & Ma (30)
        D(θ) ~ -g.ℋ*(1-3*w)*θ - ẇ/(1+w)*θ + k^2*δP/((1+w)*ρ) - k^2*σ + k^2*g.Ψ + f/(ρ+P) # Bertschinger & Ma (30) with additional momentum exchange (ρ+P)θ′ = … + f
        Δ ~ δ + 3*g.ℋ*(1+w)*θ/k^2
        u ~ θ / k
        u̇ ~ D(u)
        σ ~ _σ
    ])
    adiabatic && push!(eqs, cₛ² ~ cₐ², δP ~ cₐ²*ρ*δ)
    !interact && push!(eqs, f ~ 0)
    ieqs = [
        δ ~ -3//2 * (1+w) * g.Ψ # adiabatic: δᵢ/(1+wᵢ) == δⱼ/(1+wⱼ) (https://cmb.wintherscoming.no/theory_initial.php#adiabatic) # TODO: match CLASS with higher-order (for photons)? https://github.com/lesgourg/class_public/blob/22b49c0af22458a1d8fdf0dd85b5f0840202551b/source/perturbations.c#L5631-L5632
        θ ~ 1//2 * (k^2*τ) * g.Ψ # τ ≈ 1/ℋ # TODO: include σ ≠ 0 # solve u′ + ℋ(1-3w)u = w/(1+w)*kδ + kΨ with Ψ=const, IC for δ, Φ=-Ψ, ℋ=H₀√(Ωᵣ₀)/a after converting ′ -> d/da by gathering terms with u′ and u in one derivative using the trick to multiply by exp(X(a)) such that X′(a) will "match" the terms in front of u
    ]
    return System(eqs, τ, vars, [pars; k]; initialization_eqs = ieqs, initial_conditions, name, kwargs...)
end

"""
    matter(g; name = :m, kwargs...)

Create a particle species for matter (with equation of state `w ~ 0`) in the spacetime with metric `g`.
"""
function matter(g; name = :m, kwargs...)
    description = "Matter"
    return species_constant_eos(g, 0; adiabatic = true, name, description, kwargs...)
end

"""
    radiation(g; name = :r, kwargs...)

Create a particle species for radiation (with equation of state `w ~ 1/3`) in the spacetime with metric `g`.
"""
function radiation(g; name = :r, kwargs...)
    r = species_constant_eos(g, 1//3; name, kwargs...) |> complete
    pars = @parameters begin
        T₀, [description = "Temperature today (in K)"]
    end
    vars = @variables begin
        T(τ), [description = "Temperature"] # TODO: define in constant_eos? https://physics.stackexchange.com/questions/650508/whats-the-relation-between-temperature-and-scale-factor-for-arbitrary-eos-1
    end
    eqs = [T ~ T₀ / g.a]
    description = "Radiation"
    return extend(r, System(eqs, τ, vars, pars; name); description)
end

"""
    effective_species(g, species; effective_name = "", kwargs...)

Create an effective "read-only" species for several given `species` with metric `g`.
Additive properties (like ``ρ``, ``P``, ``δρ`` and ``δP``) are summed, and used to express non-additive properties (like ``w`` and ``δ``).
"""
function effective_species(g, species; effective_name = "", kwargs...)
    scope = ParentScope
    pars = @parameters begin
        Ω₀ = scope(sum(s.Ω₀ for s in species)), [description = "Reduced background density today"]
    end
    vars = @variables begin
        w(τ), [description = "Equation of state"]
        ρ(τ), [description = "Background density"]
        P(τ), [description = "Background pressure"]
        δ(τ, k), [description = "Overdensity (gauge-dependent)"]
        Δ(τ, k), [description = "Overdensity (gauge-independent)"]
        θ(τ, k), [description = "Velocity divergence"]
        δP(τ, k), [description = "Pressure perturbation"]
    end
    eqs = [
        ρ ~ scope(sum(s.ρ for s in species))
        P ~ scope(sum(s.P for s in species))
        w ~ P / ρ
        δ ~ scope(sum(s.δ*s.ρ for s in species)) / ρ
        θ ~ scope(sum((1+s.w)*s.ρ*s.θ for s in species)) / (ρ + P)
        Δ ~ scope(sum(s.ρ*s.Δ for s in species)) / ρ
        δP ~ scope(sum(s.δP for s in species))
    ]
    description = "Effective species for " * join(nameof.(species), '+')
    if !isempty(effective_name)
        description = "$effective_name ($(lowercasefirst(description)))"
    end
    return System(eqs, τ, vars, pars; description, kwargs...)
end
