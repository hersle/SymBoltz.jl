"""
    ΛCDM(;
        lmax = 6,
        recombination = true,
        reionization = true,
        acceleration = false,
        g = metric(),
        G = general_relativity(g; acceleration),
        γ = photons(g; lmax),
        ν = massless_neutrinos(g; lmax),
        h = massive_neutrinos(g; lmax),
        c = cold_dark_matter(g; name = :c),
        b = baryons(g; recombination, reionization, name = :b),
        K = nothing,
        Λ = cosmological_constant(g),
        I = harrison_zeldovich(g; name = :I),
        name = :ΛCDM,
        kwargs...
    )

Create a ΛCDM model.
"""
function ΛCDM(;
    lmax = 6,
    recombination = true,
    reionization = true,
    acceleration = false,
    g = metric(),
    G = general_relativity(g; acceleration),
    γ = photons(g; lmax),
    ν = massless_neutrinos(g; lmax),
    h = massive_neutrinos(g; lmax),
    c = cold_dark_matter(g; name = :c),
    b = baryons(g; recombination, reionization, name = :b),
    K = nothing,
    Λ = cosmological_constant(g),
    I = harrison_zeldovich(g; name = :I),
    name = :ΛCDM,
    kwargs...
)
    species = filter(have, [γ, ν, c, b, h, K, Λ])
    pars = @parameters C
    vars = @variables fν(t) S0(t) S1(t) S_SW(t) S_ISW(t) S_Dop(t) S_pol(t)
    defs = Dict(
        C => 1//2,
        g.Ψ => 20C / (15 + 4fν) # Φ found from solving initialization system
    )
    have(ν) && have(γ) && merge!(defs, Dict(
        ν.T₀ => (4/11)^(1/3) * γ.T₀, # note: CLASS uses fudged 0.71611 ≠ (4/11)^(1/3)
        ν.Ω₀ => ν.Neff * 7/8 * (4/11)^(4/3) * γ.Ω₀
    ))
    have(h) && have(γ) && merge!(defs, Dict(
        h.T₀ => (4/11)^(1/3) * γ.T₀, # note: CLASS uses fudged 0.71611 ≠ (4/11)^(1/3)
    ))
    eqs0 = [
        G.ρ ~ sum(s.ρ for s in species)
        G.P ~ sum(s.P for s in species)
        b.rec.ρb ~ b.ρ * (H100*g.h)^2/GN # kg/m³ (convert from H₀=1 units to SI units)
        b.rec.Tγ ~ γ.T
        fν ~ sum(have(s) ? s.ρ : 0 for s in [ν, h]) / sum(s.ρ for s in [ν, h, γ] if have(s))
    ] .|> O(ϵ^0)
    eqs1 = [
        G.δρ ~ sum(s.δ * s.ρ for s in species) # total energy density perturbation
        G.δP ~ sum(s.δ * s.ρ * s.cₛ² for s in species) # total pressure perturbation
        G.Π ~ sum((s.ρ + s.P) * s.σ for s in species)
        b.θinteraction ~ -b.rec.τ̇ * 4*γ.ρ/(3*b.ρ) * (γ.θ - b.θ) # k^2*b.cₛ²*b.δ already added in baryons() # TODO: define some common interaction type, e.g. momentum transfer # TODO: would love to write something like interaction = thompson_scattering(γ, b)
        γ.τ̇ ~ b.rec.τ̇
        γ.θb ~ b.θ
        S_SW ~ b.rec.v * (γ.δ/4 + g.Ψ + γ.Π/16)
        S_ISW ~ exp(-b.rec.τ) * D(g.Ψ + g.Φ) |> expand_derivatives
        S_Dop ~ D(b.rec.v*b.u) / k |> expand_derivatives
        S_pol ~ 3/(16*k) * D(b.rec.v*γ.Π) |> expand_derivatives
        S0 ~ S_SW + S_ISW + S_Dop
        S1 ~ S_pol
    ] .|> O(ϵ^1)
    # TODO: do various initial condition types (adiabatic, isocurvature, ...) from here?
    # TODO: automatically solve for initial conditions following e.g. https://arxiv.org/pdf/1012.0569 eq. (1)?
    defs = merge(defs, Ω₀_defaults(G, species))
    description = "Standard cosmological constant and cold dark matter cosmological model"
    connections = ODESystem([eqs0; eqs1], t, vars, [pars; k]; defaults = defs, name, description)
    components = filter(!isnothing, [g; G; species; I])
    M = compose(connections, components...)
    return complete(M; flatten = false, split = false)
end

"""
    RMΛ(;
        acceleration = false,
        adiabatic = true,
        g = metric(),
        r = radiation(g; adiabatic),
        m = matter(g; adiabatic),
        Λ = cosmological_constant(g; adiabatic),
        K = nothing,
        G = general_relativity(g; acceleration),
        name = :RMΛ, kwargs...
    )

Create a simple model with pure non-interacting radiation, matter and cosmological constant.
"""
function RMΛ(;
    acceleration = false,
    adiabatic = true,
    g = metric(),
    r = radiation(g; adiabatic),
    m = matter(g; adiabatic),
    Λ = cosmological_constant(g; adiabatic),
    K = nothing,
    G = general_relativity(g; acceleration),
    name = :RMΛ, kwargs...
)
    species = filter(have, [r, m, K, Λ])
    eqs0 = [
        G.ρ ~ sum(s.ρ for s in species)
        G.P ~ sum(s.P for s in species)
    ] .|> O(ϵ^0)
    eqs1 = [
        G.δρ ~ sum(s.δ * s.ρ for s in species) # total energy density perturbation
        G.δP ~ sum(s.δ * s.ρ * s.cₛ² for s in species) # total pressure perturbation
        G.Π ~ sum((s.ρ + s.P) * s.σ for s in species)
    ] .|> O(ϵ^1)
    defs = Dict(
        g.Ψ => 20 // 15
    )
    defs = merge(defs, Ω₀_defaults(G, species))
    connections = ODESystem([eqs0; eqs1], t, [], [k]; defaults = defs, name)
    M = compose(connections, g, G, species...)
    return complete(M; flatten = false, split = false)
end

"""
    QCDM(v; name = :QCDM, kwargs...)

Create a ΛCDM model, but with the quintessence scalar field in the potential `v` as dark energy instead of the cosmological constant.
"""
function QCDM(v; name = :QCDM, kwargs...)
    M = ΛCDM()
    Q = quintessence(M.g, v)
    return ΛCDM(Λ = Q; name, kwargs...)
end

function GRΛCDM(args...; kwargs...)
    return ΛCDM(args...; kwargs...)
end

"""
    BDΛCDM(; name = :BDΛCDM, kwargs...)

Create a ΛCDM model, but with the Brans-Dicke theory of gravity instead of General Relativity.
"""
function BDΛCDM(; acceleration = false, name = :BDΛCDM, kwargs...)
    M = GRΛCDM(; Λanalytical = true, kwargs...)
    G = brans_dicke(M.g; acceleration)
    return ΛCDM(; G = G, Λanalytical = true, name, kwargs...)
end

"""
    BDRMΛ(; name = :BDRMΛ, kwargs...)

Create a RMΛ model, but with the Brans-Dicke theory of gravity instead of General Relativity.
"""
function BDRMΛ(; acceleration = false, name = :BDRMΛ, kwargs...)
    M = RMΛ(; Λanalytical = true, kwargs...)
    G = brans_dicke(M.g; acceleration)
    return RMΛ(; G = G, Λanalytical = true, name, kwargs...)
end

"""
    w0waCDM(; name = :w0waCDM, kwargs...)

Create a ΛCDM model, but with w₀wₐ-parametrized dark energy instead of the cosmological constant.
"""
function w0waCDM(; name = :w0waCDM, Λanalytical = true, kwargs...)
    M = ΛCDM(; kwargs...)
    X = w0wa(M.g; analytical = Λanalytical)
    return ΛCDM(Λ = X; name, Λanalytical, kwargs...)
end

function Ω₀_defaults(G::ODESystem, species)
    defs = Dict()
    if all(map(s -> :Ω₀ in Symbol.(parameters(s)), species)) && startswith(ModelingToolkit.description(G), "General relativity")
        s = species[end]
        if !haskey(defs, s.Ω₀)
            push!(defs, s.Ω₀ => 1 - sum(s′.Ω₀ for s′ in species[begin:end-1])) # TODO: do for all species, or do sum(s.Ω₀) ~ 1, when parameter initialization is workign
        end
    end
    return defs
end
