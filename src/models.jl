"""
    ΛCDM(;
        recombination = true,
        reionization = true,
        g = metric(),
        G = general_relativity(g),
        γ = photons(g),
        ν = massless_neutrinos(g),
        h = massive_neutrinos(g),
        c = cold_dark_matter(g; name = :c),
        b = baryons(g; recombination, reionization, name = :b),
        Λ = cosmological_constant(g),
        kwargs...
    )

Create a ΛCDM model.
"""
function ΛCDM(;
    lmax = 6,
    recombination = true,
    reionization = true,
    acceleration = false,
    Λanalytical = false,
    g = metric(),
    G = general_relativity(g; acceleration),
    γ = photons(g; lmax),
    ν = massless_neutrinos(g; lmax),
    h = massive_neutrinos(g; lmax),
    c = cold_dark_matter(g; name = :c),
    b = baryons(g; recombination, reionization, name = :b),
    Λ = cosmological_constant(g; analytical = Λanalytical),
    I = harrison_zeldovich(g; name = :I),
    name = :ΛCDM,
    kwargs...
)
    species = filter(have, [γ, ν, c, b, h, Λ])
    pars = @parameters C fν
    vars = @variables S(t) S_SW(t) S_ISW(t) S_Dop(t) S_pol(t)
    defs = Dict(
        C => 1//2,
        ϵ => 1,
        g.Ψ => 20C / (15 + 4fν) # Φ found from solving initialization system
    )
    have(ν) && have(γ) && merge!(defs, Dict(
        ν.T₀ => #=(ν.Neff/3)^(1/4) *=# (4/11)^(1/3) * γ.T₀, # TODO: CLASS uses something slightly different
        ν.Ω₀ => ν.Neff * 7/8 * (4/11)^(4/3) * γ.Ω₀
    ))
    have(h) && have(γ) && merge!(defs, Dict(
        h.T₀ => #=(ν.Neff/3)^(1/4) *=# (4/11)^(1/3) * γ.T₀, # TODO: CLASS uses something slightly different
    ))
    push!(defs, fν => sum(s.ρ for s in [ν, h] if have(s)) / sum(s.ρ for s in [ν, h, γ] if have(s)))
    eqs0 = [
        G.ρ ~ sum(s.ρ for s in species)
        G.P ~ sum(s.P for s in species)
        b.rec.ρb ~ b.ρ * g.H₀^2/GN # kg/m³ (convert from H₀=1 units to SI units)
        b.rec.Tγ ~ γ.T
    ] .|> O(ϵ^0)
    eqs1 = [
        G.δρ ~ sum(s.δ * s.ρ for s in species) # total energy density perturbation
        G.δP ~ sum(s.δ * s.ρ * s.cₛ² for s in species) # total pressure perturbation
        G.Π ~ sum((s.ρ + s.P) * s.σ for s in species)
        b.θinteraction ~ -b.rec.τ̇ * 4*γ.ρ/(3*b.ρ) * (γ.θ - b.θ) # k^2*b.cₛ²*b.δ already added in baryons() # TODO: define some common interaction type, e.g. momentum transfer # TODO: would love to write something like interaction = thompson_scattering(γ, b)
        γ.τ̇ ~ b.rec.τ̇
        γ.θb ~ b.θ
        S_SW ~ b.rec.v * (γ.δ/4 + g.Ψ)
        S_ISW ~ exp(-b.rec.τ) * D(g.Ψ + g.Φ)
        S_Dop ~ (b.rec.v̇*b.u + b.rec.v*D(b.u)) / k
        S_pol ~ 3/(4*k^2) * (b.rec.v̈*G.Π + 2*b.rec.v̇*D(γ.Π) + 0*b.rec.v*D(D(γ.Π)))
        S ~ S_SW + S_ISW + S_Dop + 0*S_pol # TODO: include all terms
    ] .|> O(ϵ^1)
    # TODO: do various IC types (adiabatic, isocurvature, ...) from here?
    initE = !Λanalytical
    if Λanalytical
        push!(defs, species[end].Ω₀ => 1 - sum(s.Ω₀ for s in species[begin:end-1])) # TODO: unsafe outside GR
    end
    description = "Standard cosmological constant and cold dark matter cosmological model"
    connections = ODESystem([eqs0; eqs1], t, vars, [pars; k]; defaults = defs, name, description)
    M = compose(connections, g, G, species..., I)
    return CosmologyModel(M; initE, kwargs...)
end

"""
    RMΛ(;
        acceleration = false,
        adiabatic = true,
        g = metric(),
        r = radiation(g; adiabatic),
        m = matter(g; adiabatic),
        Λ = cosmological_constant(g; adiabatic),
        G = general_relativity(g; acceleration),
        name = :RMΛ, kwargs...
    )

Create a simple model with pure non-interacting radiation, matter and cosmological constant.
"""
function RMΛ(;
    acceleration = false,
    adiabatic = true,
    Λanalytical = false,
    g = metric(),
    r = radiation(g; adiabatic),
    m = matter(g; adiabatic),
    Λ = cosmological_constant(g; adiabatic, analytical = Λanalytical),
    G = general_relativity(g; acceleration),
    name = :RMΛ, kwargs...
)
    species = filter(have, [r, m, Λ])
    eqs0 = [
        G.ρ ~ sum(s.ρ for s in species)
        G.P ~ sum(s.P for s in species)
    ] .|> O(ϵ^0)
    eqs1 = [
        G.δρ ~ sum(s.δ * s.ρ for s in species) # total energy density perturbation
        G.δP ~ sum(s.δ * s.ρ * s.cₛ² for s in species) # total pressure perturbation
        G.Π ~ sum((s.ρ + s.P) * s.σ for s in species)
    ] .|> O(ϵ^1)
    defs = [
        g.Ψ => 20 // 15,
        ϵ => 1
    ]
    connections = ODESystem([eqs0; eqs1], t, [], [k]; defaults = defs, name)
    M = compose(connections, g, G, species...)
    return CosmologyModel(M; spline_thermo = false)
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
