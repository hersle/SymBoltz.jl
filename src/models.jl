"""
    ΛCDM(;
        recombination = true,
        g = metric(),
        G = general_relativity(g),
        γ = photons(g),
        ν = massless_neutrinos(g),
        h = massive_neutrinos(g),
        c = cold_dark_matter(g; name = :c),
        b = baryons(g; recombination, name = :b),
        Λ = cosmological_constant(g),
        kwargs...
    )

Create a ΛCDM model.
"""
function ΛCDM(;
    lmax = 6,
    recombination = true,
    acceleration = false,
    Λanalytical = false,
    g = metric(),
    G = general_relativity(g; acceleration),
    γ = photons(g; lmax),
    ν = massless_neutrinos(g; lmax),
    h = massive_neutrinos(g; lmax),
    c = cold_dark_matter(g; name = :c),
    b = baryons(g; recombination, name = :b),
    Λ = cosmological_constant(g; analytical = Λanalytical),
    name = :ΛCDM,
    kwargs...
)
    species = [γ, ν, c, b, h, Λ]
    pars = @parameters C fν
    defs = Dict(
        ν.T0 => (ν.Neff/3)^(1/4) * (4/11)^(1/3) * γ.T0,
        ν.Ω0 => (ν.Neff/3) * 7/8 * (4/11)^(4/3) * γ.Ω0,
        h.T0 => (ν.Neff/3)^(1/4) * (4/11)^(1/3) * γ.T0, # same as for massless neutrinos # TODO: are the massive neutrino density parameters correct?
        h.Ω0_massless => 7/8 * (h.T0/γ.T0)^4 * γ.Ω0, # Ω0 for corresponding massless neutrinos # TODO: reconcile with class? https://github.com/lesgourg/class_public/blob/ae99bcea1cd94994228acdfaec70fa8628ae24c5/source/background.c#L1561
        k => NaN, # make background shut up # TODO: avoid
        fν => ν.ρ0 / (ν.ρ0 + ν.ρ0),
        C => 0.48, # TODO: why does ≈ 0.48 give better agreement with CLASS? # TODO: phi set here? https://github.com/lesgourg/class_public/blob/ae99bcea1cd94994228acdfaec70fa8628ae24c5/source/perturbations.c#L5713
        g.Ψ => 20C / (15 + 4fν), # Φ found from solving initialization system # TODO: is this correct when having both massless and massive neutrinos?
        ϵ => 1 # TODO: remove
    )
    eqs0 = [
        G.ρ ~ sum(s.ρ for s in species) # TODO: only if G has ρ
        G.P ~ sum(s.P for s in species) # TODO: only if G has P
        b.rec.ρb ~ b.ρ * g.H0^2/GN # kg/m³ (convert from H0=1 units to SI units)
        b.rec.Tγ ~ γ.T
    ] .|> O(ϵ^0)
    eqs1 = [
        G.δρ ~ sum(s.δ * s.ρ for s in species) # total energy density perturbation
        G.δP ~ sum(s.δ * s.ρ * s.cs² for s in species) # total pressure perturbation
        G.Π ~ sum((s.ρ + s.P) * s.σ for s in species)
        b.θinteraction ~ k^2*b.cs²*b.δ + -b.rec.τ̇ * 4*γ.ρ/(3*b.ρ) * (γ.θ - b.θ) # TODO: define some common interaction type, e.g. momentum transfer # TODO: would love to write something like interaction = thompson_scattering(γ, b)
        γ.τ̇ ~ b.rec.τ̇
        γ.θb ~ b.θ
    ] .|> O(ϵ^1)
    # TODO: do various IC types (adiabatic, isocurvature, ...) from here?
    initE = !Λanalytical
    if Λanalytical
        push!(defs, species[end].Ω0 => 1 - sum(s.Ω0 for s in species[begin:end-1])) # TODO: unsafe outside GR
    end
    connections = ODESystem([eqs0; eqs1], t, [], [pars; k]; defaults = defs, name)
    M = compose(connections, g, G, species...)
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
    species = [r, m, Λ]
    eqs0 = [
        G.ρ ~ sum(s.ρ for s in species) # TODO: only if G has ρ
        G.P ~ sum(s.P for s in species) # TODO: only if G has P
    ] .|> O(ϵ^0)
    eqs1 = [
        G.δρ ~ sum(s.δ * s.ρ for s in species) # total energy density perturbation
        G.δP ~ sum(s.δ * s.ρ * s.cs² for s in species) # total pressure perturbation
        G.Π ~ sum((s.ρ + s.P) * s.σ for s in species)
    ] .|> O(ϵ^1)
    defs = [
        g.Ψ => 20 / 15, # TODO: put to what?
        ϵ => 1 # TODO: remove
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
