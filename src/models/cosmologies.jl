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
    pars = @parameters begin
        C, [description = "Initial conditions integration constant"]
        τ0, [description = "Conformal time today"]
    end
    vars = @variables begin
        χ(τ), [description = "Conformal lookback time from today"]
        fν(τ), [description = "Neutrino-to-radiation density fraction"]
        ST0(τ, k), [description = "Temperature source function weighted against spherical Bessel function"]
        ST1(τ, k), [description = "Temperature source function weighted against spherical Bessel function 1st derivative"]
        ST0_SW(τ, k),
        ST0_ISW(τ, k), ST1_ISW(τ, k),
        ST0_Doppler(τ, k), ST1_Doppler(τ, k),
        ST0_polarization(τ, k), ST1_polarization(τ, k), ST2_polarization(τ, k)
    end
    defs = Dict(
        C => 1//2,
        τ0 => NaN
    )
    ics = [
        g.Ψ ~ 20C / (15 + 4fν) # Φ found from solving initialization system
    ]
    eqs = Ω₀_eqs(G, species)
    have(ν) && have(γ) && push!(eqs,
        ν.T₀ ~ (4/11)^(1/3) * γ.T₀, # note: CLASS uses fudged 0.71611 ≠ (4/11)^(1/3)
        ν.Ω₀ ~ ν.Neff * 7/8 * (4/11)^(4/3) * γ.Ω₀,
    )
    have(h) && have(γ) && push!(eqs,
        h.T₀ ~ (4/11)^(1/3) * γ.T₀, # note: CLASS uses fudged 0.71611 ≠ (4/11)^(1/3)
    )
    append!(eqs, [
        G.ρ ~ sum(s.ρ for s in species)
        G.P ~ sum(s.P for s in species)
        b.rec.ρb ~ b.ρ * (H100*g.h)^2/GN # kg/m³ (convert from H₀=1 units to SI units)
        b.rec.Tγ ~ γ.T
        fν ~ sum(have(s) ? s.ρ : 0 for s in [ν, h]) / sum(s.ρ for s in [ν, h, γ] if have(s))
        χ ~ τ0 - τ
    ])
    append!(eqs, [
        G.δρ ~ sum(s.δ * s.ρ for s in species) # total energy density perturbation
        G.δP ~ sum(s.δ * s.ρ * s.cₛ² for s in species) # total pressure perturbation
        G.Π ~ sum((1 + s.w) * s.ρ * s.σ for s in species) # TODO: factor 2/3 or 3/2? See e.g. https://arxiv.org/pdf/astro-ph/9506072 bottom of page 10? Check all models.
        b.θinteraction ~ -b.rec.κ̇ * 4*γ.ρ/(3*b.ρ) * (γ.θ - b.θ) # k^2*b.cₛ²*b.δ already added in baryons() # TODO: define some common interaction type, e.g. momentum transfer # TODO: would love to write something like interaction = thompson_scattering(γ, b)
        γ.κ̇ ~ b.rec.κ̇
        γ.θb ~ b.θ

        ST0_SW ~ b.rec.v * (γ.δ/4 + g.Ψ + γ.Π/16)
        ST0_ISW ~ exp(-b.rec.κ) * (g.Ψ̇ + g.Φ̇)
        ST0_Doppler ~ D(b.rec.v*b.u) / k |> expand_derivatives
        ST1_Doppler ~ b.rec.v*b.u
        ST0_polarization ~ 3/(16*k^2) * D(D(b.rec.v*γ.Π)) |> expand_derivatives
        ST1_polarization ~ 3/(16*k) * D(b.rec.v*γ.Π) |> expand_derivatives
        ST2_polarization ~ 3/16 * b.rec.v*γ.Π
        ST0 ~ ST0_SW + ST0_ISW + ST0_Doppler + ST0_polarization
        ST1 ~ 0
    ])
    # TODO: do various initial condition types (adiabatic, isocurvature, ...) from here?
    # TODO: automatically solve for initial conditions following e.g. https://arxiv.org/pdf/1012.0569 eq. (1)?
    description = "Standard cosmological constant and cold dark matter cosmological model"
    connections = System(eqs, τ, vars, [pars; k]; defaults = defs, initialization_eqs = ics, name, description)
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
        I = harrison_zeldovich(g; name = :I),
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
    I = harrison_zeldovich(g; name = :I),
    name = :RMΛ, kwargs...
)
    vars = @variables begin
        χ(τ), [description = "Conformal lookback time from today"]
    end
    pars = @parameters begin
        τ0, [description = "Conformal time today"]
    end
    species = filter(have, [r, m, K, Λ])
    eqs = [
        G.ρ ~ sum(s.ρ for s in species)
        G.P ~ sum(s.P for s in species)
        χ ~ τ0 - τ

        G.δρ ~ sum(s.δ * s.ρ for s in species) # total energy density perturbation
        G.δP ~ sum(s.δ * s.ρ * s.cₛ² for s in species) # total pressure perturbation
        G.Π ~ sum((s.ρ + s.P) * s.σ for s in species)
    ]
    defs = Dict(
        g.Ψ => 20 // 15,
        τ0 => NaN
    )
    append!(eqs, Ω₀_eqs(G, species)) # parameter equations
    connections = System(eqs, τ, vars, [pars; k]; defaults = defs, name)
    components = filter(!isnothing, [g; G; species; I])
    M = compose(connections, components...)
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

function Ω₀_eqs(G::System, species)
    if all(map(s -> :Ω₀ in Symbol.(full_parameters(s)), species)) && startswith(ModelingToolkit.description(G), "General relativity")
        return [species[end].Ω₀ ~ 1 - sum(s′.Ω₀ for s′ in species[begin:end-1])] # TODO: do for all species, or do sum(s.Ω₀) ~ 1, when parameter initialization is workign
    else
        return Equation[]
    end
end
