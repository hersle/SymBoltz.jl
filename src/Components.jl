ϵ = only(GlobalScope.(@parameters ϵ)) # perturbative expansion parameter

"""
    metric(; name = :g, kwargs...)

Create a symbolic component for the perturbed FLRW spacetime metric in the conformal Newtonian gauge with signature

    a^2 * diag(-1-2*ϵ*Ψ, +1 - 2*ϵ*Φ, +1 - 2*ϵ*Φ, +1 - 2*ϵ*Φ)
"""
function metric(; name = :g, kwargs...)
    vars = a, ℰ, E, H, ℋ, Φ, Ψ, g11, g22, h11, h22 = GlobalScope.(@variables a(t) ℰ(t) E(t) H(t) ℋ(t) Φ(t) Ψ(t) g11(t) g22(t) h11(t) h22(t))
    pars = H0, h = GlobalScope.(@parameters H0 h)
    defs = [
        H0 => H100 * h
        h => H0 / H100
    ]
    return ODESystem([
        D(a) ~ ℰ * a # ℰ = ℋ/ℋ0 = ℋ/H0
        ℰ ~ a * E # E = H/H0
        ℋ ~ ℰ * H0
        H ~ E * H0
        g11 ~ -a^2
        g22 ~ +a^2
        h11 * ϵ ~ -2 * a^2 * Ψ * ϵ
        h22 * ϵ ~ -2 * a^2 * Φ * ϵ
    ], t, vars, pars; defaults=defs, name, kwargs...)
end

"""
    general_relativity(g; name = :G, kwargs...)

Create a symbolic component for the general relativistic (GR) theory of gravity in the spacetime with the metric `g`.
"""
function general_relativity(g; name = :G, kwargs...)
    vars = @variables ρ(t) ρcrit(t) δρ(t) Π(t)
    eqs0 = [
        g.E ~ √(8π/3 * ρ) # Friedmann equation
        ρcrit ~ 3/8π * g.E^2 # critical density (H² = 8πG/3 * ρcrit)
    ] .|> O(ϵ^0)
    eqs1 = [
        D(g.Φ) ~ -4π/3*g.a^2/g.ℰ*δρ - k^2/(3*g.ℰ)*g.Φ - g.ℰ*g.Ψ
        k^2 * (g.Φ - g.Ψ) ~ 12π * g.a^2 * Π
    ] .|> O(ϵ^1)
    guesses = [ρ => 1]
    return ODESystem([eqs0; eqs1], t, vars, []; guesses, name, kwargs...)
end

"""
    species_constant_eos(g, w, cs² = w, ẇ = 0, _σ = 0; θinteract = false, kwargs...)

Create a symbolic component for a particle species with equation of state `w ~ P/ρ` in the spacetime with the metric `g`.
"""
function species_constant_eos(g, w, cs² = w, ẇ = 0, _σ = 0; analytical = true, θinteract = false, kwargs...)
    @assert ẇ == 0 && _σ == 0 # TODO: relax (need to include in ICs)
    pars = analytical ? (@parameters ρ0 Ω0) : [] # TODO: split pars0, pars1
    vars = @variables ρ(t) P(t) δ(t) θ(t) Δ(t) θinteraction(t) σ(t) # TODO: split vars0, vars1
    eqs0 = [
        P ~ w * ρ # equation of state
        analytical ? (ρ ~ ρ0 * g.a^(-3*(1+w))) : (D(ρ) ~ -3 * g.ℰ * (ρ + P)) # alternative derivative: D(ρ) ~ -3 * g.ℰ * (ρ + P)
    ] .|> O(ϵ^0)
    eqs1 = [
        D(δ) ~ -(1+w)*(θ-3*D(g.Φ)) - 3*g.ℰ*(cs²-w)*δ # Bertschinger & Ma (30) with Φ -> -Φ; or Baumann (4.4.173) with Φ -> -Φ
        D(θ) ~ -g.ℰ*(1-3*w)*θ - ẇ/(1+w)*θ + cs²/(1+w)*k^2*δ - k^2*σ + k^2*g.Ψ + θinteraction # Bertschinger & Ma (30) with θ = kv
        Δ ~ δ + 3(1+w) * g.ℰ/θ # Baumann (4.2.144) with v -> -u
        σ ~ _σ
    ] .|> O(ϵ^1)
    ics1 = [
        δ ~ -3/2 * (1+w) * g.Ψ # adiabatic: δᵢ/(1+wᵢ) == δⱼ/(1+wⱼ) (https://cmb.wintherscoming.no/theory_initial.php#adiabatic)
        θ ~ 1/2 * (k^2*t) * g.Ψ # # TODO: include σ ≠ 0 # solve u′ + ℋ(1-3w)u = w/(1+w)*kδ + kΨ with Ψ=const, IC for δ, Φ=-Ψ, ℋ=H₀√(Ωᵣ₀)/a after converting ′ -> d/da by gathering terms with u′ and u in one derivative using the trick to multiply by exp(X(a)) such that X′(a) will "match" the terms in front of u
    ] .|> O(ϵ^1)
    defs = analytical ? [ρ0 => 3/8π * Ω0] : Dict()
    !θinteract && push!(eqs1, (θinteraction ~ 0) |> O(ϵ^1))
    return ODESystem([eqs0; eqs1], t, vars, pars; initialization_eqs=ics1, defaults=defs, kwargs...)
end

"""
    matter(g; name = :m, kwargs...)

Create a particle species for matter (with equation of state `w ~ 0`) in the spacetime with metric `g`.
"""
function matter(g; name = :m, kwargs...)
    return species_constant_eos(g, 0; name, kwargs...)
end

"""
    radiation(g; name = :r, kwargs...)

Create a particle species for radiation (with equation of state `w ~ 1/3`) in the spacetime with metric `g`.
"""
function radiation(g; name = :r, kwargs...)
    r = species_constant_eos(g, 1//3; name, kwargs...)
    pars = @parameters T0
    vars = @variables T(t) # TODO: define in constant_eos? https://physics.stackexchange.com/questions/650508/whats-the-relation-between-temperature-and-scale-factor-for-arbitrary-eos-1
    return extend(r, ODESystem([T ~ T0 / g.a], t, vars, pars; name, kwargs...))
end

"""
    cosmological_constant(g; name = :Λ, kwargs...)

Create a particle species for the cosmological constant (with equation of state `w ~ -1`) in the spacetime with metric `g`.
"""
function cosmological_constant(g; name = :Λ, analytical = false, kwargs...)
    Λ = species_constant_eos(g, -1; name, analytical, kwargs...) |> thermodynamics |> complete # discard ill-defined perturbations
    vars = @variables δ(t) θ(t) σ(t)
    return extend(Λ, ODESystem([δ ~ 0, θ ~ 0, σ ~ 0] .|> O(ϵ^1), t, vars, []; name, kwargs...)) # manually set perturbations to zero
end

"""
    cosmological_constant(g; name = :Λ, kwargs...)

Create a particle species for photons in the spacetime with metric `g`.
"""
function photons(g; polarization=true, lmax=6, name = :γ, kwargs...)
    γ = radiation(g; name, kwargs...) |> thermodynamics |> complete # prevent namespacing in extension below

    vars = @variables F0(t) F(t)[1:lmax] δ(t) θ(t) σ(t) τ̇(t) θb(t) Π(t) G0(t) G(t)[1:lmax]
    defs = [
        γ.T0 => (15/π^2 * γ.ρ0 * g.H0^2/GN * ħ^3*c^5)^(1/4) / kB
    ]
    eqs1 = [
        D(F0) ~ -k*F[1] + 4*D(g.Φ)
        D(F[1]) ~ k/3*(F0-2*F[2]+4*g.Ψ) - 4/3 * τ̇/k * (θb - θ)
        D(F[2]) ~ 2/5*k*F[1] - 3/5*k*F[3] + 9/10*τ̇*F[2] - 1/10*τ̇*(G0+G[2])
        [D(F[l]) ~ k/(2*l+1) * (l*F[l-1] - (l+1)*F[l+1]) + τ̇*F[l] for l in 3:lmax-1]... # TODO: Π in last term here?
        D(F[lmax]) ~ k*F[lmax-1] - (lmax+1) / t * F[lmax] + τ̇ * F[lmax] # TODO: assumes lmax ≥ ???
        δ ~ F0
        θ ~ 3/4*k*F[1]
        σ ~ F[2]/2
        Π ~ F[2] + G0 + G[2]
    ] .|> O(ϵ^1)
    ics1 = [
        δ ~ -2 * g.Ψ # Dodelson (7.89)
        θ ~ 1/2 * (k^2*t) * g.Ψ # Dodelson (7.95)
        F[2] ~ 0 # (polarization ? -8/15 : -20/45) * k/dτ * Θ[1], # depends on whether polarization is included
        [F[l] ~ 0 #=-l/(2*l+1) * k/dτ * Θ[l-1]=# for l in 3:lmax]...
    ] .|> O(ϵ^1)
    if polarization
        union!(eqs1, [
            D(G0) ~ k * (-G[1]) - τ̇ * (-G0 + Π/2)
            D(G[1]) ~ k/3 * (G0 - 2*G[2]) - τ̇ * (-G[1])
            [D(G[l]) ~ k/(2*l+1) * (l*G[l-1] - (l+1)*G[l+1]) - τ̇ * (-G[l] + Π/10*δkron(l,2)) for l in 2:lmax-1]... # TODO: collect all equations here once G[0] works
            D(G[lmax]) ~ k*G[lmax-1] - (lmax+1) / t * G[lmax] + τ̇ * G[lmax]
        ] .|> O(ϵ^1))
        union!(ics1, [
            G0 ~ 0 #5/4 * Θ[2],
            G[1] ~ 0 #-1/4 * k/dτ * Θ[2],
            G[2] ~ 0 #1/4 * Θ[2],
            [G[l] ~ 0 #=-l/(2*l+1) * k/dτ * ΘP[l-1]=# for l in 3:lmax]...    
        ] .|> O(ϵ^1))
    else
        union!(eqs1, [G0 ~ 0, collect(G .~ 0)...] .|> O(ϵ^1)) # pin to zero
    end
    return extend(γ, ODESystem(eqs1, t, vars, []; initialization_eqs=ics1, defaults=defs, name, kwargs...))
end

"""
    massless_neutrinos(g; lmax=6, name = :ν, kwargs...)

Create a particle species for massless neutrinos in the spacetime with metric `g`.
"""
function massless_neutrinos(g; lmax=6, name = :ν, kwargs...)
    ν = radiation(g; name, kwargs...) |> thermodynamics |> complete

    vars = @variables F0(t) F(t)[1:lmax+1] δ(t) θ(t) σ(t)
    pars = @parameters Neff
    defs = [
        ν.T0 => NaN # TODO: use this
        Neff => 3.046
    ]
    eqs1 = [
        D(F0) ~ -k*F[1] + 4*D(g.Φ)
        D(F[1]) ~ k/3*(F0-2*F[2]+4*g.Ψ)
        [D(F[l]) ~ k/(2*l+1) * (l*F[l-1] - (l+1)*F[l+1]) for l in 2:lmax]...
        F[lmax+1] ~ (2*lmax+1) / (k*t) * F[lmax] - F[lmax-1]
        δ ~ F0
        θ ~ 3/4*k*F[1]
        σ ~ F[2]/2
    ] .|> O(ϵ^1)
    ics1 = [
        δ ~ -2 * g.Ψ # adiabatic: δᵢ/(1+wᵢ) == δⱼ/(1+wⱼ) (https://cmb.wintherscoming.no/theory_initial.php#adiabatic)
        θ ~ 1/2 * (k^2*t) * g.Ψ
        σ ~ 1/15 * (k*t)^2 * g.Ψ # TODO: how to set ICs consistently with Ψ, Π and Θν2?
        [F[l] ~ 0 #=1/(2*l+1) * k*t * Θ[l-1]=# for l in 3:lmax]...
    ] .|> O(ϵ^1)
    return extend(ν, ODESystem(eqs1, t, vars, pars; initialization_eqs=ics1, defaults=defs, name, kwargs...))
end

# TODO: use vector equations and simplify loops
"""
    massive_neutrinos(g; nx=5, lmax=4, name = :h, kwargs...)

Create a particle species for massive neutrinos in the spacetime with metric `g`.
"""
function massive_neutrinos(g; nx=5, lmax=4, name = :h, kwargs...)
    pars = @parameters Ω0_massless ρ0_massless Ω0 ρ0 m T0 y0
    vars = @variables ρ(t) T(t) y(t) P(t) w(t) δ(t) σ(t) θ(t) ψ0(t)[1:nx] ψ(t)[1:nx,1:lmax+1]

    f0(x) = 1 / (exp(x) + 1) # TODO: why not exp(E)?
    dlnf0_dlnx(x) = -x / (1 + exp(-x))
    x, W = gauss(x -> x^2 * f0(x), nx, 0.0, 1e3) # Gaussian quadrature weights, reduced momentum bins x = q*c / (kB*T0) # these points give accurate integral for Iρmν in the background, at least # TODO: ok for perturbations? # TODO: also include common x^2 factor in weighting?
    ∫dx_x²_f0(f) = sum(collect(f) .* W) # a function that approximates the weighted integral ∫dx*x^2*f(x)*f0(x)

    E(x, y) = √(x^2 + y^2)
    Iρ(y) = ∫dx_x²_f0(@. E(x, y)) # Iρ(0) = 7π^4/120
    IP(y) = ∫dx_x²_f0(@. x^2 / E(x, y)) # IP(0) = Iρ(0)
    
    eqs0 = [
        T ~ T0 / g.a # TODO: move into radiation?
        y ~ m*c^2 / (kB*T)
        ρ ~ ρ0_massless/g.a^4 * Iρ(y) / Iρ(0) # have ρ = Cρ * Iρ(y) / a⁴, so Cρ = ρ0 * 1⁴ / Iρ(y0) # TODO: div by Iρ(0) or Iρ(y0)?
        P ~ 1/3 * ρ0_massless/g.a^4 * IP(y) / Iρ(0) # have P = CP * IP(y) / a⁴, and in the early universe Iρ(y→0) → IP(y→0) and P/ρ = CP * IP(y) / (Cρ * Iρ(y)) → CP/Cρ → 1/3, so CP = Cρ/3 # TODO: div by Iρ(0) or Iρ(y0)?
        w ~ P / ρ
    ] .|> O(ϵ^0)
    eqs1 = [
        δ ~ ∫dx_x²_f0(@. E(x, y)*ψ0) / ∫dx_x²_f0(@. E(x, y))
        # TODO: θ
        σ ~ (2/3) * ∫dx_x²_f0(@. x^2/E(x,y)*ψ[:,2]) / (∫dx_x²_f0(@. E(x,y)) + 1/3*∫dx_x²_f0(@. x^2/E(x,y)))
    ] .|> O(ϵ^1)
    defs = [
        Ω0 => Ω0_massless * Iρ(y0) / Iρ(0) # ≈ Ω0_massless * (3ζ(3)/2)/(7π^4/120) * y0 for y0 → ∞
        ρ0 => 3/8π * Ω0
        ρ0_massless => 3/8π * Ω0_massless
        m => 0.02 * eV/c^2 # one massive neutrino with this mass # TODO: specify by user
        y0 => m*c^2 / (kB*T0)
    ]
    ics1 = []
    for i in 1:nx
        union!(eqs1, [
            D(ψ0[i]) ~ -k * x[i]/E(x[i],y) * ψ[i,1] - D(g.Φ) * dlnf0_dlnx(x[i])
            D(ψ[i,1]) ~ k/3 * x[i]/E(x[i],y) * (ψ0[i] - 2*ψ[i,2]) - k/3 * E(x[i],y)/x[i] * g.Ψ * dlnf0_dlnx(x[i])
            [D(ψ[i,l]) ~ k/(2*l+1) * x[i]/E(x[i],y) * (l*ψ[i,l-1] - (l+1)*ψ[i,l+1]) for l in 2:lmax]...
            ψ[i,lmax+1] ~ (2*lmax+1) * E(x[i],y)/x[i] * ψ[i,lmax] / (k*t) - ψ[i,lmax-1]
        ] .|> O(ϵ^1))
        union!(ics1, [
            ψ0[i] ~ -1/4 * (-2*g.Ψ) * dlnf0_dlnx(x[i])
            ψ[i,1] ~ -1/(3*k) * E(x[i],y)/x[i] * (1/2*(k^2*t)*g.Ψ) * dlnf0_dlnx(x[i])
            ψ[i,2] ~ -1/2 * (1/15*(k*t)^2*g.Ψ) * dlnf0_dlnx(x[i])
            [ψ[i,l] ~ 0 for l in 3:lmax] # TODO: full ICs
        ] .|> O(ϵ^1))
    end
    return ODESystem([eqs0; eqs1], t, vars, pars; initialization_eqs=ics1, defaults=defs, name, kwargs...)
end

"""
    baryons(g; recombination=true, name = :b, kwargs...)

Create a particle species for baryons in the spacetime with metric `g`.
"""
function baryons(g; recombination=true, name = :b, kwargs...)
    b = matter(g; θinteract=true, name, kwargs...)
    if recombination
        @named rec = thermodynamics_recombination_recfast(g)
    else
        vars = @variables dτ(t) ρb(t) Tγ(t)
        @named rec = ODESystem([dτ ~ 0], t, vars, [])
    end
    b = compose(b, rec)
    return b
end

function background(sys)
    sys = thermodynamics(sys)
    sys = extend(sys, ODESystem([], t; initialization_eqs = [sys.g.E ~ 1], name = sys.name)) # initialize with H == H0 today
    return replace(sys, sys.b.rec => ODESystem([], t; name = :rec))
end

function thermodynamics(sys)
    return transform((sys, _) -> extract_order(sys, [0]), sys)
end

function perturbations(sys; spline_thermo=true)
    if spline_thermo
        @named rec = thermodynamics_recombination_splined()
        sys = replace(sys, sys.b.rec => rec)
    end
    return transform((sys, _) -> extract_order(sys, [0, 1]), sys)
end

"""
    ΛCDM(;
        recombination = true,
        g = metric(),
        G = general_relativity(g),
        γ = photons(g),
        ν = massless_neutrinos(g),
        h = massive_neutrinos(g),
        c = matter(g; name = :c),
        b = baryons(g; recombination, name = :b),
        Λ = cosmological_constant(g),
        kwargs...
    )

Create a complete ΛCDM model.
"""
function ΛCDM(;
    recombination = true,
    g = metric(),
    G = general_relativity(g),
    γ = photons(g),
    ν = massless_neutrinos(g),
    h = massive_neutrinos(g),
    c = matter(g; name = :c),
    b = baryons(g; recombination, name = :b),
    Λ = cosmological_constant(g),
    name = :ΛCDM,
    kwargs...
)
    species = [γ, ν, c, b, h, Λ]
    pars = @parameters C fν
    defs = [
        ν.Ω0 => (ν.Neff/3) * 7/8 * (4/11)^(4/3) * γ.Ω0
        h.T0 => (ν.Neff/3)^(1/4) * (4/11)^(1/3) * γ.T0 # same as for massless neutrinos # TODO: are the massive neutrino density parameters correct?
        h.Ω0_massless => 7/8 * (h.T0/γ.T0)^4 * γ.Ω0 # Ω0 for corresponding massless neutrinos # TODO: reconcile with class? https://github.com/lesgourg/class_public/blob/ae99bcea1cd94994228acdfaec70fa8628ae24c5/source/background.c#L1561
        k => NaN # make background shut up # TODO: avoid
        fν => ν.ρ0 / (ν.ρ0 + ν.ρ0)
        C => 0.48 # TODO: why does ≈ 0.48 give better agreement with CLASS? # TODO: phi set here? https://github.com/lesgourg/class_public/blob/ae99bcea1cd94994228acdfaec70fa8628ae24c5/source/perturbations.c#L5713
        g.Ψ => 20C / (15 + 4fν) # Φ found from solving initialization system # TODO: is this correct when having both massless and massive neutrinos?
    ]
    eqs0 = [
        G.ρ ~ sum(s.ρ for s in species)
        b.rec.ρb ~ b.ρ * g.H0^2/GN # kg/m³ (convert from H0=1 units to SI units)
        b.rec.Tγ ~ γ.T
    ] .|> O(ϵ^0)
    eqs1 = [
        G.δρ ~ sum(s.δ * s.ρ for s in species) # total energy density perturbation
        G.Π ~ sum((s.ρ + s.P) * s.σ for s in species)
        b.θinteraction ~ #=g.k^2*csb²*bar.δ +=# -b.rec.dτ * 4*γ.ρ/(3*b.ρ) * (γ.θ - b.θ) # TODO: enable csb² when it seems stable... # TODO: define some common interaction type, e.g. momentum transfer # TODO: would love to write something like interaction = thompson_scattering(γ, b)
        γ.τ̇ ~ b.rec.dτ
        γ.θb ~ b.θ
    ] .|> O(ϵ^1)
    # TODO: do various IC types (adiabatic, isocurvature, ...) from here?
    connections = ODESystem([eqs0; eqs1], t, [], [pars; k]; defaults=defs, name)
    M = compose(connections, g, G, species...) |> complete
    return CosmologyModel(M; kwargs...)
end