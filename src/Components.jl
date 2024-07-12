ϵ = only(GlobalScope.(@parameters ϵ)) # perturbative expansion parameter

function metric(; kwargs...)
    a, ℰ, E, H, ℋ, Φ, Ψ = GlobalScope.(@variables a(t) ℰ(t) E(t) H(t) ℋ(t) Φ(t) Ψ(t)) # TODO: more natural way to connect them?
    H0, h = GlobalScope.(@parameters H0 h)
    return ODESystem([
        ℰ ~ D(a) / a # ℰ = ℋ/ℋ0
        E ~ ℰ / a # E = H/H0
        ℋ ~ ℰ * H0
        H ~ E * H0
    ], t, [a, ℰ, E, H, ℋ, Φ, Ψ], [H0, h]; defaults = [H0 => H100 * h], kwargs...)
end

function gravity(g; kwargs...)
    @variables ρ(t) ρcrit(t) δρ(t) Π(t)
    eqs0 = [
        D(g.a) ~ √(8π/3 * ρ) * g.a^2 # Friedmann equation
        ρcrit ~ 3/8π * g.E^2 # critical density (H² = 8πG/3 * ρcrit)
    ] .|> O(ϵ^0)
    eqs1 = [
        D(g.Φ) ~ -4π/3*g.a^2/g.ℰ*δρ - k^2/(3*g.ℰ)*g.Φ - g.ℰ*g.Ψ
        k^2 * (g.Φ - g.Ψ) ~ 12π * g.a^2 * Π
    ] .|> O(ϵ^1)
    return ODESystem([eqs0; eqs1], t, [ρ, ρcrit, δρ, Π], []; kwargs...)
end

function species_constant_eos(g, w, cs² = w, ẇ = 0, _σ = 0; θinteract = false, kwargs...)
    @assert ẇ == 0 && _σ == 0 # TODO: relax (need to include in ICs)
    pars = @parameters Ω0 ρ0 # TODO: split pars0, pars1
    vars = @variables ρ(t) P(t) δ(t) θ(t) Δ(t) θinteraction(t) σ(t) # TODO: split vars0, vars1
    eqs0 = [
        P ~ w * ρ # equation of state
        ρ ~ ρ0 * g.a^(-3*(1+w)) # alternative derivative: D(ρ) ~ -3 * g.ℰ * (ρ + P)
    ] .|> O(ϵ^0)
    eqs1 = [
        D(δ) ~ -(1+w)*(θ-3*D(g.Φ)) - 3*g.ℰ*(cs²-w)*δ # Bertschinger & Ma (30) with Φ -> -Φ; or Baumann (4.4.173) with Φ -> -Φ
        D(θ) ~ -g.ℰ*(1-3*w)*θ - ẇ/(1+w)*θ + cs²/(1+w)*k^2*δ - k^2*σ + k^2*g.Ψ + θinteraction # Bertschinger & Ma (30) with θ = kv
        Δ ~ δ + 3(1+w) * g.ℰ/θ # Baumann (4.2.144) with v -> -u
        σ ~ _σ
    ] .|> O(ϵ^1)
    ics1 = [
        δ ~ -3/2 * (1+w) * g.Ψ # adiabatic: δᵢ/(1+wᵢ) == δⱼ/(1+wⱼ) (https://cmb.wintherscoming.no/theory_initial.php#adiabatic)
        θ ~ 1/2 * (k^2*t) * g.Ψ # # TODO: fix with θ: -1/2 * k*t * g.Ψ # TODO: include σ ≠ 0 # solve u′ + ℋ(1-3w)u = w/(1+w)*kδ + kΨ with Ψ=const, IC for δ, Φ=-Ψ, ℋ=H₀√(Ωᵣ₀)/a after converting ′ -> d/da by gathering terms with u′ and u in one derivative using the trick to multiply by exp(X(a)) such that X′(a) will "match" the terms in front of u
    ] .|> O(ϵ^1)
    defs = [
        ρ0 => 3/8π*Ω0
    ]
    !θinteract && push!(eqs1, (θinteraction ~ 0) |> O(ϵ^1))
    return ODESystem([eqs0; eqs1], t, vars, pars; initialization_eqs=ics1, defaults=defs, kwargs...)
end

function matter(g; kwargs...)
    return species_constant_eos(g, 0; kwargs...)
end

function radiation(g; kwargs...)
    return species_constant_eos(g, 1//3; kwargs...)
end

function cosmological_constant(g; kwargs...)
    Λ = species_constant_eos(g, -1; kwargs...)
    Λ = background(Λ) # discard nonexisting perturbations
    @variables δ(t) σ(t)
    @parameters ρ0 Ω0 # TODO: avoid re-creating
    defs = [ρ0 => 3/8π*Ω0]
    Λ = extend(ODESystem([δ ~ 0, σ ~ 0] .|> O(ϵ^1), t, [δ, σ], [ρ0, Ω0]; defaults=defs, name=:Λ), complete(Λ)) # no perturbations
    return Λ
end

function photons(g; polarization=true, lmax=6, kwargs...)
    @parameters T0
    @variables T(t) F0(t) F(t)[1:lmax] δ(t) θ(t) σ(t) τ̇(t) θb(t) Π(t) G0(t) G(t)[1:lmax]

    γ = background(radiation(g; kwargs...))
    γ = complete(γ) # prevent namespacing in extension below
    γ = extend(γ, ODESystem([
        T ~ T0 / g.a # alternative derivative: D(Tγ) ~ -1*Tγ * g.ℰ
    ], t, [T], [T0]; defaults = [
        T0 => (15/π^2 * γ.ρ0 * g.H0^2/GN * ħ^3*c^5)^(1/4) / kB
    ], name=:γ))

    # perturbations
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
        (polarization ? [
            D(G0) ~ k * (-G[1]) - τ̇ * (-G0 + Π/2)
            D(G[1]) ~ k/3 * (G0 - 2*G[2]) - τ̇ * (-G[1])
            [D(G[l]) ~ k/(2*l+1) * (l*G[l-1] - (l+1)*G[l+1]) - τ̇ * (-G[l] + Π/10*δkron(l,2)) for l in 2:lmax-1]... # TODO: collect all equations here once G[0] works
            D(G[lmax]) ~ k*G[lmax-1] - (lmax+1) / t * G[lmax] + τ̇ * G[lmax]
        ] : [
            G0 ~ 0, collect(G .~ 0)... # pin to zero
        ])...
    ] .|> O(ϵ^1)
    ics1 = [
        δ ~ -2 * g.Ψ # Dodelson (7.89)
        θ ~ 1/2 * (k^2*t) * g.Ψ # Dodelson (7.95)
        F[2] ~ 0 # (polarization ? -8/15 : -20/45) * k/dτ * Θ[1], # depends on whether polarization is included # TODO: move to initialization_eqs?
        [F[l] ~ 0 #=-l/(2*l+1) * k/dτ * Θ[l-1]=# for l in 3:lmax]...
        G0 ~ 0 #5/4 * Θ[2],
        G[1] ~ 0 #-1/4 * k/dτ * Θ[2],
        G[2] ~ 0 #1/4 * Θ[2],
        [G[l] ~ 0 #=-l/(2*l+1) * k/dτ * ΘP[l-1]=# for l in 3:lmax]...
    ] .|> O(ϵ^1)
    γ = extend(γ, ODESystem(eqs1, t, [γ.ρ, δ, θ, σ, τ̇, θb], [T0]; initialization_eqs=ics1, kwargs...))
    return γ
end

function massless_neutrinos(g; lmax=6, kwargs...)
    ν = radiation(g; kwargs...) |> background |> complete

    vars = @variables F0(t) F(t)[1:lmax+1] δ(t) θ(t) σ(t)
    pars = @parameters Neff
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
    defs = [Neff => 3.046]
    return extend(ν, ODESystem(eqs1, t, vars, pars; initialization_eqs=ics1, defaults=defs, kwargs...))
end

# TODO: use vector equations and simplify loops
function massive_neutrinos(g; nx=5, lmax=4, kwargs...)
    pars = @parameters Ω0_massless ρ0_massless Ω0 ρ0 m T0 y0
    vars = @variables ρ(t) T(t) y(t) P(t) w(t) δ(t) σ(t) θ(t) ψ0(t)[1:nx] ψ(t)[1:nx,1:lmax+1]

    f0(x) = 1 / (exp(x) + 1) # TODO: why not exp(E)?
    dlnf0_dlnx(x) = -x / (1 + exp(-x))
    E(x, y) = √(x^2 + y^2)
    x, W = gauss(x -> x^2 * f0(x), nx, 0.0, 1e3) # Gaussian quadrature weights, reduced momentum bins x = q*c / (kB*T0) # these points give accurate integral for Iρmν in the background, at least # TODO: ok for perturbations? # TODO: also include common x^2 factor in weighting?
    ∫dx_x²_f0(f) = sum(collect(f) .* W) # a function that approximates the weighted integral ∫dx*x^2*f(x)*f0(x)
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
        σ ~ (2/3) * ∫dx_x²_f0(@. x^2/E(x,y)*ψ[:,2]) / (∫dx_x²_f0(@. E(x,y)) + 1/3*∫dx_x²_f0(@. x^2/E(x,y)))
    ] .|> O(ϵ^1)
    defs = [
        Ω0 => Ω0_massless * Iρ(y0)/Iρ(0) # ≈ Ω0_massless * (3ζ(3)/2)/(7π^4/120) * y0 for y0 → ∞
        ρ0 => 3/8π * Ω0
        ρ0_massless => 3/8π * Ω0_massless
        m => 0.02 * eV/c^2 # one massive neutrino with this mass
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
            [ψ[i,l] ~ 0 for l in 3:lmax] # TODO: proper ICs    
        ] .|> O(ϵ^1))
    end
    return ODESystem([eqs0; eqs1], t, vars, pars; initialization_eqs=ics1, defaults=defs, kwargs...)
end

function baryons(g; recombination=true, kwargs...)
    b = matter(g; θinteract=true, kwargs...)
    if recombination
        @named rec = thermodynamics_recombination_recfast(g)
        b = compose(b, rec)
    end
    return b
end

function transform(f::Function, sys::ODESystem; fullname=string(sys.name))
    subs = [transform(f, sub; fullname = (ModelingToolkit.iscomplete(sys) ? "" : fullname * "₊") * string(sub.name)) for sub in sys.systems]
    sys = f(sys, fullname)
    return compose(sys, subs)
end

#=
function basename(sys::ODESystem)
    name = ModelingToolkit.get_name(sys) |> string
    prefix = ModelingToolkit.get_namespace(sys) * "₊"
    return chopprefix(name, prefix) |> Symbol
end
=#

function replace(sys::ODESystem, old_new_subsys::Pair{ODESystem, ODESystem})
    old_subsys, new_subsys = old_new_subsys # unpack
    fullname_target = ModelingToolkit.get_name(old_subsys) |> string
    return transform((sys, fullname) -> (fullname == fullname_target ? new_subsys : identity(sys)), sys)
end

# for testing: transform(identity, sys) should do no harm to a system
function identity(sys)
    eqs = ModelingToolkit.get_eqs(sys)
    ieqs = ModelingToolkit.get_initialization_eqs(sys)
    vars = ModelingToolkit.get_unknowns(sys)
    pars = ModelingToolkit.get_ps(sys)
    defs = ModelingToolkit.get_defaults(sys)
    guesses = ModelingToolkit.get_guesses(sys)
    return ODESystem(eqs, t, vars, pars; initialization_eqs=ieqs, defaults=defs, guesses=guesses, name=sys.name)
end

O(x, ϵⁿ) = x * ϵⁿ
O(eq::Equation, ϵⁿ) = O(eq.lhs, ϵⁿ) ~ O(eq.rhs, ϵⁿ)
O(ϵⁿ) = x -> O(x, ϵⁿ)

function extract_order(expr, order)
    if order == 0
        return substitute(expr, ϵ => 0)
    else
        expr = Differential(ϵ)(expr) |> expand_derivatives # differentiate away one power of ϵ^order -> order*ϵ^(order-1)
        expr = expr / order # remove prefactor from differentiation
        return extract_order(expr, order - 1)
    end
end

function extract_order(eq::Equation, order)
    return extract_order(eq.lhs, order) ~ extract_order(eq.rhs, order)
end

function extract_order(sys::ODESystem, orders)
    eqs = ModelingToolkit.get_eqs(sys)
    ieqs = ModelingToolkit.get_initialization_eqs(sys)
    vars = ModelingToolkit.get_unknowns(sys)
    pars = ModelingToolkit.get_ps(sys)
    defs = ModelingToolkit.get_defaults(sys)
    guesses = ModelingToolkit.get_guesses(sys)

    # extract requested orders
    eqs = vcat((extract_order.(eqs, order) for order in orders)...)
    ieqs = vcat((extract_order.(ieqs, order) for order in orders)...)

    # remove resulting trivial equations
    eqs = filter(eq -> eq != (0 ~ 0), eqs)
    ieqs = filter(eq -> eq != (0 ~ 0), ieqs)

    sys0 = ODESystem(eqs, t, vars, pars; initialization_eqs=ieqs, defaults=defs, guesses=guesses, name=sys.name)
    return sys0
end

function background(sys)
    return transform((sys, _) -> extract_order(sys, [0]), sys)
end

function perturbations(sys; spline_thermo=true)
    if spline_thermo
        @named rec = thermodynamics_recombination_splined()
        sys = replace(sys, sys.b.rec => rec)
    end
    return transform((sys, _) -> extract_order(sys, [0, 1]), sys)
end

function ΛCDM(; kwargs...)
    @named g = metric()
    @named G = gravity(g)
    @named γ = photons(g)
    @named ν = massless_neutrinos(g)
    @named h = massive_neutrinos(g)
    @named c = matter(g)
    @named b = baryons(g)
    @named Λ = cosmological_constant(g)
    species = [γ, ν, c, b, h, Λ]
    initialization_eqs = [
        g.a ~ √(γ.Ω0 + ν.Ω0 + h.Ω0_massless) * t # analytical radiation-dominated solution # TODO: write t ~ 1/g.ℰ ?
    ]
    @parameters C
    defaults = [
        species[end].Ω0 => 1 - sum(s.Ω0 for s in species[begin:end-1]) # TODO: solve nonlinear system # TODO: any combination of all but one species
        ν.Ω0 => (ν.Neff/3) * 7/8 * (4/11)^(4/3) * γ.Ω0
        h.T0 => (ν.Neff/3)^(1/4) * (4/11)^(1/3) * γ.T0 # same as for massless neutrinos # TODO: are the massive neutrino density parameters correct?
        h.Ω0_massless => 7/8 * (h.T0/γ.T0)^4 * γ.Ω0 # Ω0 for corresponding massless neutrinos # TODO: reconcile with class? https://github.com/lesgourg/class_public/blob/ae99bcea1cd94994228acdfaec70fa8628ae24c5/source/background.c#L1561
        k => NaN # make background shut up # TODO: avoid
        # # TODO: fν => bg.neu.ρ0 / (bg.neu.ρ0 + bg.ph.ρ0)
        C => 0.48 # TODO: why does ≈ 0.48 give better agreement with CLASS? # TODO: phi set here? https://github.com/lesgourg/class_public/blob/ae99bcea1cd94994228acdfaec70fa8628ae24c5/source/perturbations.c#L5713
        g.Ψ => 20C / (15 #=+ 4fν=#) # Φ found from solving initialization system # TODO: is this correct when having both massless and massive neutrinos?
    ]
    eqs0 = [
        G.ρ ~ sum(s.ρ for s in species)
        b.rec.ρb ~ b.ρ * g.H0^2/GN # kg/m³ (convert from H0=1 units to SI units)
        b.rec.Tγ ~ γ.T
    ] .|> O(ϵ^0)
    eqs1 = [
        G.δρ ~ sum(s.δ * s.ρ for s in species) # total energy density perturbation
        G.Π ~ sum((s.ρ + s.P) * s.σ for s in species)
        b.θinteraction ~ #=g.k^2*csb²*bar.δ +=# -b.rec.dτ * 4*γ.ρ/(3*b.ρ) * (γ.θ - b.θ) # TODO: enable csb² when it seems stable... # TODO: define some common interaction type, e.g. momentum transfer
        γ.τ̇ ~ b.rec.dτ
        γ.θb ~ b.θ
    ] .|> O(ϵ^1)
    connections = ODESystem([eqs0; eqs1], t, [], [C, k]; initialization_eqs, defaults, kwargs...)
    return compose(connections, g, G, species...)
end

#=
# TODO: would love something like this to work
function ΛCDM()
    g = metric()
    G = gravity_GR()
    γ = photons()
    b = baryons()
    c = cold_dark_matter()
    h = massive_neutrinos()
    Λ = cosmological_constant()
    interaction = thompson_scattering(γ, b) # TODO: how???
    species = [γ, b, c, h, Λ]
    return cosmology(g, G, species)
end
=#