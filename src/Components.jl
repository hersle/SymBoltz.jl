ϵ = only(GlobalScope.(@parameters ϵ)) # perturbative expansion parameter

"""
    metric(; name = :g, kwargs...)

Create a symbolic component for the perturbed FLRW spacetime metric in the conformal Newtonian gauge with signature

    a^2 * diag(-1-2*ϵ*Ψ, +1 - 2*ϵ*Φ, +1 - 2*ϵ*Φ, +1 - 2*ϵ*Φ)
"""
function metric(; name = :g, kwargs...)
    vars = a, ℰ, E, H, ℋ, Φ, Ψ, Φ′, Ψ′, g11, g22, h11, h22, z = GlobalScope.(@variables a(t) ℰ(t) E(t) H(t) ℋ(t) Φ(t) Ψ(t) Φ′(t) Ψ′(t) g11(t) g22(t) h11(t) h22(t) z(t))
    pars = H0, h = GlobalScope.(@parameters H0 h)
    defs = [
        H0 => H100 * h
        h => H0 / H100
    ]
    return ODESystem([
        z ~ 1/a - 1
        ℰ ~ D(a) / a # ℰ = ℋ/ℋ0 = ℋ/H0
        E ~ ℰ / a # E = H/H0
        ℋ ~ ℰ * H0
        H ~ E * H0
        g11 ~ -a^2
        g22 ~ +a^2
        h11 * ϵ ~ -2 * a^2 * Ψ * ϵ
        h22 * ϵ ~ -2 * a^2 * Φ * ϵ
        Φ′ * ϵ ~ D(Φ) * ϵ # TODO: why unstable at early times?
        Ψ′ * ϵ ~ D(Ψ) * ϵ # TODO: why unstable at early times?
    ], t, vars, pars; defaults = defs, name, kwargs...)
end

"""
    general_relativity(g; name = :G, kwargs...)

Create a symbolic component for the general relativistic (GR) theory of gravity in the spacetime with the metric `g`.
"""
function general_relativity(g; acceleration = false, name = :G, kwargs...)
    @variables ρ(t) P(t) ρcrit(t) δρ(t) Π(t) Δ(t)
    @parameters Δfac = 0
    a = g.a
    F1 = D(a)^2 ~ 8*Num(π)/3 * ρ * a^4 # Friedmann constraint equation
    F2 = D(D(a)) ~ D(a)^2/(2*a) * (1 - 3*P/ρ) # Friedmann acceleration equation (alternatively D(D(a)) ~ D(a)^2/a - 4*Num(π)/3 * (ρ + 3*P) * a^3)
    if acceleration
        eqs0 = [
            F2.lhs ~ F2.rhs + Δfac*Δ # use constraint damping # TODO: incorporate sign(∂Δ/∂a′), but positive in GR?
            Δ ~ F1.lhs - F1.rhs # violation of Friedmann constraint
        ] .|> O(ϵ^0)
        ics0 = [F1] .|> O(ϵ^0)
    else
        eqs0 = [
            F1.lhs^(1/2) ~ +√(F1.rhs) # "normal" Friedmann equation
            Δ ~ F2.lhs - F2.rhs # violation of acceleration equation
        ] .|> O(ϵ^0)
        ics0 = []
    end
    eqs1 = [
        D(g.Φ) ~ -4*Num(π)/3*a^2/g.ℰ*δρ - k^2/(3*g.ℰ)*g.Φ - g.ℰ*g.Ψ
        k^2 * (g.Φ - g.Ψ) ~ 12*Num(π) * a^2 * Π
    ] .|> O(ϵ^1)
    guesses = [ρ => 1, D(a) => +1]
    return ODESystem([eqs0; eqs1], t; initialization_eqs = ics0, guesses, name, kwargs...)
end

# TODO: potential
# TODO: mass
# TODO: optionally solve F1 directly for D(a) > 0 with quadratic formula
"""
    brans_dicke(g; name = :G, kwargs...)

"""
function brans_dicke(g; name = :G, acceleration = false, kwargs...)
    pars = @parameters ω
    vars = @variables ρ(t) P(t) ϕ(t) δρ(t) Π(t) Δ(t) G(t)
    a = g.a
    F1 = D(a)^2 ~ 8*Num(π)/3*ρ*a^4/ϕ - D(a)*a*D(ϕ)/ϕ + ω/6*a^2*(D(ϕ)/ϕ)^2
    F2 = D(D(a)) ~ D(a)^2/(2*a) - 4*Num(π)*a^3*P/ϕ - ω/4*a*(D(ϕ)/ϕ)^2 - D(a)/2*D(ϕ)/ϕ - a/2*D(D(ϕ))/ϕ
    KG = D(D(ϕ)) ~ 8*Num(π)/(2*ω+3) * a^2 * (ρ-3*P) - 2*D(a)/a*D(ϕ)
    eqs0 = [
        KG
        G ~ (2*ω+4) / (2*ω+3) / ϕ # effective gravitational constant
    ] .|> O(ϵ^0)
    ics0 = []
    if acceleration
        append!(eqs0, [
            F2
            Δ ~ F1.lhs - F1.rhs # violation of Friedmann contraint
        ] .|> O(ϵ^0))
        push!(ics0, F1 |> O(ϵ^0))
    else
        append!(eqs0, [
            D(a) ~ -a/2*D(ϕ)/ϕ + √((a/2*D(ϕ)/ϕ)^2 + 8*Num(π)/3*ρ*a^4/ϕ + ω/6*(a*D(ϕ)/ϕ)^2) # TODO: Symbolics.symbolic_solve F1 for ȧ
            Δ ~ F2.lhs - F2.rhs # violation of acceleration equation
        ] .|> O(ϵ^0))
    end
    # TODO: perturbations
    eqs1 = [
        D(g.Φ) ~ δρ
        g.Ψ ~ g.Φ + Π
    ] .|> O(ϵ^1)
    guesses = [ρ => 1.0, D(g.a) => +1.0]
    return ODESystem([eqs0; eqs1], t, vars, pars; name, initialization_eqs = ics0, guesses, kwargs...) # TODO: don't pass vars and pars
end

"""
    species_constant_eos(g, w, cs² = w, ẇ = 0, _σ = 0; θinteract = false, kwargs...)

Create a symbolic component for a particle species with equation of state `w ~ P/ρ` in the spacetime with the metric `g`.
"""
function species_constant_eos(g, w, ẇ = 0, _σ = 0; analytical = true, θinteract = false, kwargs...)
    @assert ẇ == 0 && _σ == 0 # TODO: relax (need to include in ICs)
    pars = analytical ? (@parameters ρ0 Ω0) : []
    vars = @variables ρ(t) P(t) Ω(t) δ(t) θ(t) Δ(t) θinteraction(t) σ(t) cs²(t) u(t) u′(t)
    eqs0 = [
        P ~ w * ρ # equation of state
        analytical ? (ρ ~ ρ0 * g.a^(-3*(1+w))) : (D(ρ) ~ -3 * g.ℰ * (ρ + P)) # alternative derivative: D(ρ) ~ -3 * g.ℰ * (ρ + P)
        Ω ~ 8*Num(π)/3 * ρ
    ] .|> O(ϵ^0)
    eqs1 = [
        D(δ) ~ -(1+w)*(θ-3*D(g.Φ)) - 3*g.ℰ*(cs²-w)*δ # Bertschinger & Ma (30) with Φ -> -Φ; or Baumann (4.4.173) with Φ -> -Φ
        D(θ) ~ -g.ℰ*(1-3*w)*θ - ẇ/(1+w)*θ + cs²/(1+w)*k^2*δ - k^2*σ + k^2*g.Ψ + θinteraction # Bertschinger & Ma (30) with θ = kv
        u ~ θ / k # velocity
        u′ ~ D(u)
        Δ ~ δ + 3(1+w) * g.ℰ/θ # Baumann (4.2.144) with v -> -u
        σ ~ _σ
    ] .|> O(ϵ^1)
    ics1 = [
        δ ~ -3/2 * (1+w) * g.Ψ # adiabatic: δᵢ/(1+wᵢ) == δⱼ/(1+wⱼ) (https://cmb.wintherscoming.no/theory_initial.php#adiabatic)
        θ ~ 1/2 * (k^2*t) * g.Ψ # # TODO: include σ ≠ 0 # solve u′ + ℋ(1-3w)u = w/(1+w)*kδ + kΨ with Ψ=const, IC for δ, Φ=-Ψ, ℋ=H₀√(Ωᵣ₀)/a after converting ′ -> d/da by gathering terms with u′ and u in one derivative using the trick to multiply by exp(X(a)) such that X′(a) will "match" the terms in front of u
    ] .|> O(ϵ^1)
    defs = analytical ? [ρ0 => 3/(8*Num(π)) * Ω0] : Dict()
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
function photons(g; polarization = true, lmax = 6, name = :γ, kwargs...)
    lmax >= 3 || error("Need lmax >= 3")
    γ = radiation(g; name, kwargs...) |> thermodynamics |> complete # prevent namespacing in extension below

    vars = @variables F0(t) F(t)[1:lmax] δ(t) θ(t) σ(t) τ̇(t) θb(t) Π(t) G0(t) G(t)[1:lmax]
    defs = [
        γ.Ω0 => π^2/15 * (kB*γ.T0)^4 / (ħ^3*c^5) * 8π*GN / (3*g.H0^2)
    ]
    eqs1 = [
        D(F0) ~ -k*F[1] + 4*D(g.Φ)
        D(F[1]) ~ k/3*(F0-2*F[2]+4*g.Ψ) - 4/3 * τ̇/k * (θb - θ)
        D(F[2]) ~ 2/5*k*F[1] - 3/5*k*F[3] + 9/10*τ̇*F[2] - 1/10*τ̇*(G0+G[2])
        [D(F[l]) ~ k/(2*l+1) * (l*F[l-1] - (l+1)*F[l+1]) + τ̇*F[l] for l in 3:lmax-1]... # TODO: Π in last term here?
        D(F[lmax]) ~ k*F[lmax-1] - (lmax+1) / t * F[lmax] + τ̇ * F[lmax]
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
        append!(eqs1, [
            D(G0) ~ k * (-G[1]) - τ̇ * (-G0 + Π/2)
            D(G[1]) ~ k/3 * (G0 - 2*G[2]) - τ̇ * (-G[1])
            [D(G[l]) ~ k/(2*l+1) * (l*G[l-1] - (l+1)*G[l+1]) - τ̇ * (-G[l] + Π/10*δkron(l,2)) for l in 2:lmax-1]... # TODO: collect all equations here once G[0] works
            D(G[lmax]) ~ k*G[lmax-1] - (lmax+1) / t * G[lmax] + τ̇ * G[lmax]
        ] .|> O(ϵ^1))
        append!(ics1, [
            G0 ~ 0 #5/4 * Θ[2],
            G[1] ~ 0 #-1/4 * k/dτ * Θ[2],
            G[2] ~ 0 #1/4 * Θ[2],
            [G[l] ~ 0 #=-l/(2*l+1) * k/dτ * ΘP[l-1]=# for l in 3:lmax]...    
        ] .|> O(ϵ^1))
    else
        append!(eqs1, [G0 ~ 0, collect(G .~ 0)...] .|> O(ϵ^1)) # pin to zero
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
    x, W = gauss(x -> x^2 * f0(x), nx, 0.0, 1e3) # Gaussian quadrature weights, reduced momentum bins x = q*c / (kB*T0) # these points give accurate integral for Iρmν in the background, at least # TODO: ok for perturbations?
    ∫dx_x²_f0(f) = sum(collect(f) .* W) # a function that approximates the weighted integral ∫dx*x^2*f(x)*f0(x)

    E(x, y) = √(x^2 + y^2)
    Iρ(y) = ∫dx_x²_f0(@. E(x, y)) # Iρ(0) = 7π^4/120
    IP(y) = ∫dx_x²_f0(@. x^2 / E(x, y)) # IP(0) = Iρ(0)
    
    eqs0 = [
        T ~ T0 / g.a
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
        ρ0 => 3/(8*Num(π)) * Ω0
        ρ0_massless => 3/8π * Ω0_massless
        m => 0.02 * eV/c^2 # one massive neutrino with this mass # TODO: specify by user
        y0 => m*c^2 / (kB*T0)
    ]
    ics1 = []
    for i in 1:nx
        append!(eqs1, [
            D(ψ0[i]) ~ -k * x[i]/E(x[i],y) * ψ[i,1] - D(g.Φ) * dlnf0_dlnx(x[i])
            D(ψ[i,1]) ~ k/3 * x[i]/E(x[i],y) * (ψ0[i] - 2*ψ[i,2]) - k/3 * E(x[i],y)/x[i] * g.Ψ * dlnf0_dlnx(x[i])
            [D(ψ[i,l]) ~ k/(2*l+1) * x[i]/E(x[i],y) * (l*ψ[i,l-1] - (l+1)*ψ[i,l+1]) for l in 2:lmax]...
            ψ[i,lmax+1] ~ (2*lmax+1) * E(x[i],y)/x[i] * ψ[i,lmax] / (k*t) - ψ[i,lmax-1]
        ] .|> O(ϵ^1))
        append!(ics1, [
            ψ0[i] ~ -1/4 * (-2*g.Ψ) * dlnf0_dlnx(x[i])
            ψ[i,1] ~ -1/(3*k) * E(x[i],y)/x[i] * (1/2*(k^2*t)*g.Ψ) * dlnf0_dlnx(x[i])
            ψ[i,2] ~ -1/2 * (1/15*(k*t)^2*g.Ψ) * dlnf0_dlnx(x[i])
            [ψ[i,l] ~ 0 for l in 3:lmax] # TODO: full ICs
        ] .|> O(ϵ^1))
    end
    return ODESystem([eqs0; eqs1], t, vars, pars; initialization_eqs=ics1, defaults=defs, name, kwargs...)
end

"""
    cold_dark_matter(g; name = :c, kwargs...)

Create a particle species for cold dark matter in the spacetime with metric `g`.
"""
function cold_dark_matter(g; name = :c, kwargs...)
    c = matter(g; name, kwargs...) |> complete
    c = extend(c, ODESystem([c.cs² ~ 0] .|> O(ϵ^1), t, [], []; name))
    return c
end

"""
    baryons(g; recombination=true, name = :b, kwargs...)

Create a particle species for baryons in the spacetime with metric `g`.
"""
function baryons(g; recombination=true, name = :b, kwargs...)
    b = matter(g; θinteract=true, name, kwargs...) |> complete
    if recombination # TODO: dont add recombination system when recombination = false
        @named rec = thermodynamics_recombination_recfast(g)
    else
        vars = @variables τ(t) τ̇(t) ρb(t) Tγ(t) cs²(t)
        @named rec = ODESystem([τ ~ 0, τ̇ ~ 0, cs² ~ 0], t, vars, [])
    end
    b = extend(b, ODESystem([b.cs² ~ rec.cs²] .|> O(ϵ^1), t, [], []; name))
    b = compose(b, rec)
    return b
end

"""
    quintessence(g, v; name = :ϕ, kwargs...)

Create a species for a quintessence scalar field with potential `v` in the spacetime with metric `g`.
"""
function quintessence(g, v; name = :Q, kwargs...)
    @variables ϕ(t) ρ(t) P(t) w(t) δ(t) σ(t) V(t) V′(t) V″(t) K(t) m²(t) ϵs(t) ηs(t)
    ∂_∂ϕ = Differential(ϕ)
    eqs0 = [
        V ~ v(ϕ)
        V′ ~ ∂_∂ϕ(v(ϕ)) |> expand_derivatives |> simplify
        V″ ~ ∂_∂ϕ(∂_∂ϕ(v(ϕ))) |> expand_derivatives |> simplify
        K ~ (D(ϕ)/g.a)^2 / 2 # ϕ̇²/2 = (ϕ′/a)²/2
        D(D(ϕ)) ~ -2 * g.ℰ * D(ϕ) - g.a^2 * V′ # with cosmic time: ϕ̈ + 3*E*ϕ̇ + V′ = 0
        ρ ~ K + V
        P ~ K - V
        w ~ P / ρ
        m² ~ V″
        ϵs ~ (V′/V)^2 / (16*Num(π)) # 1st slow roll parameter
        ηs ~ (V″/V) / (8*Num(π)) # 2nd slow roll parameter
    ] .|> O(ϵ^0)
    eqs1 = [ # TODO: perturbations
        δ ~ 0
        σ ~ 0
    ] .|> O(ϵ^1)
    return ODESystem([eqs0; eqs1], t; name, kwargs...)
end

function background(sys)
    sys = thermodynamics(sys)
    return replace(sys, sys.b.rec => ODESystem([], t; name = :rec))
end

function thermodynamics(sys)
    return transform((sys, _) -> extract_order(sys, [0]), sys)
end

function perturbations(sys; spline_thermo = true)
    if spline_thermo && !all([eq.rhs === 0 for eq in equations(sys.b.rec)])
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
        c = cold_dark_matter(g; name = :c),
        b = baryons(g; recombination, name = :b),
        Λ = cosmological_constant(g),
        kwargs...
    )

Create a ΛCDM model.
"""
function ΛCDM(;
    recombination = true,
    acceleration = false,
    Λanalytical = false,
    g = metric(),
    G = general_relativity(g; acceleration),
    γ = photons(g),
    ν = massless_neutrinos(g),
    h = massive_neutrinos(g),
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
        G.Π ~ sum((s.ρ + s.P) * s.σ for s in species)
        b.θinteraction ~ k^2*b.cs²*b.δ + -b.rec.τ̇ * 4*γ.ρ/(3*b.ρ) * (γ.θ - b.θ) # TODO: define some common interaction type, e.g. momentum transfer # TODO: would love to write something like interaction = thompson_scattering(γ, b)
        γ.τ̇ ~ b.rec.τ̇
        γ.θb ~ b.θ
    ] .|> O(ϵ^1)
    # TODO: do various IC types (adiabatic, isocurvature, ...) from here?
    connections = ODESystem([eqs0; eqs1], t, [], [pars; k]; defaults=defs, name)
    M = compose(connections, g, G, species...)
    initE = !all([:Ω0 in Symbol.(parameters(s)) for s in species])
    if !initE # add default if all species are analytical
        defs = merge(Dict(s.Ω0 => 1 - (sum(s′.Ω0 for s′ in species if s′ != s)) for s in species), defaults(M)) # TODO: generalize to non-GR
        M = extend(M, ODESystem([], t; defaults = defs, name))
    end
    return CosmologyModel(M; initE, kwargs...)
end

function parameters_Planck18(M::CosmologyModel)
    # https://arxiv.org/pdf/1807.06209, Table 5
    h = 0.6736
    return [
        M.g.h => h
        M.γ.T0 => 2.7255
        M.c.Ω0 => 0.1200 / h^2
        M.b.Ω0 => 0.02237 / h^2
        M.b.rec.Yp => 0.2454

        # TODO: how to handle
        # 1) backwards with default for Λ.Ω0? --> skip E = 1
        # 2) forwards with no default?        --> shooting
        #M.Λ.Ω => 1 - sum(s.Ω0 for s in [M.γ, M.ν, M.h, M.c, M.b])
    ]
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
