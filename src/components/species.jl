"""
    species_constant_eos(g, w, cₛ² = w, ẇ = 0, _σ = 0; θinteract = false, kwargs...)

Create a symbolic component for a particle species with equation of state `w ~ P/ρ` in the spacetime with the metric `g`.
"""
function species_constant_eos(g, _w, ẇ = 0, _σ = 0; analytical = true, θinteract = false, adiabatic = false, name = :s, kwargs...)
    @assert ẇ == 0 && _σ == 0 # TODO: relax (need to include in ICs)
    pars = analytical ? (@parameters ρ₀ Ω₀) : []
    vars = @variables w(t) ρ(t) P(t) Ω(t) δ(t) θ(t) Δ(t) θinteraction(t) σ(t) cₛ²(t) u(t) u̇(t)
    eqs0 = [
        w ~ _w # equation of state
        P ~ w * ρ # equation of state
        analytical ? (ρ ~ ρ₀ * g.a^(-3*(1+w))) : (D(ρ) ~ -3 * g.ℰ * (ρ + P)) # alternative derivative: D(ρ) ~ -3 * g.ℰ * (ρ + P)
        Ω ~ 8*Num(π)/3 * ρ
    ] .|> O(ϵ^0)
    eqs1 = [
        D(δ) ~ -(1+w)*(θ-3*D(g.Φ)) - 3*g.ℰ*(cₛ²-w)*δ # Bertschinger & Ma (30) with Φ -> -Φ; or Baumann (4.4.173) with Φ -> -Φ
        D(θ) ~ -g.ℰ*(1-3*w)*θ - ẇ/(1+w)*θ + cₛ²/(1+w)*k^2*δ - k^2*σ + k^2*g.Ψ + θinteraction # Bertschinger & Ma (30) with θ = kv
        u ~ θ / k # velocity
        u̇ ~ D(u)
        Δ ~ δ + 3(1+w) * g.ℰ/θ # Baumann (4.2.144) with v -> -u
        σ ~ _σ
    ] .|> O(ϵ^1)
    adiabatic && push!(eqs1, O(ϵ^1)(cₛ² ~ w))
    ics1 = [
        δ ~ -3//2 * (1+w) * g.Ψ # adiabatic: δᵢ/(1+wᵢ) == δⱼ/(1+wⱼ) (https://cmb.wintherscoming.no/theory_initial.php#adiabatic) # TODO: match CLASS with higher-order (for photons)? https://github.com/lesgourg/class_public/blob/22b49c0af22458a1d8fdf0dd85b5f0840202551b/source/perturbations.c#L5631-L5632
        θ ~ 1//2 * (k^2*t) * g.Ψ # # TODO: include σ ≠ 0 # solve u′ + ℋ(1-3w)u = w/(1+w)*kδ + kΨ with Ψ=const, IC for δ, Φ=-Ψ, ℋ=H₀√(Ωᵣ₀)/a after converting ′ -> d/da by gathering terms with u′ and u in one derivative using the trick to multiply by exp(X(a)) such that X′(a) will "match" the terms in front of u
    ] .|> O(ϵ^1)
    defs = analytical ? [ρ₀ => 3/(8*Num(π)) * Ω₀] : Dict()
    !θinteract && push!(eqs1, (θinteraction ~ 0) |> O(ϵ^1))
    return ODESystem([eqs0; eqs1], t, vars, pars; initialization_eqs=ics1, defaults=defs, name, kwargs...)
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
    r = species_constant_eos(g, 1//3; name, kwargs...) |> complete
    pars = @parameters T₀
    vars = @variables T(t) # TODO: define in constant_eos? https://physics.stackexchange.com/questions/650508/whats-the-relation-between-temperature-and-scale-factor-for-arbitrary-eos-1
    return extend(r, ODESystem([T ~ T₀ / g.a], t, vars, pars; name))
end

"""
    cosmological_constant(g; name = :Λ, kwargs...)

Create a particle species for the cosmological constant (with equation of state `w ~ -1`) in the spacetime with metric `g`.
"""
function cosmological_constant(g; name = :Λ, analytical = false, kwargs...)
    description = "Cosmological constant"
    Λ = species_constant_eos(g, -1; name, analytical, description, kwargs...) |> thermodynamics |> complete # discard ill-defined perturbations
    vars = @variables δ(t) θ(t) σ(t)
    description = "Cosmological constant"
    return extend(Λ, ODESystem([δ ~ 0, θ ~ 0, σ ~ 0, Λ.cₛ² ~ -1] .|> O(ϵ^1), t, vars, []; name); description) # manually set perturbations to zero
end

"""
    w0wa(g; kwargs...)

Create a particle species for the w₀-wₐ dark energy (CPL) parametrization in the spacetime with metric `g`.
"""
function w0wa(g; name = :X, analytical = false, kwargs...)
    pars = @parameters w0 wa cₛ²
    vars = @variables ρ(t) P(t) w(t) ẇ(t) δ(t) θ(t) σ(t)
    # TODO: generate equations with a generic species_eos function
    eqs0 = [
        w ~ w0 + wa * (1 - g.a) # equation of state
        ẇ ~ D(w)
        P ~ w * ρ # pressure
    ] .|> SymBoltz.O(ϵ^0) # O(ϵ⁰) multiplies all equations by 1 (no effect, but see step 5)
    if analytical
        append!(pars, @parameters ρ₀ Ω₀)
        push!(eqs0, ρ ~ ρ₀ * abs(g.a)^(-3 * (1 + w0 + wa)) * exp(-3 * wa * (1 - g.a))) # energy density # TODO: get rid of abs
        defaults = Dict(ρ₀ => 3/(8*Num(π)) * Ω₀)
    else
        push!(eqs0, D(ρ) ~ -3 * g.ℰ * ρ * (1 + w))
        defaults = Dict()
    end
    eqs1 = [
        D(δ) ~ -(1 + w) * (θ - 3*D(g.Φ)) - 3 * g.ℰ * (cₛ² - w) * δ # energy overdensity
        D(θ) ~ -g.ℰ * (1 - 3*w) * θ - D(w) / (1 + w) * θ + cₛ² / (1 + w) * k^2 * δ - k^2 * σ + k^2 * g.Ψ # momentum
        σ ~ 0 # shear stress
    ] .|> SymBoltz.O(ϵ^1) # O(ϵ¹) multiplies all equations by ϵ, marking them as perturbation equations
    ics1 = [
        δ ~ -3//2 * (1+w) * g.Ψ
        θ ~ 1//2 * (k^2*t) * g.Ψ
    ] .|> SymBoltz.O(ϵ^1)
    description = "w₀wₐ (CPL) dark energy"
    return ODESystem([eqs0; eqs1], t, vars, pars; initialization_eqs=ics1, defaults, name, description, kwargs...)
end

"""
    photons(g; polarization = true, lmax = 6, name = :γ, kwargs...)

Create a particle species for photons in the spacetime with metric `g`.
"""
function photons(g; polarization = true, lmax = 6, name = :γ, kwargs...)
    lmax >= 3 || error("Need lmax >= 3")
    description = "Photons"
    γ = radiation(g; name, description, kwargs...) |> thermodynamics |> complete # prevent namespacing in extension below

    vars = @variables F(t)[0:lmax] δ(t) θ(t) σ(t) τ̇(t) θb(t) Π(t) Π̇(t) G(t)[0:lmax]
    defs = [
        γ.Ω₀ => π^2/15 * (kB*γ.T₀)^4 / (ħ^3*c^5) * 8π*GN / (3*g.H₀^2)
    ]
    eqs1 = [
        # Bertschinger & Ma (64) with anₑσₜ -> -τ̇
        D(F[0]) ~ -k*F[1] + 4*D(g.Φ)
        D(F[1]) ~ k/3*(F[0]-2*F[2]+4*g.Ψ) - 4//3 * τ̇/k * (θb - θ)
        D(F[2]) ~ k/5*(2*F[1] - 3*F[3]) + 9//10*τ̇*F[2] - 1//10*τ̇*(G[0]+G[2])
        [D(F[l]) ~ k/(2*l+1) * (l*F[l-1] - (l+1)*F[l+1]) + τ̇*F[l] for l in 3:lmax-1]...
        D(F[lmax]) ~ k*F[lmax-1] - (lmax+1) / t * F[lmax] + τ̇ * F[lmax]
        δ ~ F[0]
        θ ~ 3*k*F[1]/4
        σ ~ F[2]/2
        Π ~ F[2] + G[0] + G[2]
        Π̇ ~ D(Π)
        γ.cₛ² ~ 1//3
    ] .|> O(ϵ^1)
    ics1 = [
        δ ~ -2 * g.Ψ # Dodelson (7.89)
        θ ~ 1//2 * (k^2*t) * g.Ψ # Dodelson (7.95)
        F[2] ~ 0 # (polarization ? -8/15 : -20/45) * k/dτ * Θ[1], # depends on whether polarization is included
        [F[l] ~ 0 #=-l/(2*l+1) * k/dτ * Θ[l-1]=# for l in 3:lmax]...
    ] .|> O(ϵ^1)
    if polarization
        append!(eqs1, [
            D(G[0]) ~ k * (-G[1]) - τ̇ * (-G[0] + Π/2)
            [D(G[l]) ~ k/(2*l+1) * (l*G[l-1] - (l+1)*G[l+1]) - τ̇ * (-G[l] + Π/10*δkron(l,2)) for l in 1:lmax-1]...
            D(G[lmax]) ~ k*G[lmax-1] - (lmax+1) / t * G[lmax] + τ̇ * G[lmax]
        ] .|> O(ϵ^1))
        append!(ics1, [
            G[0] ~ 0 #5/4 * Θ[2],
            G[1] ~ 0 #-1/4 * k/dτ * Θ[2],
            G[2] ~ 0 #1/4 * Θ[2],
            [G[l] ~ 0 #=-l/(2*l+1) * k/dτ * ΘP[l-1]=# for l in 3:lmax]...    
        ] .|> O(ϵ^1))
    else
        append!(eqs1, [collect(G .~ 0)...] .|> O(ϵ^1)) # pin to zero
    end
    description = "Photon radiation"
    return extend(γ, ODESystem(eqs1, t, vars, []; initialization_eqs=ics1, defaults=defs, name, kwargs...); description)
end

"""
    massless_neutrinos(g; lmax = 6, name = :ν, kwargs...)

Create a particle species for massless neutrinos in the spacetime with metric `g`.
"""
function massless_neutrinos(g; lmax = 6, name = :ν, kwargs...)
    description = "Massless neutrinos"
    ν = radiation(g; name, description, kwargs...) |> thermodynamics |> complete

    vars = @variables F(t)[0:lmax+1] δ(t) θ(t) σ(t)
    pars = @parameters Neff
    eqs1 = [
        D(F[0]) ~ -k*F[1] + 4*D(g.Φ)
        D(F[1]) ~ k/3*(F[0]-2*F[2]+4*g.Ψ)
        [D(F[l]) ~ k/(2*l+1) * (l*F[l-1] - (l+1)*F[l+1]) for l in 2:lmax]...
        F[lmax+1] ~ (2*lmax+1) / (k*t) * F[lmax] - F[lmax-1]
        δ ~ F[0]
        θ ~ 3*k*F[1]/4
        σ ~ F[2]/2
        ν.cₛ² ~ 1//3
    ] .|> O(ϵ^1)
    ics1 = [
        δ ~ -2 * g.Ψ # adiabatic: δᵢ/(1+wᵢ) == δⱼ/(1+wⱼ) (https://cmb.wintherscoming.no/theory_initial.php#adiabatic)
        θ ~ 1//2 * (k^2*t) * g.Ψ
        σ ~ 1//15 * (k*t)^2 * g.Ψ
        [F[l] ~ 0 #=1/(2*l+1) * k*t * Θ[l-1]=# for l in 3:lmax]...
    ] .|> O(ϵ^1)
    description = "Massless neutrinos"
    return extend(ν, ODESystem(eqs1, t, vars, pars; initialization_eqs=ics1, name, kwargs...); description)
end

# TODO: use vector equations and simplify loops
"""
    massive_neutrinos(g; nx = 5, lmax = 4, name = :h, kwargs...)

Create a particle species for massive neutrinos in the spacetime with metric `g`.
"""
function massive_neutrinos(g; nx = 5, lmax = 4, name = :h, kwargs...)
    pars = @parameters m T₀ x[1:nx] W[1:nx] # Ω₀ ρ₀ y₀
    vars = @variables ρ(t) Ω(t) T(t) y(t) P(t) w(t) cₛ²(t) δ(t) σ(t) θ(t) u(t) E(t)[1:nx] ψ(t)[1:nx,0:lmax+1] In(t) Iρ(t) IP(t) Iδρ(t)

    f₀(x) = 1 / (exp(x) + 1) # not exp(E); distribution function is "frozen in"; see e.g. Dodelson exercise 3.9
    dlnf₀_dlnx(x) = -x / (1 + exp(-x))
    ∫dx_x²_f₀(f) = sum(collect(f .* W)) # a function that approximates the weighted integral ∫dx*x^2*f(x)*f₀(x)

    eqs0 = [
        T ~ T₀ / g.a
        y ~ m*c^2 / (kB*T)
        In ~ ∫dx_x²_f₀(1)
        Iρ ~ ∫dx_x²_f₀(E)
        IP ~ ∫dx_x²_f₀(x.^2 ./ E)
        ρ ~ 2/(2*π^2) * (kB*T)^4 / (ħ*c)^3 * Iρ / (g.H₀^2*c^2/GN) # compute g/(2π²ħ³) * ∫dp p² √((pc)² + (mc²)²) / (exp(pc/(kT)) + 1) with dimensionless x = pc/(kT) and degeneracy factor g = 2
        P ~ 2/(6*π^2) * (kB*T)^4 / (ħ*c)^3 * IP / (g.H₀^2*c^2/GN) # compute g/(6π²ħ³) * ∫dp p⁴ / √((pc)² + (mc²)²) / (exp(pc/(kT)) + 1) with dimensionless x = pc/(kT) and degeneracy factor g = 2
        w ~ P / ρ
        Ω ~ 8*Num(π)/3 * ρ
    ] .|> O(ϵ^0)
    eqs1 = [
        Iδρ ~ ∫dx_x²_f₀(E .* ψ[:,0])
        δ ~ Iδρ / Iρ
        u ~ ∫dx_x²_f₀(x .* ψ[:,1]) / (Iρ + IP/3)
        θ ~ u * k
        σ ~ (2//3) * ∫dx_x²_f₀(x.^2 ./ E .* ψ[:,2]) / (Iρ + IP/3)
        cₛ² ~ ∫dx_x²_f₀(x.^2 ./ E .* ψ[:,0]) / Iδρ # TODO: numerator ψ[:,0] or ψ[:,2]?
    ] .|> O(ϵ^1)
    # TODO: use defaults to set Ω₀ parameters etc. (e.g. with Λanalytical)?
    ics1 = []
    for i in 1:nx
        push!(eqs0, O(ϵ^0)(E[i] ~ √(x[i]^2 + y^2)))
        append!(eqs1, [
            D(ψ[i,0]) ~ -k * x[i]/E[i] * ψ[i,1] - D(g.Φ) * dlnf₀_dlnx(x[i])
            D(ψ[i,1]) ~ k/3 * x[i]/E[i] * (ψ[i,0] - 2*ψ[i,2]) - k/3 * E[i]/x[i] * g.Ψ * dlnf₀_dlnx(x[i])
            [D(ψ[i,l]) ~ k/(2*l+1) * x[i]/E[i] * (l*ψ[i,l-1] - (l+1)*ψ[i,l+1]) for l in 2:lmax]...
            ψ[i,lmax+1] ~ (2*lmax+1) * E[i]/x[i] * ψ[i,lmax] / (k*t) - ψ[i,lmax-1]
        ] .|> O(ϵ^1))
        append!(ics1, [
            ψ[i,0] ~ -1//4 * (-2*g.Ψ) * dlnf₀_dlnx(x[i])
            ψ[i,1] ~ -1/(3*k) * E[i]/x[i] * (1/2*(k^2*t)*g.Ψ) * dlnf₀_dlnx(x[i])
            ψ[i,2] ~ -1//2 * (1//15*(k*t)^2*g.Ψ) * dlnf₀_dlnx(x[i])
            [ψ[i,l] ~ 0 for l in 3:lmax] # TODO: full ICs
        ] .|> O(ϵ^1))
    end

    # compute numerical reduced momenta x = q*c / (kB*T) and Gaussian quadrature weights, and map the symbolic parameters to them
    _x, _W = gauss(x -> x^2 * f₀(x), nx, 0.0, 1e3)
    defs = [collect(x .=> _x); collect(W .=> _W)]

    description = "Massive neutrino"
    return ODESystem([eqs0; eqs1], t, vars, pars; initialization_eqs=ics1, defaults=defs, name, description, kwargs...)
end

"""
    cold_dark_matter(g; name = :c, kwargs...)

Create a particle species for cold dark matter in the spacetime with metric `g`.
"""
function cold_dark_matter(g; name = :c, kwargs...)
    c = matter(g; name, kwargs...) |> complete
    description = "Cold dark matter"
    c = extend(c, ODESystem([c.cₛ² ~ 0] .|> O(ϵ^1), t, [], []; name); description)
    return c
end

"""
    baryons(g; recombination=true, name = :b, kwargs...)

Create a particle species for baryons in the spacetime with metric `g`.
"""
function baryons(g; recombination = true, reionization = true, name = :b, kwargs...)
    description = "Baryons"
    b = matter(g; θinteract=true, name, description, kwargs...) |> complete
    if recombination # TODO: dont add recombination system when recombination = false
        @named rec = thermodynamics_recombination_recfast(g; reionization)
    else
        vars = @variables τ(t) τ̇(t) ρb(t) Tγ(t) cₛ²(t)
        @named rec = ODESystem([τ ~ 0, τ̇ ~ 0, cₛ² ~ 0], t, vars, [])
    end
    description = "Baryonic matter"
    b = extend(b, ODESystem([b.cₛ² ~ rec.cₛ²] .|> O(ϵ^1), t, [], []; name); description)
    b = compose(b, rec)
    return b
end

"""
    quintessence(g, v; name = :ϕ, kwargs...)

Create a species for a quintessence scalar field with potential `v` in the spacetime with metric `g`.
"""
function quintessence(g, v, v′, v′′; name = :Q, kwargs...)
    @variables ϕ(t) ρ(t) P(t) w(t) δ(t) σ(t) V(t) V′(t) V′′(t) K(t) m²(t) ϵs(t) ηs(t) cₛ²(t)
    eqs0 = [
        V ~ v
        V′ ~ v′
        V′′ ~ v′′
        K ~ (D(ϕ)/g.a)^2 / 2 # ϕ̇²/2 = (ϕ′/a)²/2
        D(D(ϕ)) ~ -2 * g.ℰ * D(ϕ) - g.a^2 * V′ # with cosmic time: ϕ̈ + 3*E*ϕ̇ + V′ = 0
        ρ ~ K + V
        P ~ K - V
        w ~ P / ρ
        m² ~ V′′
        ϵs ~ (V′/V)^2 / (16*Num(π)) # 1st slow roll parameter
        ηs ~ (V′′/V) / (8*Num(π)) # 2nd slow roll parameter
    ] .|> O(ϵ^0)
    eqs1 = [ # TODO: perturbations
        δ ~ 0
        σ ~ 0
        cₛ² ~ 0
    ] .|> O(ϵ^1)
    description = "Quintessence dark energy"
    return ODESystem([eqs0; eqs1], t; name, description, kwargs...)
end
function quintessence(g, v; name = :Q, kwargs...)
    @variables ϕ(t)
    ∂_∂ϕ = Differential(ϕ)
    v′ = ∂_∂ϕ(v(ϕ)) |> expand_derivatives |> simplify
    v′′ = ∂_∂ϕ(∂_∂ϕ(v(ϕ))) |> expand_derivatives |> simplify
    return quintessence(g, v(ϕ), v′, v′′; name, kwargs...)
end
