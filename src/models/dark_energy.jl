"""
    cosmological_constant(g; name = :Λ, kwargs...)

Create a particle species for the cosmological constant (with equation of state `w ~ -1`) in the spacetime with metric `g`.
"""
function cosmological_constant(g; name = :Λ, analytical = true, kwargs...)
    description = "Cosmological constant"
    Λ = species_constant_eos(g, -1; name, analytical, adiabatic = true, description, kwargs...) |> background |> complete # discard ill-defined perturbations
    vars = @variables begin
        δ(τ, k), [description = "Overdensity"]
        θ(τ, k), [description = "Velocity divergence"]
        σ(τ, k), [description = "Shear stress"]
    end
    eqs = [δ ~ 0, θ ~ 0, σ ~ 0]
    description = "Cosmological constant"
    return extend(Λ, System(eqs, τ, vars, []; name); description) # manually set perturbations to zero
end

"""
    w0wa(g; kwargs...)

Create a particle species for the w₀-wₐ dark energy (CPL) parametrization in the spacetime with metric `g`.
"""
function w0wa(g; name = :X, analytical = false, kwargs...)
    pars = @parameters begin
        w0, [description = "Equation of state today"]
        wa, [description = "Equation of state evolution"]
        cₛ², [description = "Rest-frame speed of sound squared"]
    end
    vars = @variables begin
        ρ(τ), [description = "Background density"]
        P(τ), [description = "Background pressure"]
        w(τ), [description = "Equation of state"]
        ẇ(τ), [description = "Equation of state derivative"]
        cₐ²(τ), [description = "Adiabatic speed of sound squared"]
        δ(τ, k), [description = "Overdensity"]
        θ(τ, k), [description = "Velocity divergence"]
        σ(τ, k), [description = "Shear stress"]
    end
    # TODO: generate equations with a generic species_eos function
    eqs = [
        w ~ w0 + wa * (1 - g.a)
        ẇ ~ D(w)
        P ~ w * ρ
    ]
    if analytical
        append!(pars, @parameters Ω₀)
        push!(eqs, ρ ~ 3/(8*Num(π))*Ω₀ * abs(g.a)^(-3 * (1 + w0 + wa)) * exp(-3 * wa * (1 - g.a))) # energy density # TODO: get rid of abs
    else
        push!(eqs, D(ρ) ~ -3 * g.ℰ * ρ * (1 + w))
    end
    append!(eqs, [
        # Following https://arxiv.org/pdf/1002.1311 section II
        cₐ² ~ w - ẇ/(3*g.ℰ*(1+w))
        D(δ) ~ -(1+w)*(θ-3*D(g.Φ)) - 3*g.ℰ*(cₛ²-w)*δ - 9*(g.ℰ/k)^2*(1+w)*(cₛ²-cₐ²)*θ
        D(θ) ~ -g.ℰ*(1-3*cₛ²)*θ + cₛ²/(1+w)*k^2*δ - k^2*σ + k^2*g.Ψ
        σ ~ 0
    ])
    ics = [
        δ ~ -3//2 * (1+w) * g.Ψ # adiabatic ICs, see e.g. https://arxiv.org/abs/1004.5509 eq. (3.17)
        θ ~ 1//2 * (k^2*τ) * g.Ψ # τ ≈ 1/ℰ; adiabatic ICs, see e.g. https://arxiv.org/abs/1004.5509 eq. (3.18)
    ]
    description = "w₀wₐ (CPL) dark energy"
    return System(eqs, τ, vars, pars; initialization_eqs=ics, name, description, kwargs...)
end

"""
    quintessence(g, v; name = :ϕ, kwargs...)

Create a species for a quintessence scalar field with potential `v` in the spacetime with metric `g`.
"""
function quintessence(g, v, v′, v′′; name = :Q, kwargs...)
    @variables begin
        ϕ(τ), [description = "Background scalar field"]
        ρ(τ), [description = "Effective background density"]
        P(τ), [description = "Effective background pressure"]
        w(τ), [description = "Equation of state"]
        δ(τ, k), [description = "Overdensity"]
        σ(τ, k), [description = "Shear stress"]
        V(τ), [description = "Potential of scalar field"]
        V′(τ), [description = "Potential derivative wrt. scalar field"]
        V′′(τ), [description = "Potential 2nd derivative wrt. scalar field"]
        K(τ), [description = "Effective kinetic energy"]
        m²(τ), [description = "Effective mass"]
        ϵs(τ), [description = "1st slow roll parameter"]
        ηs(τ), [description = "2nd slow roll parameter"]
        cₛ²(τ), [description = "Speed of sound squared"]
    end
    eqs = [
        V ~ v
        V′ ~ v′
        V′′ ~ v′′
        K ~ (D(ϕ)/g.a)^2 / 2 # ϕ̇²/2 = (ϕ′/a)²/2
        D(D(ϕ)) ~ -2 * g.ℰ * D(ϕ) - g.a^2 * V′ # with cosmic time: ϕ̈ + 3*E*ϕ̇ + V′ = 0
        ρ ~ K + V
        P ~ K - V
        w ~ P / ρ
        m² ~ V′′
        ϵs ~ (V′/V)^2 / (16*Num(π))
        ηs ~ (V′′/V) / (8*Num(π))

        δ ~ 0
        σ ~ 0
        cₛ² ~ 0
    ]
    description = "Quintessence dark energy"
    return System(eqs, τ; name, description, kwargs...)
end
function quintessence(g, v; name = :Q, kwargs...)
    @variables begin
        ϕ(τ), [description = "Background scalar field"]
    end
    ∂_∂ϕ = Differential(ϕ)
    v′ = ∂_∂ϕ(v(ϕ)) |> expand_derivatives |> simplify
    v′′ = ∂_∂ϕ(∂_∂ϕ(v(ϕ))) |> expand_derivatives |> simplify
    return quintessence(g, v(ϕ), v′, v′′; name, kwargs...)
end
