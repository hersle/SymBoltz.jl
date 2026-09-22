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
        δ(τ, k), [description = "Overdensity (gauge-dependent)"]
        Δ(τ, k), [description = "Overdensity (gauge-independent)"]
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
        push!(eqs, D(ρ) ~ -3 * g.ℋ * ρ * (1 + w))
    end
    append!(eqs, [
        # Following https://arxiv.org/pdf/1002.1311 section II
        cₐ² ~ w - ẇ/(3*g.ℋ*(1+w))
        D(δ) ~ -(1+w)*(θ-3*D(g.Φ)) - 3*g.ℋ*(cₛ²-w)*δ - 9*(g.ℋ/k)^2*(1+w)*(cₛ²-cₐ²)*θ
        D(θ) ~ -g.ℋ*(1-3*cₛ²)*θ + cₛ²/(1+w)*k^2*δ - k^2*σ + k^2*g.Ψ
        Δ ~ δ + 3*g.ℋ*(1+w)*θ/k^2
        σ ~ 0
    ])
    ieqs = [
        δ ~ -3//2 * (1+w) * g.Ψ # adiabatic ICs, see e.g. https://arxiv.org/abs/1004.5509 eq. (3.17)
        θ ~ 1//2 * (k^2*τ) * g.Ψ # τ ≈ 1/ℋ; adiabatic ICs, see e.g. https://arxiv.org/abs/1004.5509 eq. (3.18)
    ]
    description = "w₀wₐ (CPL) dark energy"
    return System(eqs, τ, vars, pars; initialization_eqs = ieqs, name, description, kwargs...)
end

"""
    quintessence(g, V; name = :Q, kwargs...)

Create a species with a quintessence scalar field in the potential `V(ϕ)` in the spacetime with metric `g`.
The reduced density today `Ω₀` is imposed as a shooting condition, so one parameter (e.g. `ϕini` or one in `V`) must be shot for in `CosmologyProblem`.
"""
function quintessence(g, V; name = :Q, kwargs...)
    v = V # the potential function; V is redefined as a variable below
    pars = @parameters begin
        Ω₀, [description = "Reduced background density today"]
        ϕini, [description = "Initial scalar field"]
        ϕ′ini = 0.0, [description = "Initial scalar field conformal time derivative"]
    end
    vars = @variables begin
        ϕ(τ), [description = "Background scalar field"]
        ρ(τ), [description = "Effective background density"]
        P(τ), [description = "Effective background pressure"]
        w(τ), [description = "Equation of state"]
        δϕ(τ, k), [description = "Scalar field perturbation"]
        δρ(τ, k), [description = "Density perturbation"]
        δP(τ, k), [description = "Pressure perturbation"]
        δ(τ, k), [description = "Overdensity"]
        σ(τ, k), [description = "Shear stress"]
        V(τ), [description = "Potential of scalar field"]
        V′(τ), [description = "Potential derivative wrt. scalar field"]
        V′′(τ), [description = "Potential 2nd derivative wrt. scalar field"] # ″ does not show properly in Mathjax
        K(τ), [description = "Effective kinetic energy"]
        m²(τ), [description = "Effective mass"]
        ϵs(τ), [description = "1st slow roll parameter"]
        ηs(τ), [description = "2nd slow roll parameter"]
        cₛ²(τ, k), [description = "Effective speed of sound squared"]
    end
    ∂_∂ϕ = Differential(ϕ)
    eqs = [
        V ~ v(ϕ)
        V′ ~ simplify(expand_derivatives(∂_∂ϕ(v(ϕ))))
        V′′ ~ simplify(expand_derivatives(∂_∂ϕ(∂_∂ϕ(v(ϕ)))))
        K ~ (D(ϕ)/g.a)^2 / 2 # ϕ̇²/2 = (ϕ′/a)²/2
        D(D(ϕ)) ~ -2 * g.ℋ * D(ϕ) - g.a^2 * V′ # with cosmic time: ϕ̈ + 3*H*ϕ̇ + V′ = 0
        ρ ~ K + V
        P ~ K - V
        w ~ P / ρ
        m² ~ V′′
        ϵs ~ (V′/V)^2 / (16*Num(π))
        ηs ~ (V′′/V) / (8*Num(π))

        # perturbed Klein-Gordon equation and effective fluid perturbations in the conformal Newtonian gauge
        D(D(δϕ)) ~ -2*g.ℋ*D(δϕ) - (k^2 + g.a^2*V′′)*δϕ + D(ϕ)*(g.Ψ̇ + 3*g.Φ̇) - 2*g.a^2*V′*g.Ψ
        δρ ~ (D(ϕ)*D(δϕ) - D(ϕ)^2*g.Ψ) / g.a^2 + V′*δϕ
        δP ~ (D(ϕ)*D(δϕ) - D(ϕ)^2*g.Ψ) / g.a^2 - V′*δϕ
        δ ~ δρ / ρ
        cₛ² ~ δP / δρ # so δP = cₛ² * δρ
        σ ~ 0 # no anisotropic stress
    ]
    append!(pars, filter(ModelingToolkit.isparameter, Symbolics.get_variables(v(ϕ)))) # collect any additional parameters in the potential function
    initial_conditions = [ϕ => ϕini]
    initialization_eqs = [
        D(ϕ) ~ ϕ′ini
        δϕ ~ D(ϕ)*g.Ψ*τ/2 # adiabatic ICs with same time shift as photons: δϕ/ϕ′ = δργ/ργ′ = Ψ/(2ℋ) with ℋ ≈ 1/τ
        D(δϕ) ~ (D(D(ϕ))*τ + D(ϕ))*g.Ψ/2 # time derivative of the above with Ψ constant on superhorizon scales
    ]
    constraints = [ρ ~ 3/(8*Num(π)) * Ω₀] # density today
    description = "Quintessence dark energy"
    return System(eqs, τ, vars, pars; initial_conditions, initialization_eqs, constraints, name, description, kwargs...)
end
