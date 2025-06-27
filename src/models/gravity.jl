"""
    general_relativity(g; name = :G, kwargs...)

Create a symbolic component for the general relativistic (GR) theory of gravity in the spacetime with the metric `g`.
"""
function general_relativity(g; acceleration = false, name = :G, kwargs...)
    vars = @variables begin
        ρ(τ), [description = "Total background density"]
        P(τ), [description = "Total background pressure"]
        ρcrit(τ), [description = "Critical background density"]
        δρ(τ, k), [description = "Total density perturbation"]
        δP(τ, k), [description = "Total pressure perturbation"]
        Π(τ, k), [description = "Anisotropic stress perturbation"]
        F₁(τ), [description = "1st Friedmann equation residual"]
        F₂(τ), [description = "2nd Friedmann equation residual"]
    end
    a = g.a
    if acceleration
        pars = @parameters begin
            Δfac = 0, [description = "Constraint dampening factor of 1st Friedmann equation"]
        end
        eqs = [
            D(D(a)) ~ D(a)^2/(2*a)*(1-3*P/ρ) + Δfac*F₁ # use constraint damping (alternatively D(D(a)) ~ D(a)^2/a - 4*Num(π)/3*(ρ+3*P)*a^3)) # TODO: incorporate sign(∂Δ/∂a′), but positive in GR?
            F₁ ~ D(a)^2 - 8*Num(π)/3*ρ*a^4 # violation of (squared) Friedmann constraint equation
            F₂ ~ 0 # we are enforcing the 2nd Friedmann equation
        ]
        ics = [F₁ ~ 0]
    else
        pars = []
        eqs = [
            D(a) ~ √(8*Num(π)/3 * ρ) * a^2 # "normal" Friedmann equation (+√ of (squared) constraint equation)
            F₁ ~ 0 # we are enforcing the 1st Friedmann equation
            F₂ ~ D(D(a)) - D(a)^2/(2*a)*(1-3*P/ρ) # violation of acceleration equation
        ]
        ics = []
    end
    append!(eqs, [
        D(g.Φ) ~ -4*Num(π)/3*a^2/g.ℰ*δρ - k^2/(3*g.ℰ)*g.Φ - g.ℰ*g.Ψ
        k^2 * (g.Φ - g.Ψ) ~ 12*Num(π) * a^2 * Π
        g.Ψ̇ ~ D(g.Ψ)
        g.Φ̇ ~ D(g.Φ)
    ])
    guesses = [ρ => 0.1, D(a) => +1]
    description = "General relativity gravity"
    return System(eqs, τ, vars, pars; initialization_eqs = ics, guesses, name, description, kwargs...)
end

# TODO: potential
# TODO: mass
"""
    brans_dicke(g; name = :G, kwargs...)

Create a symbolic component for the Brans-Dicke (BD) theory of gravity in the spacetime with the metric `g`.
"""
function brans_dicke(g; name = :G, acceleration = false, kwargs...)
    pars = @parameters begin
        ω, [description = "Brans-Dicke coupling constant"] 
    end
    vars = @variables begin
        ρ(τ), [description = "Total background density"]
        P(τ), [description = "Total background pressure"]
        ϕ(τ), [description = "Brans-Dicke background scalar field"]
        δϕ(τ, k), [description = "Brans-Dicke scalar field perturbation"]
        δρ(τ, k), [description = "Total density perturbation"]
        δP(τ, k), [description = "Total pressure perturbation"]
        Π(τ, k), [description = "Anisotropic stress perturbation"]
        F₁(τ), [description = "1st Friedmann equation residual"]
        F₂(τ), [description = "2nd Friedmann equation residual"]
        G(τ), [description = "Effective gravitational strength"]
    end
    @unpack a, ℰ, Φ, Ψ = g # shorthand
    F1 = D(a)^2 ~ 8*Num(π)/3*ρ*a^4/ϕ - D(a)*a*D(ϕ)/ϕ + ω/6*a^2*(D(ϕ)/ϕ)^2
    F2 = D(D(a)) ~ D(a)^2/(2*a) - 4*Num(π)*a^3*P/ϕ - ω/4*a*(D(ϕ)/ϕ)^2 - D(a)/2*D(ϕ)/ϕ - a/2*D(D(ϕ))/ϕ
    KG = D(D(ϕ)) ~ 8*Num(π)/(2*ω+3) * a^2 * (ρ-3*P) - 2*D(a)/a*D(ϕ)
    eqs = [
        KG
        G ~ (2*ω+4) / (2*ω+3) / ϕ # effective gravitational constant
    ]
    ics = []
    if acceleration
        append!(eqs, [
            F2
            F₂ ~ 0
            F₁ ~ F1.lhs - F1.rhs # violation of Friedmann contraint
        ])
        push!(ics, F1)
    else
        append!(eqs, [
            D(a) ~ -a/2*D(ϕ)/ϕ + √((a/2*D(ϕ)/ϕ)^2 + 8*Num(π)/3*ρ*a^4/ϕ + ω/6*(a*D(ϕ)/ϕ)^2) # solve quadratic equation for ȧ # TODO: Symbolics.symbolic_solve
            F₁ ~ 0
            F₂ ~ F2.lhs - F2.rhs # violation of acceleration equation
        ])
    end
    append!(eqs, [
        # https://arxiv.org/pdf/2006.04273 equations (D.20, D.22, D.25) with Φ → Ψ and Ψ → -Φ
        D(Φ) ~ -(8*Num(π)*a^2*(2*Ψ*ρ+δρ) + 2*k^2*Φ*ϕ - (3*ℰ^2+k^2)*δϕ - 3*ℰ*D(δϕ) - ω/2*(δϕ*(D(ϕ)/ϕ)^2-2*D(ϕ)*D(δϕ))) / (6*ℰ*ϕ + 3*D(ϕ)) # (μ,ν) = (0,0)
        Ψ ~ Φ - δϕ/ϕ - 12*Num(π)/ϕ * a^2 * Π / k^2 # (μ,ν) = (i,j), i ≠ j, add matter stress divided by ϕ from field equations # TODO: correct?
        D(D(δϕ)) ~ -8*Num(π)/(3+2*ω)*a^2*(3*δP-δρ) + 2*D(D(ϕ))*Ψ - k^2*δϕ - 2*ℰ*D(δϕ) + D(ϕ)*(D(Ψ)+3*D(ϕ)) + 4*ℰ*Ψ*D(ϕ) # perturbed Klein-Gordon equation
        g.Ψ̇ ~ D(g.Ψ)
        g.Φ̇ ~ D(g.Φ)
    ])
    ics = append!(ics, [
        δϕ ~ 0.0 # TODO: set properly
        D(δϕ) ~ 0.0 # TODO: set properly
    ])
    guesses = [ρ => 1.0, D(g.a) => +1.0]
    description = "Brans-Dicke gravity"
    return System(eqs, τ, vars, pars; name, description, initialization_eqs = ics, guesses, kwargs...) # TODO: don't pass vars and pars
end
