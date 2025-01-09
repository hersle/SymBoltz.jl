"""
    general_relativity(g; name = :G, kwargs...)

Create a symbolic component for the general relativistic (GR) theory of gravity in the spacetime with the metric `g`.
"""
function general_relativity(g; acceleration = false, name = :G, kwargs...)
    vars = @variables ρ(t) P(t) ρcrit(t) δρ(t) δP(t) Π(t) Δ(t)
    pars = @parameters Δfac = 0
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
    guesses = [ρ => 0.1, D(a) => +1]
    description = "General relativity gravity"
    return ODESystem([eqs0; eqs1], t, vars, pars; initialization_eqs = ics0, guesses, name, description, kwargs...)
end

# TODO: potential
# TODO: mass
"""
    brans_dicke(g; name = :G, kwargs...)

Create a symbolic component for the Brans-Dicke (BD) theory of gravity in the spacetime with the metric `g`.
"""
function brans_dicke(g; name = :G, acceleration = false, kwargs...)
    pars = @parameters ω
    vars = @variables ρ(t) P(t) ϕ(t) δϕ(t) δρ(t) δP(t) Π(t) Δ(t) G(t)
    @unpack a, ℰ, Φ, Ψ = g # shorthand
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
            D(a) ~ -a/2*D(ϕ)/ϕ + √((a/2*D(ϕ)/ϕ)^2 + 8*Num(π)/3*ρ*a^4/ϕ + ω/6*(a*D(ϕ)/ϕ)^2) # solve quadratic equation for ȧ # TODO: Symbolics.symbolic_solve
            Δ ~ F2.lhs - F2.rhs # violation of acceleration equation
        ] .|> O(ϵ^0))
    end
    eqs1 = [
        # https://arxiv.org/pdf/2006.04273 equations (D.20, D.22, D.25) with Φ → Ψ and Ψ → -Φ
        D(Φ) ~ -(8*Num(π)*a^2*(2*Ψ*ρ+δρ) + 2*k^2*Φ*ϕ - (3*ℰ^2+k^2)*δϕ - 3*ℰ*D(δϕ) - ω/2*(δϕ*(D(ϕ)/ϕ)^2-2*D(ϕ)*D(δϕ))) / (6*ℰ*ϕ + 3*D(ϕ)) # (μ,ν) = (0,0)
        Ψ ~ Φ - δϕ/ϕ - 12*Num(π)/ϕ * a^2 * Π / k^2 # (μ,ν) = (i,j), i ≠ j, add matter stress divided by ϕ from field equations # TODO: correct?
        D(D(δϕ)) ~ -8*Num(π)/(3+2*ω)*a^2*(3*δP-δρ) + 2*D(D(ϕ))*Ψ - k^2*δϕ - 2*ℰ*D(δϕ) + D(ϕ)*(D(Ψ)+3*D(ϕ)) + 4*ℰ*Ψ*D(ϕ) # perturbed Klein-Gordon equation
    ] .|> O(ϵ^1)
    ics1 = [
        δϕ ~ 0.0 # TODO: set properly
        D(δϕ) ~ 0.0 # TODO: set properly
    ] .|> O(ϵ^1)
    guesses = [ρ => 1.0, D(g.a) => +1.0]
    description = "Brans-Dicke gravity"
    return ODESystem([eqs0; eqs1], t, vars, pars; name, description, initialization_eqs = [ics0; ics1], guesses, kwargs...) # TODO: don't pass vars and pars
end
