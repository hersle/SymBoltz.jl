using ModelingToolkit
using DifferentialEquations
using Plots

# TODO: try EnsembleProblems https://docs.sciml.ai/DiffEqDocs/dev/features/ensemble/

@variables a
@variables ρr(a) ρm(a) ρΛ(a) ρ(a)
@variables Pr(a) Pm(a) PΛ(a)
@variables Ωr(a) Ωm(a) ΩΛ(a)
@variables H(a)
eqs = [
    # Friedmann equation
    ρ ~ ρr + ρm + ρΛ
    H ~ √(8*π/3 * ρ)

    # equations of state
    Pr ~ ρr/3
    Pm ~ 0
    PΛ ~ -ρΛ

    # continuity equation
    Differential(a)(ρr) ~ -3 / a * (ρr + Pr)
    Differential(a)(ρm) ~ -3 / a * (ρm + Pm)
    Differential(a)(ρΛ) ~ -3 / a * (ρΛ + PΛ)

    # dimensionless density parameters
    Ωr ~ ρr / ρ
    Ωm ~ ρm / ρ
    ΩΛ ~ ρΛ / ρ
]
@mtkbuild universe = ODESystem(eqs, a)
prob = ODEProblem(
    universe, [Ωr => 1e-5, Ωm => 0.3, ΩΛ => 1 - Ωr - Ωm], (1.0, 1e-8), [];
    guesses = [ρr => 1, ρm => 1, ρΛ => 1, ρ => 1]
)
sol = solve(prob)

plot(log10.(sol[a]), [sol[Ωr], sol[Ωm], sol[ΩΛ]])
plot(log10.(sol[a]), log10.(sol[H] / sol[H][1]))