using ModelingToolkit
using DifferentialEquations
using Plots
Plots.default(label=nothing, markershape=:pixel)

# independent variable: scale factor
@variables a
Da = Differential(a)
aini, atoday = 1e-8, 1e0

# background variables
@variables ρr(a) ρm(a) ρΛ(a) ρ(a)
@variables Pr(a) Pm(a) PΛ(a)
@variables Ωr(a) Ωm(a) ΩΛ(a)
@variables H(a)

bg_eqs = [
    # Friedmann equation
    ρ ~ ρr + ρm + ρΛ
    H ~ √(8*π/3 * ρ)

    # equations of state
    Pr ~ ρr/3
    Pm ~ 0
    PΛ ~ -ρΛ

    # continuity equation
    Da(ρr) ~ -3/a * (ρr + Pr)
    Da(ρm) ~ -3/a * (ρm + Pm)
    Da(ρΛ) ~ -3/a * (ρΛ + PΛ)

    # dimensionless density parameters
    Ωr ~ ρr / ρ
    Ωm ~ ρm / ρ
    ΩΛ ~ ρΛ / ρ
]
@mtkbuild bg = ODESystem(bg_eqs, a)
bg_ics = [Ωr => 1e-5, Ωm => 0.3, ΩΛ => 1 - Ωr - Ωm] # TODO: enforce sum(Ω) == 1 in equations?
bg_prob = ODEProblem(bg, bg_ics, (atoday, aini), []; guesses = [ρr => 1, ρm => 1, ρΛ => 1, ρ => 1]) # TODO: integrate from aini; make a "shooting interface" for different parameters
bg_sol = solve(bg_prob)

@parameters k
@variables Φ(a) # TODO: Ψ ≠ Φ
@variables Θr0(a) Θr1(a)
@variables δm(a) um(a)
@variables δρr(a) δρm(a) δρ(a) Δm(a) # TODO: δr and ur?

pert_eqs = [
    Da(Θr0) + k/(a^2*H)*Θr1 ~ -Da(Φ) # Dodelson (8.10)
    Da(Θr1) - k/(3*a^2*H)*Θr0 ~ -k/(3*a^2*H)*Φ # Dodelson (8.11)
    Da(δm) + k/(a^2*H)*um ~ -3*Da(Φ) # Dodelson (8.12) with i*uc -> um
    Da(um) + um/a ~ -k/(a^2*H)*Φ # Dodelson (8.13) with i*uc -> um
    Da(Φ) ~ (4π*a^2*(ρm*δm + 4*ρr*Θr0) - k^2*Φ - 3*a^2*H^2*Φ) / (3*a^3*H^2) # Dodelson (8.14) # TODO: write in more natural form?
    Δm ~ k^2*Φ / (4π*a^2*ρm) # gauge-invariant overdensity (from Poisson equation)
]
ρrini, ρmini, ρΛini, Hini = bg_sol(aini; idxs = [ρr, ρm, ρΛ, H])
@mtkbuild pert = ODESystem([bg_eqs; pert_eqs], a)
pert_ics = [
    ρr => ρrini, ρm => ρmini, ρΛ => ρΛini, # from background solution
    Φ => 1, # arbitrary normalization (from primordial curvature power spectrum?)
    Θr0 => Φ/2, # Dodelson (7.89)
    δm => 3*Θr0, # Dodelson (7.94)
    Θr1 => -k*Φ/(6*aini*Hini), # Dodelson (7.95) # TODO: can replace aini -> a when this is fixed: https://github.com/SciML/ModelingToolkit.jl/issues/2543
    um => 3*Θr1 # Dodelson (7.95)
]
pert_prob = ODEProblem(pert, pert_ics, (aini, atoday), [k => 0])

P0(k) = 1 / k^3 # primordial power spectrum
ks = 10 .^ range(-3, +6, length=30)
pert_probs = EnsembleProblem(pert_prob; prob_func = (prob, i, _) -> remake(prob; p = [k => ks[i]]))
pert_sols = solve(pert_probs, Rodas5P(), trajectories=length(ks), reltol=1e-3) # TODO: use different EnsembleAlgorithm https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#Stiff-Problems
Δms = map(sol -> sol(atoday; idxs=Δm), pert_sols)
Ps = @. P0(ks) * Δms^2 # present power spectrum

p1 = plot(log10.(bg_sol[a]), [bg_sol[Ωr], bg_sol[Ωm], bg_sol[ΩΛ]]; xlabel="lg(a)", label=["Ωr" "Ωm" "ΩΛ"])
p2 = plot(log10.(bg_sol[a]), log10.(bg_sol[H] / bg_sol(atoday; idxs=H)); xlabel="lg(a)", ylabel="H/H₀")
p3 = plot() # TODO
p4 = plot(map(sol -> (log10.(sol[a]), sol[Φ]), pert_sols); xlabel="lg(a)", ylabel="Φ")
p5 = plot(map(sol -> (log10.(sol[a]), log10.(sol[δm])), pert_sols); xlabel="lg(a)", ylabel="δm")
p6 = plot(log10.(ks), log10.(Ps); xlabel="k/H0", ylabel="lg(P)")
p = plot(p1, p2, p3, p4, p5, p6, layout=(2,3), size=(1200, 600), margin=20*Plots.px)
display(p)