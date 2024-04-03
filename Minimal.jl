using ModelingToolkit
using DifferentialEquations
using ForwardDiff
using FiniteDiff
using Plots
Plots.default(label=nothing, markershape=:pixel)

# TODO: shooting method https://docs.sciml.ai/DiffEqDocs/stable/tutorials/bvp_example/ (not supported by ModelingToolkit: https://github.com/SciML/ModelingToolkit.jl/issues/924, https://discourse.julialang.org/t/boundary-value-problem-with-modellingtoolkit-or-diffeqoperators/57656)

# independent variable: scale factor
@variables a
Da = Differential(a)
aini, atoday = 1e-8, 1e0

# background variables
@variables E(a) ρ(a) ρr(a) ρm(a) ρΛ(a) Pr(a) Pm(a) PΛ(a) Ωr(a) Ωm(a) ΩΛ(a)
bg_eqs = [
    # Friedmann equation
    ρ ~ ρr + ρm + ρΛ # energy densities and pressures are in units of 3H₀²c² / 8πG
    E ~ √(ρ) # H/H₀

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
@mtkbuild bg = ODESystem(bg_eqs, a, [ρr, ρm, ρΛ], [])
bg_prob = ODEProblem(bg, [ρr, ρm, ρΛ] .=> NaN, (atoday, aini), []) # TODO: integrate from aini; make a "shooting interface" for different parameters
function solve_background(ρr0, ρm0) # TODO: define parameter struct?
    ρΛ0 = 1 - ρr0 - ρm0
    ics = ModelingToolkit.varmap_to_vars([ρr => ρr0, ρm => ρm0, ρΛ => ρΛ0], unknowns(bg))  # https://docs.sciml.ai/ModelingToolkit/dev/basics/FAQ/#Transforming-value-maps-to-arrays
    prob = remake(bg_prob; u0=ics) # https://github.com/SciML/SymbolicIndexingInterface.jl/issues/59
    return solve(prob)
end

# perturbation variables
@parameters k # really k*c/H0
@variables Φ(a) Θr0(a) Θr1(a) δm(a) um(a) Δm(a) # TODO: Ψ ≠ Φ
pert_eqs = [
    Da(Θr0) + k/(a^2*E)*Θr1 ~ -Da(Φ) # Dodelson (8.10)
    Da(Θr1) - k/(3*a^2*E)*Θr0 ~ -k/(3*a^2*E)*Φ # Dodelson (8.11)
    Da(δm) + k/(a^2*E)*um ~ -3*Da(Φ) # Dodelson (8.12) with i*uc -> um
    Da(um) + um/a ~ -k/(a^2*E)*Φ # Dodelson (8.13) with i*uc -> um
    Da(Φ) ~ (3/2*a^2*(ρm*δm + 4*ρr*Θr0) - k^2*Φ - 3*a^2*E^2*Φ) / (3*a^3*E^2) # Dodelson (8.14) # TODO: write in more natural form?
    Δm ~ k^2*Φ / (3/2*a^2*ρm) # gauge-invariant overdensity (from Poisson equation)
]
@mtkbuild pert = ODESystem([bg_eqs; pert_eqs], a, [ρr, ρm, ρΛ, Φ, Θr0, δm, Θr1, um], [k])
pert_prob = ODEProblem(pert, [ρr, ρm, ρΛ, Φ, Θr0, δm, Θr1, um] .=> NaN, (aini, atoday), [k => NaN]) # TODO: use remake

function solve_perturbations(kval, ρr0, ρm0)
    bg_sol = solve_background(ρr0, ρm0)
    ρrini, ρmini, ρΛini, Eini = bg_sol(aini; idxs = [ρr, ρm, ρΛ, E]) # integrate background from atoday back to aini
    Φini = 1 # arbitrary normalization (from primordial curvature power spectrum?)
    Θr0ini = Φini/2 # Dodelson (7.89)
    δmini = 3*Θr0ini # Dodelson (7.94)
    Θr1ini = -kval*Φini/(6*aini*Eini) # Dodelson (7.95) # TODO: replace aini -> a when this is fixed? https://github.com/SciML/ModelingToolkit.jl/issues/2543
    umini = 3*Θr1ini # Dodelson (7.95)
    prob = remake(pert_prob; u0=[ρrini, ρmini, ρΛini, Φini, Θr0ini, Θr1ini, δmini, umini], p=[kval]) # TODO: guarantee order!!
    return solve(prob, KenCarp4(), reltol=1e-7) # KenCarp4 and Kvaerno5 works well # TODO: use different EnsembleAlgorithm https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#Stiff-Problems
end

function P(kval, ρr0, ρm0)
    P0 = 1 / kval^3 # primordial power spectrum
    pert_sol = solve_perturbations(kval, ρr0, ρm0)
    return P0 * pert_sol(atoday; idxs=Δm)^2
end

if true
    ρr0 = 1e-5
    ρm0 = 0.3
    as = 10 .^ range(log10(aini), log10(atoday), length=400)
    k0 = 1 / 2997.92458 # h/Mpc
    ks = 10 .^ range(-4, +2, length=50) / k0 # in code units of k0 = H0/c

    bg_sol = solve_background(ρr0, ρm0)
    pert_sols = [solve_perturbations(kval, ρr0, ρm0) for kval in ks] # TODO: use EnsembleProblem again

    # TODO: add plot recipe!
    p1 = plot(log10.(as), reduce(vcat, bg_sol.(as; idxs=[Ωr,Ωm,ΩΛ])'); xlabel="lg(a)", ylabel="Ω", label=["Ωr" "Ωm" "ΩΛ"], legend=:left)
    p2 = plot(log10.(as), log10.(bg_sol.(as; idxs=E) / bg_sol(atoday; idxs=E)); xlabel="lg(a)", ylabel="lg(H/H0)")
    p3 = plot(; xlabel="???", ylabel="???") # TODO
    p4 = plot(log10.(as), [pert_sol.(as; idxs=Φ) for pert_sol in pert_sols]; xlabel="lg(a)", ylabel="Φ")
    p5 = plot(log10.(as), [log10.(pert_sol.(as; idxs=δm)) for pert_sol in pert_sols]; xlabel="lg(a)", ylabel="lg(δm)")
    
    lgP(lgk) = log10(P(10^lgk, ρr0, ρm0))
    lgPs = lgP.(log10.(ks))
    dlgP_dlgks_autodiff = ForwardDiff.derivative.(lgP, log10.(ks))
    dlgP_dlgks_findiff = FiniteDiff.finite_difference_derivative.(lgP, log10.(ks))
    p6 = plot(; xlabel="lg(k/(h/Mpc))", legend=:bottom)
    plot!(p6, log10.(ks*k0), [lgPs, dlgP_dlgks_findiff, dlgP_dlgks_autodiff]; label=["lg(P)" "d lg(P) / d lg(k) (fin. diff.)" "d lg(P) / d lg(k) (auto. diff.)"])

    p = plot(p1, p2, p3, p4, p5, p6, layout=(2,3), size=(1200, 600), margin=20*Plots.px) 
    display(p)
end