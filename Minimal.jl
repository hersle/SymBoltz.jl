using ModelingToolkit
using DifferentialEquations
using ForwardDiff
using FiniteDiff
using Plots
Plots.default(label=nothing, markershape=:pixel)

# TODO: shooting method https://docs.sciml.ai/DiffEqDocs/stable/tutorials/bvp_example/ (not supported by ModelingToolkit: https://github.com/SciML/ModelingToolkit.jl/issues/924, https://discourse.julialang.org/t/boundary-value-problem-with-modellingtoolkit-or-diffeqoperators/57656)
# TODO: register thermodynamics functions: https://docs.sciml.ai/ModelingToolkit/stable/tutorials/ode_modeling/#Specifying-a-time-variable-forcing-function
# TODO: make simpler Cosmology interface
# TODO: compare runtime for finite vs. dlgP_dlgks_autodiff
# TODO: compare accuracy with class
# TODO: non-linear: higher-order perturbations vs halofit vs N-body?
# TODO: CMB power spectrum
# TODO: baryons: Recfast -> Recfast++ -> CosmoRec -> HyRec -> HyRec-2: call out, or integrate equations into my code to make use of my background calculation?
# TODO: composable models, generate equations
# TODO: fix remake/replace etc. with new ModelingToolkit version
# TODO: modified gravity: (coupled) quintessence, Brans-Dicke, DGP, parametrized framework, EFT of LSS, ...
# TODO: analytical solutions of e.g. background ρ evolution
# TODO: GPU-parallellized EnsembleProblem
# TODO: use b = ln(a) as independent variable?

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
@variables Φ(a) Ψ(a) δρ(a) Θr0(a) Θr1(a) δm(a) um(a) Δm(a)
pert_eqs = [
    # radiation perturbations (density & velocity)
    Da(Θr0) + k/(a^2*E)*Θr1 ~ -Da(Φ) # Dodelson (5.67) or (8.10)
    Da(Θr1) - k/(3*a^2*E)*Θr0 ~ k/(3*a^2*E)*Ψ # Dodelson (5.67) or (8.11)

    # matter perturbations (density & velocity)
    Da(δm) + k/(a^2*E)*um ~ -3*Da(Φ) # Dodelson (5.69) or (8.12) with i*uc -> um
    Da(um) + um/a ~ k/(a^2*E)*Ψ # Dodelson (5.70) or (8.13) with i*uc -> um

    # gravity
    δρ ~ ρm*δm + 4*ρr*Θr0 # total energy density perturbation
    Da(Φ) ~ (3/2*a^2*δρ - k^2*Φ - 3*a^2*E^2*Φ) / (3*a^3*E^2) # Dodelson (8.14) # TODO: write in more natural form?
    Δm ~ k^2*Φ / (3/2*a^2*ρm) # gauge-invariant overdensity (from Poisson equation)

    # anisotropic stress
    Ψ ~ -Φ # TODO: relax
]
@mtkbuild pert = ODESystem([bg_eqs; pert_eqs], a, [ρr, ρm, ρΛ, Φ, Θr0, δm, Θr1, um], [k])
pert_prob = ODEProblem(pert, [ρr, ρm, ρΛ, Φ, Θr0, δm, Θr1, um] .=> NaN, (aini, atoday), [k => NaN]; jac=true) # TODO: use remake

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

# power spectra
P0(k, As) = As / k^3
P(k, ρr0, ρm0, As) = P0(k, As) * solve_perturbations(k, ρr0, ρm0)(atoday; idxs=Δm)^2

if true
    ρr0 = 1e-5
    ρm0 = 0.3
    As = 2e-9
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
    
    # compute derivatives of output power spectrum with respect to input parameters
    derivatives = [FiniteDiff.finite_difference_derivative, ForwardDiff.derivative]
    p6 = plot(log10.(ks*k0), [[derivative.(lgk -> log10(P(10^lgk, ρr0, ρm0, As)), log10.(ks)) for derivative in derivatives]..., log10.(P.(ks, ρr0, ρm0, As)/k0^3)]; xlabel="lg(k/(h/Mpc))", label=["d lg(P) / d lg(k) (auto. diff.)" "d lg(P) / d lg(k) (fin. diff.)" "lg(P/(Mpc/h)³)"], legend=:bottomleft)
    p7 = plot(log10.(ks*k0), [derivative(lgρr0 -> log10.(P.(ks, 10^lgρr0, ρm0, As)), log10(ρr0)) for derivative in derivatives]; xlabel="lg(k/(h/Mpc))", ylabel="d lg(P) / d lg(Ωr0)", labels=["fin. diff." "auto. diff."])
    p8 = plot(log10.(ks*k0), [derivative(lgρm0 -> log10.(P.(ks, ρr0, 10^lgρm0, As)), log10(ρm0)) for derivative in derivatives]; xlabel="lg(k/(h/Mpc))", ylabel="d lg(P) / d lg(Ωm0)", labels=["fin. diff." "auto. diff."])
    p9 = plot(log10.(ks*k0), [derivative(lgAs  -> log10.(P.(ks, ρr0, ρm0, 10^lgAs)), log10(As)) for derivative in derivatives]; xlabel="lg(k/(h/Mpc))", ylabel="d lg(P) / d lg(As)", labels=["fin. diff." "auto. diff."], ylims=(0, 2))

    p = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, layout=(3,3), size=(1200, 900), margin=20*Plots.px) 
    display(p)
end