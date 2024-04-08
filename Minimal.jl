using ModelingToolkit
using DifferentialEquations
using ForwardDiff
using FiniteDiff
using Plots
Plots.default(label=nothing, markershape=:pixel)

# TODO: avoid this; make recombination work naturally with Unitful units
using PhysicalConstants, Unitful, UnitfulAstro
const c  = PhysicalConstants.CODATA2018.c_0 / u"m/s"
const h  = PhysicalConstants.CODATA2018.h / u"J*s"
const ħ  = PhysicalConstants.CODATA2018.ħ / u"J*s"
const kB = PhysicalConstants.CODATA2018.k_B / u"J/K"
const G  = PhysicalConstants.CODATA2018.G / u"m^3/kg/s^2"
const α  = PhysicalConstants.CODATA2018.α # fine structure constant, ≈ 1/137
const me = PhysicalConstants.CODATA2018.m_e / u"kg"
const mp = PhysicalConstants.CODATA2018.m_p / u"kg"
const mH = 1.6737236e-27 # kg
const σT = PhysicalConstants.CODATA2018.σ_e / u"m^2"
const km  = 1u"km/m"  |> NoUnits
const pc  = 1u"pc/m"  |> NoUnits
const Mpc = 1u"Mpc/m" |> NoUnits
const eV  = 1u"eV/J"  |> NoUnits
const EHion = 13.59844 * eV

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
    # TODO: replace with analytical solutions
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

# thermodynamics / recombination variables
# TODO: just merge with background?
@parameters fb H0 T0
@variables Xe(a) XH(a) Xp(a) α2(a) β(a) λe(a) ρb(a) nb(a) np(a) ne(a) nH(a) H(a) T(a) dτ(a) R(a)
thermo_eqs = [
    H ~ E * H0
    T ~ T0 / a # TODO: diff eq for temperature evolution?

    # TODO: add Peebles' corrections & He? see Dodelson exercise 4.7
    Da(Xe) ~ (XH*β - Xe*Xp*nb*α2) / (a*H) # same as ((1-Xe)*β - Xe^2*nb*α2) / (a*H); Xe ~ ne/nb; Dodelson (4.36) # TODO: nb or nH?
    α2 ~ 9.78 * (α*ħ/me)^2/c * √(EHion/(kB*T)) * log(EHion/(kB*T)) # Dodelson (4.38) (e⁻ + p → H + γ)
    β ~ α2 / λe^3 * exp(-EHion/(kB*T)) # Dodelson (4.37)-(4.38) (γ + H → e⁻ + p)
    λe ~ h / √(2π*me*kB*T) # electron de-Broglie wavelength

    ρb ~ fb * ρm # fb is baryon-to-matter fraction
    nb ~ ρb * 3*H0^2 / (8π*G) / (Xe*mp + (1-Xe)*mH) # ≈ ρb/mp * 3*H0^2 / (8π*G), Dodelson above (4.41) (ρb = ρp+ρH, nb=ρp+ρH)
    np ~ ne # charge neutrality
    nH ~ nb - ne # nb = nH + ne = nH + np
    ne ~ Xe * nb
    XH ~ nH / nb
    Xp ~ np / nb

    # optical depth
    dτ ~ -ne * σT * c / (a*H) # dτ = dτ/da
    R ~ 3/4 * ρb/ρr # Dodelson (5.74)

    # TODO: reionization?
]
@mtkbuild thermo = ODESystem([bg_eqs; thermo_eqs], a, [ρr, ρm, ρΛ, Xe], [H0, T0, fb])
thermo_prob = ODEProblem(thermo, [ρr, ρm, ρΛ, Xe] .=> NaN, (aini, atoday), [H0, T0, fb] .=> NaN)
function solve_thermodynamics(ρr0, ρm0, ρb0, H0)
    fb = ρb0 / ρm0; @assert fb <= 1
    T0 = (ρr0 * 15/π^2 * 3*H0^2/(8*π*G) * ħ^3*c^5)^(1/4) / kB
    bg_sol = solve_background(ρr0, ρm0)
    ρrini, ρmini, ρΛini = bg_sol(aini; idxs = [ρr, ρm, ρΛ]) # integrate background from atoday back to aini
    prob = remake(thermo_prob; u0=[ρrini, ρmini, ρΛini, 1], p=[H0, T0, fb]) # TODO: guarantee order!!
    return solve(prob, KenCarp4(), reltol=1e-9)
end

# perturbation variables
@parameters k # really k*c/H0
@variables Φ(a) Ψ(a) δρ(a) Θr0(a) Θr1(a) δc(a) δb(a) uc(a) ub(a) ρc(a) Δm(a)
pert_eqs = [
    # radiation perturbations (density & velocity)
    Da(Θr0) + k/(a^2*E)*Θr1 ~ -Da(Φ) # Dodelson (5.67) or (8.10)
    Da(Θr1) - k/(3*a^2*E)*Θr0 ~ k/(3*a^2*E)*Ψ + dτ * (Θr1 - ub/3) # Dodelson (5.67) or (8.11)

    # matter perturbations (density & velocity)
    ρc ~ ρm - ρb
    Da(δc) + k/(a^2*E)*uc ~ -3*Da(Φ) # Dodelson (5.69) or (8.12) with i*uc -> uc
    Da(uc) + uc/a ~ k/(a^2*E)*Ψ # Dodelson (5.70) or (8.13) with i*uc -> uc

    # baryon perturbations (density & velocity)
    Da(δb) + k/(a^2*E)*ub ~ -3*Da(Φ) # Dodelson (5.71) with i*ub -> ub
    Da(ub) + ub/a ~ k/(a^2*E)*Ψ + dτ/R * (ub - 3*Θr1)# Dodelson (5.72) with i*ub -> ub

    # gravity
    δρ ~ 4*ρr*Θr0 + δc*ρc + δb*ρb # total energy density perturbation
    Da(Φ) ~ (3/2*a^2*δρ - k^2*Φ - 3*a^2*E^2*Φ) / (3*a^3*E^2) # Dodelson (8.14) # TODO: write in more natural form?
    Δm ~ k^2*Φ / (3/2*a^2*ρm) # gauge-invariant overdensity (from Poisson equation)

    # anisotropic stress
    Ψ ~ -Φ # TODO: relax
]
@mtkbuild pert = ODESystem([bg_eqs; thermo_eqs; pert_eqs], a, [ρr, ρm, ρΛ, Xe, Φ, Θr0, δc, δb, Θr1, uc, ub], [H0, k, T0, fb])
pert_prob = ODEProblem(pert, [ρr, ρm, ρΛ, Xe, Φ, Θr0, δc, δb, Θr1, uc, ub] .=> NaN, (aini, atoday), [H0, k, T0, fb] .=> NaN; jac=true) # TODO: use remake

function solve_perturbations(kval, ρr0, ρm0, ρb0, H0)
    fb = ρb0 / ρm0; @assert fb <= 1 # TODO: avoid duplication thermo logic
    T0 = (ρr0 * 15/π^2 * 3*H0^2/(8*π*G) * ħ^3*c^5)^(1/4) / kB
    bg_sol = solve_background(ρr0, ρm0)
    ρrini, ρmini, ρΛini, Eini = bg_sol(aini; idxs = [ρr, ρm, ρΛ, E]) # integrate background from atoday back to aini
    Xeini = 1 # TODO: avoid duplication thermo logic
    Φini = 1 # arbitrary normalization (from primordial curvature power spectrum?)
    Θr0ini = Φini/2 # Dodelson (7.89)
    δcini = δbini = 3*Θr0ini # Dodelson (7.94)
    Θr1ini = -kval*Φini/(6*aini*Eini) # Dodelson (7.95) # TODO: replace aini -> a when this is fixed? https://github.com/SciML/ModelingToolkit.jl/issues/2543
    ucini = ubini = 3*Θr1ini # Dodelson (7.95)
    prob = remake(pert_prob; u0=[ρrini, ρmini, ρΛini, Xeini, Φini, Θr0ini, Θr1ini, δcini, ucini, δbini, ubini], p=[H0, kval, T0, fb]) # order checked with SymbolicIndexingInterface.parameter_index TODO: guarantee order!!
    return solve(prob, KenCarp4(), reltol=1e-10) # KenCarp4 and Kvaerno5 works well # TODO: use different EnsembleAlgorithm https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#Stiff-Problems
end

# power spectra
P0(k, As) = As / k^3
P(k, ρr0, ρm0, ρb0, H0, As) = P0(k, As) * solve_perturbations(k, ρr0, ρm0, ρb0, H0)(atoday; idxs=Δm)^2

if true
    ρr0 = 1e-5
    ρm0 = 0.3
    ρb0 = 0.02
    H0 = 70 * km/Mpc # s^-1
    As = 2e-9
    as = 10 .^ range(log10(aini), log10(atoday), length=400)
    k0 = 1 / 2997.92458 # h/Mpc
    ks = 10 .^ range(-4, +2, length=50) / k0 # in code units of k0 = H0/c

    bg_sol = solve_background(ρr0, ρm0)
    thermo_sol = solve_thermodynamics(ρr0, ρm0, ρb0, H0)
    pert_sols = [solve_perturbations(kval, ρr0, ρm0, ρb0, H0) for kval in ks] # TODO: use EnsembleProblem again

    # TODO: add plot recipe!
    p = plot(layout=(3,3), size=(1200, 900), margin=20*Plots.px); display(p)

    plot!(p[1], log10.(as), reduce(vcat, bg_sol.(as; idxs=[Ωr,Ωm,ΩΛ])'); xlabel="lg(a)", ylabel="Ω", label=["Ωr" "Ωm" "ΩΛ"], legend=:left); display(p)
    plot!(p[2], log10.(as), log10.(bg_sol.(as; idxs=E) / bg_sol(atoday; idxs=E)); xlabel="lg(a)", ylabel="lg(H/H0)"); display(p)
    plot!(p[3], log10.(as), log10.(thermo_sol.(as, idxs=Xe)); xlabel="lg(a)", ylabel="lg(Xe)"); display(p)
    plot!(p[4], log10.(as), [pert_sol.(as; idxs=Φ) for pert_sol in pert_sols]; xlabel="lg(a)", ylabel="Φ"); display(p)
    plot!(p[5], log10.(as), [[log10.(abs.(pert_sol.(as; idxs=δ))) for pert_sol in pert_sols] for δ in [δb,δc]]; color=[(1:length(ks))' (1:length(ks))'], xlabel="lg(a)", ylabel="lg(|δb|), lg(δc)"); display(p)

    # compute derivatives of output power spectrum with respect to input parameters
    derivatives = [FiniteDiff.finite_difference_derivative, ForwardDiff.derivative]
    plot!(p[6], log10.(ks*k0), [[derivative.(lgk -> log10(P(10^lgk, ρr0, ρm0, ρb0, H0, As)), log10.(ks)) for derivative in derivatives]..., log10.(P.(ks, ρr0, ρm0, ρb0, H0, As)/k0^3)]; xlabel="lg(k/(h/Mpc))", label=["d lg(P) / d lg(k) (auto. diff.)" "d lg(P) / d lg(k) (fin. diff.)" "lg(P/(Mpc/h)³)"], legend=:bottomleft); display(p)
    plot!(p[7], log10.(ks*k0), [derivative(lgρr0 -> log10.(P.(ks, 10^lgρr0, ρm0, ρb0, H0, As)), log10(ρr0)) for derivative in derivatives]; xlabel="lg(k/(h/Mpc))", ylabel="d lg(P) / d lg(Ωr0)", labels=["fin. diff." "auto. diff."]); display(p)
    plot!(p[8], log10.(ks*k0), [derivative(lgρm0 -> log10.(P.(ks, ρr0, 10^lgρm0, ρb0, H0, As)), log10(ρm0)) for derivative in derivatives]; xlabel="lg(k/(h/Mpc))", ylabel="d lg(P) / d lg(Ωm0)", labels=["fin. diff." "auto. diff."]); display(p)
    plot!(p[9], log10.(ks*k0), [derivative(lgAs  -> log10.(P.(ks, ρr0, ρm0, ρb0, H0, 10^lgAs)), log10(As)) for derivative in derivatives]; xlabel="lg(k/(h/Mpc))", ylabel="d lg(P) / d lg(As)", labels=["fin. diff." "auto. diff."], ylims=(0, 2)); display(p)
end