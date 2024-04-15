using ModelingToolkit
using DifferentialEquations
using SymbolicIndexingInterface
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

const solver = KenCarp4()

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
# TODO: relate parameters through parameter expressions: https://docs.sciml.ai/ModelingToolkit/stable/basics/Composition/#Variable-scope-and-parameter-expressions
# TODO: define components with @mtkmodel?

# independent variable: scale factor
# TODO: spacetime/geometry structure?
@variables a E(a)
E = GlobalScope(E)
Da = Differential(a)
aini, atoday = 1e-8, 1e0

function background_gravity_GR(; name)
    @variables ρ(a)
    return ODESystem([
        E ~ √(ρ) # H/H₀
    ], a; name)
end

function background_species_constant_eos(w; name)
    @variables ρ(a) P(a) Ω(a) ρcrit(a)
    return ODESystem([
        P ~ w*ρ
        Da(ρ) ~ -3/a * (ρ + P) # TODO: replace with analytical solutions and ρ0 parameter
        #ρ ~ ρ0 / a^(3*(1+w)) # TODO: recombination doesn't trigger with this analytical solution, at least when parametrized by ln(a). maybe it works when parametrized by a?
        Ω ~ ρ / ρcrit
    ], a; name)
end

@named rad = background_species_constant_eos(1//3)
@named mat = background_species_constant_eos(0)
@named de = background_species_constant_eos(-1)
@named grav = background_gravity_GR()

@named bg = ODESystem([
    grav.ρ ~ rad.ρ + mat.ρ + de.ρ
    rad.ρcrit ~ grav.ρ
    mat.ρcrit ~ grav.ρ
    de.ρcrit ~ grav.ρ
], a)

@named bg = compose(bg, rad, mat, de, grav) # TODO: extend?
bg_sim = structural_simplify(bg)
bg_prob = ODEProblem(bg_sim, unknowns(bg_sim) .=> NaN, (atoday, aini))
function solve_background(ρr0, ρm0)
    ρΛ0 = 1 - ρr0 - ρm0 # TODO: handle with equation between parameters once ρr0 etc. are parameters?
    prob = remake(bg_prob; u0 = [bg_sim.rad.ρ => ρr0, bg_sim.mat.ρ => ρm0, bg_sim.de.ρ => ρΛ0]) # TODO: bg.rad.ρ => ρr0 etc. doesn't work. bug?
    return solve(prob, solver, reltol=1e-8)
end

# thermodynamics / recombination variables
# TODO: just merge with background?
function thermodynamics(; name)
    @parameters fb H0 T0
    @variables ρr(a) ρm(a) Xe(a) ρb(a) nb(a) np(a) ne(a) nH(a) H(a) T(a) dτ(a) R(a) α2(a) β(a) λe(a) C(a) Λα(a) Λ2γ(a) β2(a) SahaXe(a) SahaB(a)
    return ODESystem([
        T ~ T0 / a # TODO: diff eq for temperature evolution?

        # Saha approximation # TODO: move to separate component
        SahaB ~ exp(-EHion/(kB*T)) / (nb*λe^3)
        SahaXe ~ ifelse(SahaB < 1e10, SahaB/2 * (-1 + √(1+4/SahaB)), 1 - 1/SahaB) # Taylor expansion valid when SB ≫ 1

        α2 ~ 9.78 * (α*ħ/me)^2/c * √(EHion/(kB*T)) * log(EHion/(kB*T)) # Dodelson (4.38) (e⁻ + p → H + γ)
        β ~ α2 / λe^3 * exp(-EHion/(kB*T)) # Dodelson (4.37)-(4.38) (γ + H → e⁻ + p)
        λe ~ h / √(2π*me*kB*T) # electron de-Broglie wavelength

        # Peebles' correction factor (Dodelson exercise 4.7)
        Λα ~ H * (3*EHion/(ħ*c))^3 / ((8*π)^2 * nH) # 1/s
        Λ2γ ~ 8.227 # 1/s
        β2 ~ α2 / λe^3 * exp(-EHion/(4*kB*T)) # 1/s (compute this instead of β2 = β * exp(3*EHion/(4*kB*T)) to avoid exp overflow)
        C ~ (Λ2γ + Λα) / (Λ2γ + Λα + β2)

        Da(Xe) ~ ifelse(Xe > 0.99, Da(SahaXe), C * ((1-Xe)*β - Xe^2*nb*α2) / (a*H)) # Xe ~ ne/nb; Dodelson (4.36) # TODO: nb or nH?

        ρb ~ fb * ρm # fb is baryon-to-matter fraction
        nb ~ ρb / mp
        np ~ ne # charge neutrality
        nH ~ nb - ne # nb = nH + ne = nH + np
        ne ~ Xe * nb

        # optical depth
        dτ ~ -ne * σT * c / (a*H) # dτ = dτ/da
        R ~ 3/4 * ρb/ρr # Dodelson (5.74)

        # TODO: reionization?
    ], a, [Xe, H, ρr, ρm, ρb, dτ, R], [fb, H0, T0]; name)
end
@named th = thermodynamics()
@named th_bg_conn = ODESystem([
    th.H ~ E * th.H0 # 1/s
    th.ρr ~ bg.rad.ρ * 3*th.H0^2 / (8*π*G) # kg/m³
    th.ρm ~ bg.mat.ρ * 3*th.H0^2 / (8*π*G) # kg/m³
], a)
@named th_bg = compose(th_bg_conn, th, bg)
th_sim = structural_simplify(th_bg)
th_prob = ODEProblem(th_sim, unknowns(th_sim) .=> NaN, (aini, atoday); jac=true)
function solve_thermodynamics(ρr0, ρm0, ρb0, H0)
    fb = ρb0 / ρm0; @assert fb <= 1
    T0 = (ρr0 * 15/π^2 * 3*H0^2/(8*π*G) * ħ^3*c^5)^(1/4) / kB # TODO: relate to ρr0 once that is a parameter
    ρrini, ρmini, ρΛini = solve_background(ρr0, ρm0)(aini; idxs = [bg.rad.ρ, bg.mat.ρ, bg.de.ρ]) # integrate background from atoday back to aini # TODO: avoid when ρr0 etc. are parameters
    prob = remake(th_prob; u0 = [th.Xe => 1, bg.rad.ρ => ρrini, bg.mat.ρ => ρmini, bg.de.ρ => ρΛini], p = [th.fb => fb, th.H0 => H0, th.T0 => T0])
    return solve(prob, solver, reltol=1e-8) # TODO: after switching ivar from a to b=ln(a), the integrator needs more steps. fix this?
end

@variables Φ(a) Ψ(a)
@parameters k # perturbation wavenumber # TODO: associate like pt.k
k, Φ, Ψ = GlobalScope.([k, Φ, Ψ])

function perturbations_radiation(interact=false; name)
    @variables Θ0(a) Θ1(a) δ(a)
    interaction = interact ? only(@variables interaction(a)) : 0
    return ODESystem([
        Da(Θ0) + k/(a^2*E)*Θ1 ~ -Da(Φ) # Dodelson (5.67) or (8.10)
        Da(Θ1) - k/(3*a^2*E)*Θ0 ~ k/(3*a^2*E)*Ψ + interaction # Dodelson (5.67) or (8.11)
        δ ~ 4*Θ0
    ], a; name)    
end

function perturbations_matter(interact=false; name)
    @variables δ(a) u(a)
    interaction = interact ? only(@variables interaction(a)) : 0
    return ODESystem([
        Da(δ) + k/(a^2*E)*u ~ -3*Da(Φ) # Dodelson (5.69) or (8.12) with i*uc -> uc
        Da(u) + u/a ~ k/(a^2*E)*Ψ + interaction # Dodelson (5.70) or (8.13) with i*uc -> uc
    ], a; name)
end

function perturbations_gravity(; name)
    @variables δρ(a) Δm(a) ρm(a)
    return ODESystem([
        Da(Φ) ~ (3/2*a^2*δρ - k^2*Φ - 3*(a*E)^2*Φ) / (3*a^3*E^2) # Dodelson (8.14) # TODO: write in more natural form?
        Ψ ~ -Φ # anisotropic stress # TODO: relax
        Δm ~ k^2*Φ / (3/2*a^2*ρm) # gauge-invariant overdensity (from Poisson equation) # TODO: move outside or change?
    ], a; name)
end

@named rad = perturbations_radiation(true)
@named cdm = perturbations_matter(false)
@named bar = perturbations_matter(true)
@named grav = perturbations_gravity()
@variables ρc(a) ρb(a) δρr(a) δρc(a) δρb(a) # TODO: get rid of
@named pt_th_bg_conn = ODESystem([
    ρb ~ th.fb * bg.mat.ρ
    ρc ~ bg.mat.ρ - ρb
    δρr ~ rad.δ * bg.rad.ρ
    δρc ~ cdm.δ * ρc
    δρb ~ bar.δ * ρb
    grav.δρ ~ δρr + δρc + δρb # total energy density perturbation
    grav.ρm ~ bg.mat.ρ

    # baryon-photon interactions: Compton (Thomson) scattering # TODO: define connector type?
    rad.interaction ~ -th.dτ/3    * (bar.u - 3*rad.Θ1)
    bar.interaction ~ +th.dτ/th.R * (bar.u - 3*rad.Θ1)
], a)
pt_th_bg_conn = extend(pt_th_bg_conn, th_bg_conn)

@named pt = compose(pt_th_bg_conn, rad, bar, cdm, grav, th, bg)
pt_sim = structural_simplify(pt)
pt_prob = ODEProblem(pt_sim, unknowns(pt_sim) .=> NaN, (aini, atoday); jac=true) 

function solve_perturbations(kval, ρr0, ρm0, ρb0, H0)
    println("k = $(kval*k0) Mpc/h")
    fb = ρb0 / ρm0; @assert fb <= 1 # TODO: avoid duplication thermo logic
    T0 = (ρr0 * 15/π^2 * 3*H0^2/(8*π*G) * ħ^3*c^5)^(1/4) / kB
    bg_sol = solve_background(ρr0, ρm0)
    ρrini, ρmini, ρΛini, Eini = bg_sol(aini; idxs = [bg.rad.ρ, bg.mat.ρ, bg.de.ρ, E]) # integrate background from atoday back to aini
    Φini = 1.0 # arbitrary normalization (from primordial curvature power spectrum?)
    Θr0ini = Φini/2 # Dodelson (7.89)
    δcini = δbini = 3*Θr0ini # Dodelson (7.94)
    Θr1ini = -kval*Φini/(6*aini*Eini) # Dodelson (7.95) # TODO: replace aini -> a when this is fixed? https://github.com/SciML/ModelingToolkit.jl/issues/2543
    ucini = ubini = 3*Θr1ini # Dodelson (7.95)
    prob = remake(pt_prob; u0 = [Φ => Φini, pt_sim.rad.Θ0 => Θr0ini, pt_sim.rad.Θ1 => Θr1ini, pt_sim.bar.δ => δbini, pt_sim.bar.u => ubini, pt_sim.cdm.δ => δcini, pt_sim.cdm.u => ucini, th.Xe => 1, bg.rad.ρ => ρrini, bg.mat.ρ => ρmini, bg.de.ρ => ρΛini], p = [th.fb => fb, k => kval, th.H0 => H0, th.T0 => T0])
    return solve(prob, solver, reltol=1e-8) # KenCarp4 and Kvaerno5 works well # TODO: use different EnsembleAlgorithm https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#Stiff-Problems
end

# power spectra
P0(k, As) = As / k^3
P(k, ρr0, ρm0, ρb0, H0, As) = P0(k, As) * solve_perturbations(k, ρr0, ρm0, ρb0, H0)(atoday; idxs=pt.grav.Δm)^2
P(x) = P(10 ^ x[1], 10^x[2], 10^x[3], 10^x[4], 10^x[5], 10^x[6]) # x = log10.([k, ρr0, ρm0, ρb0, H0, As])

if true
    ρr0 = 5e-5
    ρm0 = 0.3
    ρb0 = 0.02
    H0 = 70 * km/Mpc # s^-1
    As = 2e-9
    as = 10 .^ range(log10(aini), log10(atoday), length=400)
    k0 = 1 / 2997.92458 # h/Mpc
    ks = 10 .^ range(-4, +2, length=200) / k0 # in code units of k0 = H0/c

    # TODO: add plot recipe!
    p = plot(layout=(3,3), size=(1200, 900), margin=20*Plots.px); display(p)

    bg_sol = solve_background(ρr0, ρm0)
    plot!(p[1], log10.(as), reduce(vcat, bg_sol.(as; idxs=[bg_sim.rad.Ω,bg_sim.mat.Ω,bg_sim.de.Ω])'); xlabel="lg(a)", ylabel="Ω", label=["Ωr" "Ωm" "ΩΛ"], legend=:left); display(p)
    plot!(p[2], log10.(as), log10.(bg_sol.(as; idxs=E) / bg_sol(atoday; idxs=E)); xlabel="lg(a)", ylabel="lg(H/H0)"); display(p)

    th_sol = solve_thermodynamics(ρr0, ρm0, ρb0, H0)
    plot!(p[3], log10.(as), log10.(th_sol.(as, idxs=th.Xe)); xlabel="lg(a)", ylabel="lg(Xe)"); display(p)

    pt_sols = [solve_perturbations(kval, ρr0, ρm0, ρb0, H0) for kval in ks] # TODO: use EnsembleProblem again
    plot!(p[4], log10.(as), [pt_sol.(as; idxs=Φ) for pt_sol in pt_sols]; xlabel="lg(a)", ylabel="Φ"); display(p)
    plot!(p[5], log10.(as), [[log10.(abs.(pt_sol.(as; idxs=δ))) for pt_sol in pt_sols] for δ in [pt.bar.δ,pt.cdm.δ]]; color=[(1:length(ks))' (1:length(ks))'], xlabel="lg(a)", ylabel="lg(|δb|), lg(δc)"); display(p)

    x0 = [log10(ρr0), log10(ρm0), log10(ρb0), log10(H0), log10(As)]
    Ps = stack([P([log10(k), x0...]) for k in ks]) # TODO: use pert_sols
    plot!(p[6], log10.(ks*k0), log10.(Ps/k0^3); xlabel="lg(k/(h/Mpc))", label="lg(P/(Mpc/h)³)", color=3, legend=:bottomleft); display(p)

    # compute derivatives of output power spectrum with respect to input parameters
    for (derivative, name, color) in [((f, x) -> FiniteDiff.finite_difference_gradient(f, x; relstep=1e-2), "fin. diff.", 1), ((f, x) -> ForwardDiff.gradient(f, x), "auto. diff.", 2)]
        dlgP_dxs = stack([derivative(x -> log10(P(x)), [log10(k), x0...]) for k in ks])
        plot!(p[6], log10.(ks*k0), dlgP_dxs[1,:]; xlabel="lg(k/(h/Mpc))", label="d lg(P) / d lg(k) ($name)", color=color, legend=:bottomleft); display(p)
        plot!(p[7], log10.(ks*k0), dlgP_dxs[2,:]; xlabel="lg(k/(h/Mpc))", ylabel="d lg(P) / d lg(Ωr0)", color=color, label=name); display(p)
        plot!(p[8], log10.(ks*k0), dlgP_dxs[3,:]; xlabel="lg(k/(h/Mpc))", ylabel="d lg(P) / d lg(Ωm0)", color=color, label=name); display(p)
        plot!(p[9], log10.(ks*k0), dlgP_dxs[6,:]; xlabel="lg(k/(h/Mpc))", ylabel="d lg(P) / d lg(As)", color=color, label=name, ylims=(0, 2)); display(p)
    end
end