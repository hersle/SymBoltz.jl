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
# TODO: relate parameters through parameter expressions: https://docs.sciml.ai/ModelingToolkit/stable/basics/Composition/#Variable-scope-and-parameter-expressions
# TODO: define components with @mtkmodel?

# independent variable: scale factor
@variables b
Db = Differential(b)
aini, atoday = 1e-8, 1e0
bini, btoday = log(aini), log(atoday)

function background_spacetime(; name)
    @variables a(b)
    return ODESystem([
        a ~ exp(b)
    ], b; name)
end

function background_gravity_GR(; name)
    @variables E(b) ρ(b)
    return ODESystem([
        E ~ √(ρ) # H/H₀
    ], b; name)
end

function background_species_constant_eos(w; name)
    @variables ρ(b) P(b) Ω(b) ρcrit(b)
    return ODESystem([
        P ~ w*ρ
        Db(ρ) ~ -3 * (ρ + P) # TODO: replace with analytical solutions
        Ω ~ ρ / ρcrit
    ], b; name)
end

@named st = background_spacetime()
@named rad = background_species_constant_eos(1//3)
@named mat = background_species_constant_eos(0)
@named de = background_species_constant_eos(-1)
@named grav = background_gravity_GR()

@named bg = ODESystem([
    grav.ρ ~ rad.ρ + mat.ρ + de.ρ
    rad.ρcrit ~ grav.ρ
    mat.ρcrit ~ grav.ρ
    de.ρcrit ~ grav.ρ
], b)

@named bg = compose(bg, st, rad, mat, de, grav) # TODO: extend?
bg_simple = structural_simplify(bg)
bg_prob = ODEProblem(bg_simple, [rad.ρ, mat.ρ, de.ρ] .=> NaN, (btoday, bini), [])
function solve_background(ρr0, ρm0)
    ρΛ0 = 1 - ρr0 - ρm0
    prob = remake(bg_prob; u0=[ρr0, ρm0, ρΛ0]) # https://github.com/SciML/SymbolicIndexingInterface.jl/issues/59
    return solve(prob)
end

# thermodynamics / recombination variables
# TODO: just merge with background?
function thermodynamics(; name)
    @parameters fb H0 T0
    @variables a(b) ρr(b) ρm(b) Xe(b) XH(b) Xp(b) α2(b) β(b) λe(b) ρb(b) nb(b) np(b) ne(b) nH(b) H(b) T(b) dτ(b) R(b)
    return ODESystem([
        T ~ T0 / a # TODO: diff eq for temperature evolution?

        # TODO: add Peebles' corrections & He? see Dodelson exercise 4.7
        Db(Xe) ~ ((1-Xe)*β - Xe^2*nb*α2) / H # Xe ~ ne/nb; Dodelson (4.36) # TODO: nb or nH?
        α2 ~ 9.78 * (α*ħ/me)^2/c * √(EHion/(kB*T)) * log(EHion/(kB*T)) # Dodelson (4.38) (e⁻ + p → H + γ)
        β ~ α2 / λe^3 * exp(-EHion/(kB*T)) # Dodelson (4.37)-(4.38) (γ + H → e⁻ + p)
        λe ~ h / √(2π*me*kB*T) # electron de-Broglie wavelength

        ρb ~ fb * ρm # fb is baryon-to-matter fraction
        nb ~ ρb * 3*H0^2 / (8π*G) / (Xe*mp + (1-Xe)*mH) # ≈ ρb/mp * 3*H0^2 / (8π*G), Dodelson above (4.41) (ρb = ρp+ρH, nb=ρp+ρH)
        np ~ ne # charge neutrality
        nH ~ nb - ne # nb = nH + ne = nH + np
        ne ~ Xe * nb

        # optical depth
        dτ ~ -ne * σT * c / H # dτ = dτ/db
        R ~ 3/4 * ρb/ρr # Dodelson (5.74)

        # TODO: reionization?
    ], b; name)
end
@named th = thermodynamics()
@named th_bg_conn = ODESystem([
    th.a ~ bg.st.a
    th.H ~ bg.grav.E * th.H0
    th.ρr ~ bg.rad.ρ
    th.ρm ~ bg.mat.ρ
], b)
@named th_bg = compose(th_bg_conn, th, bg)
th_simple = structural_simplify(th_bg)
th_prob = ODEProblem(th_simple, [th.Xe, bg.rad.ρ, bg.mat.ρ, bg.de.ρ] .=> NaN, (bini, btoday), [th.fb, th.H0, th.T0] .=> NaN; jac=true)
function solve_thermodynamics(ρr0, ρm0, ρb0, H0)
    fb = ρb0 / ρm0; @assert fb <= 1
    T0 = (ρr0 * 15/π^2 * 3*H0^2/(8*π*G) * ħ^3*c^5)^(1/4) / kB
    bg_sol = solve_background(ρr0, ρm0)
    ρrini, ρmini, ρΛini = bg_sol(bini; idxs = [bg.rad.ρ, bg.mat.ρ, bg.de.ρ]) # integrate background from atoday back to aini
    prob = remake(th_prob; u0=[1, ρrini, ρmini, ρΛini], p=[fb, H0, T0]) # TODO: guarantee order!!
    return solve(prob, KenCarp4(), reltol=1e-15) # TODO: after switching ivar from a to b=ln(a), the integrator needs more steps. fix this?
end

@parameters k # perturbation wavenumber # TODO: associate like pt.k
k = GlobalScope(k)

function perturbations_spacetime(; name)
    @variables Φ(b) Ψ(b)
    return ODESystem(Equation[], b, [Φ, Ψ], []; name)
end

function perturbations_photons(; name)
    # TODO: represent interactions better (using fewer variables)?
    @variables Θ0(b) Θ1(b) E(b) a(b) dτ(b) ub(b) Φ(b) Ψ(b)
    return ODESystem([
        Db(Θ0) + k/(a*E)*Θ1 ~ -Db(Φ) # Dodelson (5.67) or (8.10)
        Db(Θ1) - k/(3*a*E)*Θ0 ~ k/(3*a*E)*Ψ + dτ * (Θ1 - ub/3) # Dodelson (5.67) or (8.11)    
    ], b; name)
end

function perturbations_baryons(; name)
    @variables δ(b) u(b) E(b) a(b) dτ(b) R(b) Θr1(b) Φ(b) Ψ(b)
    return ODESystem([
        Db(δ) + k/(a*E)*u ~ -3*Db(Φ) # Dodelson (5.71) with i*ub -> ub
        Db(u) + u ~ k/(a*E)*Ψ + dτ/R * (u - 3*Θr1)# Dodelson (5.72) with i*ub -> ub
    ], b; name)
end

function perturbations_cold_dark_matter(; name)
    @variables δ(b) u(b) E(b) a(b) Φ(b) Ψ(b)
    return ODESystem([
        Db(δ) + k/(a*E)*u ~ -3*Db(Φ) # Dodelson (5.69) or (8.12) with i*uc -> uc
        Db(u) + u ~ k/(a*E)*Ψ # Dodelson (5.70) or (8.13) with i*uc -> uc
    ], b; name)
end

function perturbations_gravity(; name)
    @variables δρ(b) Δm(b) E(b) a(b) ρm(b) Φ(b) Ψ(b)
    return ODESystem([
        Db(Φ) ~ (3/2*a^2*δρ - k^2*Φ - 3*(a*E)^2*Φ) / (3*(a*E)^2) # Dodelson (8.14) # TODO: write in more natural form?
        Ψ ~ -Φ # anisotropic stress # TODO: relax
        Δm ~ k^2*Φ / (3/2*a^2*ρm) # gauge-invariant overdensity (from Poisson equation) # TODO: move outside or change?
    ], b; name)
end

@named st = perturbations_spacetime()
@named rad = perturbations_photons()
@named cdm = perturbations_cold_dark_matter()
@named bar = perturbations_baryons()
@named grav = perturbations_gravity()
@variables ρc(b) # TODO: get rid of
@named pt_th_bg_conn = ODESystem([
    ρc ~ bg.mat.ρ - th.ρb
    grav.δρ ~ 4*bg.rad.ρ*rad.Θ0 + cdm.δ*ρc + bar.δ*th.ρb # total energy density perturbation

    # TODO: simplify this mess!
    rad.a ~ bg.st.a
    rad.E ~ bg.grav.E
    rad.Φ ~ st.Φ
    rad.Ψ ~ st.Ψ
    rad.ub ~ bar.u
    rad.dτ ~ th.dτ

    bar.a ~ bg.st.a
    bar.E ~ bg.grav.E
    bar.Φ ~ st.Φ
    bar.Ψ ~ st.Ψ
    bar.dτ ~ th.dτ
    bar.R ~ th.R
    bar.Θr1 ~ rad.Θ1
    
    cdm.a ~ bg.st.a
    cdm.E ~ bg.grav.E
    cdm.Φ ~ st.Φ
    cdm.Ψ ~ st.Ψ

    grav.E ~ bg.grav.E
    grav.a ~ bg.st.a
    grav.Φ ~ st.Φ
    grav.Ψ ~ st.Ψ
    grav.ρm ~ bg.mat.ρ
], b)
pt_th_bg_conn = extend(pt_th_bg_conn, th_bg_conn)

@named pt = compose(pt_th_bg_conn, st, rad, bar, cdm, grav, th, bg) # TODO
pt_simple = structural_simplify(pt)
pt_prob = ODEProblem(pt_simple, [grav.Φ, rad.Θ0, rad.Θ1, bar.δ, bar.u, cdm.δ, cdm.u, th.Xe, bg.rad.ρ, bg.mat.ρ, bg.de.ρ] .=> NaN, (bini, btoday), [th.fb, k, th.H0, th.T0] .=> NaN; jac=true) 

function solve_perturbations(kval, ρr0, ρm0, ρb0, H0)
    println("k = $(kval*k0) Mpc/h")
    fb = ρb0 / ρm0; @assert fb <= 1 # TODO: avoid duplication thermo logic
    T0 = (ρr0 * 15/π^2 * 3*H0^2/(8*π*G) * ħ^3*c^5)^(1/4) / kB
    bg_sol = solve_background(ρr0, ρm0)
    ρrini, ρmini, ρΛini, Eini = bg_sol(bini; idxs = [bg.rad.ρ, bg.mat.ρ, bg.de.ρ, bg.grav.E]) # integrate background from atoday back to aini
    Xeini = 1.0 # TODO: avoid duplication thermo logic
    Φini = 1.0 # arbitrary normalization (from primordial curvature power spectrum?)
    Θr0ini = Φini/2 # Dodelson (7.89)
    δcini = δbini = 3*Θr0ini # Dodelson (7.94)
    Θr1ini = -kval*Φini/(6*aini*Eini) # Dodelson (7.95) # TODO: replace aini -> a when this is fixed? https://github.com/SciML/ModelingToolkit.jl/issues/2543
    ucini = ubini = 3*Θr1ini # Dodelson (7.95)
    prob = remake(pt_prob; u0=[Φini, Θr0ini, Θr1ini, δbini, ubini, δcini, ucini, Xeini, ρrini, ρmini, ρΛini], p=[fb, kval, H0, T0]) # order checked with SymbolicIndexingInterface.parameter_index (on the ODEProblem!!!) TODO: guarantee order!!
    return solve(prob, KenCarp4(), reltol=1e-10) # KenCarp4 and Kvaerno5 works well # TODO: use different EnsembleAlgorithm https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#Stiff-Problems
end

# power spectra
P0(k, As) = As / k^3
P(k, ρr0, ρm0, ρb0, H0, As) = P0(k, As) * solve_perturbations(k, ρr0, ρm0, ρb0, H0)(btoday; idxs=pt.grav.Δm)^2
P(x) = P(10 ^ x[1], 10^x[2], 10^x[3], 10^x[4], 10^x[5], 10^x[6]) # x = log10.([k, ρr0, ρm0, ρb0, H0, As])

if true
    ρr0 = 1e-5
    ρm0 = 0.3
    ρb0 = 0.1
    H0 = 70 * km/Mpc # s^-1
    As = 2e-9
    bs = range(bini, btoday, length=400)
    as = exp.(bs)
    k0 = 1 / 2997.92458 # h/Mpc
    ks = 10 .^ range(-4, +2, length=200) / k0 # in code units of k0 = H0/c

    # TODO: add plot recipe!
    p = plot(layout=(3,3), size=(1200, 900), margin=20*Plots.px); display(p)

    bg_sol = solve_background(ρr0, ρm0)
    plot!(p[1], log10.(as), reduce(vcat, bg_sol.(bs; idxs=[bg.rad.Ω,bg.mat.Ω,bg.de.Ω])'); xlabel="lg(a)", ylabel="Ω", label=["Ωr" "Ωm" "ΩΛ"], legend=:left); display(p)
    plot!(p[2], log10.(as), log10.(bg_sol.(bs; idxs=bg.grav.E) / bg_sol(btoday; idxs=grav.E)); xlabel="lg(a)", ylabel="lg(H/H0)"); display(p)

    th_sol = solve_thermodynamics(ρr0, ρm0, ρb0, H0)
    plot!(p[3], log10.(as), log10.(th_sol.(bs, idxs=th.Xe)); xlabel="lg(a)", ylabel="lg(Xe)"); display(p)

    pt_sols = [solve_perturbations(kval, ρr0, ρm0, ρb0, H0) for kval in ks] # TODO: use EnsembleProblem again
    plot!(p[4], log10.(as), [pt_sol.(bs; idxs=pt.st.Φ) for pt_sol in pt_sols]; xlabel="lg(a)", ylabel="Φ"); display(p)
    plot!(p[5], log10.(as), [[log10.(abs.(pt_sol.(bs; idxs=δ))) for pt_sol in pt_sols] for δ in [pt.bar.δ,pt.cdm.δ]]; color=[(1:length(ks))' (1:length(ks))'], xlabel="lg(a)", ylabel="lg(|δb|), lg(δc)"); display(p)

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