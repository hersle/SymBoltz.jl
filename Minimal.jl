using ModelingToolkit
using DifferentialEquations
using SymbolicIndexingInterface
using DataInterpolations
using ForwardDiff, DiffResults
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
const EHe1ion = 24.58738 * eV
const EHe2ion = 54.41776 * eV

const solver = KenCarp4() # KenCarp4 and Kvaerno5 seem to work well

# TODO: shooting method https://docs.sciml.ai/DiffEqDocs/stable/tutorials/bvp_example/ (not supported by ModelingToolkit: https://github.com/SciML/ModelingToolkit.jl/issues/924, https://discourse.julialang.org/t/boundary-value-problem-with-modellingtoolkit-or-diffeqoperators/57656)
# TODO: @register_symbolic from bg -> thermo -> pert instead of reintegrating: https://docs.sciml.ai/ModelingToolkit/stable/tutorials/ode_modeling/#Specifying-a-time-variable-forcing-function
# TODO: make simpler Cosmology interface
# TODO: compare runtime for finite vs. dlgP_dlgks_autodiff
# TODO: compare accuracy with class
# TODO: non-linear: higher-order perturbations vs halofit vs N-body?
# TODO: CMB power spectrum
# TODO: baryons: Recfast -> Recfast++ -> CosmoRec -> HyRec -> HyRec-2: call out, or integrate equations into my code to make use of my background calculation?
# TODO: composable models, generate equations
# TODO: modified gravity: (coupled) quintessence, Brans-Dicke, DGP, parametrized framework, EFT of LSS, ...
# TODO: analytical solutions of e.g. background ρ evolution
# TODO: GPU-parallellized EnsembleProblem
# TODO: relate parameters through parameter expressions: https://docs.sciml.ai/ModelingToolkit/stable/basics/Composition/#Variable-scope-and-parameter-expressions
# TODO: define components with @mtkmodel?

ρr0 = 5e-5 # Ωr0
ρm0 = 0.3 # Ωm0
ρb0 = 0.02 # Ωb0
H0 = 70 * km/Mpc # s^-1
As = 2e-9
Yp = 0.25
aini, atoday = 1e-8, 1e0
as = 10 .^ range(log10(aini), log10(atoday), length=200)
k0 = 1 / 2997.92458 # h/Mpc
ks = 10 .^ range(-4, +2, length=400) / k0 # in code units of k0 = H0/c

p = plot(layout=(3,3), size=(1920, 1080), left_margin=bottom_margin=30*Plots.px); display(p) # TODO: add plot recipe!

# independent variable: scale factor
# TODO: spacetime/geometry structure?
@variables a E(a)
E = GlobalScope(E)
Da = Differential(a)

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

bg_sol = solve_background(ρr0, ρm0)
plot!(p[1], log10.(as), reduce(vcat, bg_sol.(as; idxs=[bg_sim.rad.Ω,bg_sim.mat.Ω,bg_sim.de.Ω])'); xlabel="lg(a)", ylabel="Ω", label=["Ωr" "Ωm" "ΩΛ"], legend=:left); display(p)
plot!(p[2], log10.(as), log10.(bg_sol.(as; idxs=E) / bg_sol(atoday; idxs=E)); xlabel="lg(a)", ylabel="lg(H/H0)"); display(p)

# Saha recombination (1/2)
function thermodynamics_saha(; name)
    @parameters fb H0 Yp
    @variables Xe(a) XHe₊(a) XHe₊₊(a) XH₊(a) ρb(a) nb(a) ne(a) nH(a) H(a) T(a) λe(a) R1(a) R2(a) R3(a) dτ(a)
    return ODESystem([
        Da(T) ~ -T / a # T = T0 / a # TODO: more sophisticated DE for temperature evolution

        R1 ~ 1 * exp(-EHion /(kB*T)) / (λe^3 * ne)
        R2 ~ 2 * exp(-EHe1ion/(kB*T)) / (λe^3 * ne)
        R3 ~ 4 * exp(-EHe2ion/(kB*T)) / (λe^3 * ne)
        XH₊ ~ 1 / (1 + 1/R1) # is Taylor expansion? # TODO: is this wrong? should be XH₊^2 in numerator?
        XHe₊ ~ 1 / (1 + 1/R2 + R3)
        XHe₊₊ ~ 1 / (1 + 1/R3 + 1/(R2*R3)) # TODO: XHe₊₊ ~ R3 * XHe₊
        Xe ~ XH₊ + Yp / (4*(1-Yp)) * (XHe₊ + 2*XHe₊₊) # Saha

        nb ~ ρb / mp
        nH ~ (1-Yp) * nb # TODO: correct?
        ne ~ Xe * nH
        λe ~ h / √(2π*me*kB*T) # electron de-Broglie wavelength

        dτ ~ -ne * σT * c / (a*H) # dτ = dτ/da
    ], a, [Xe, T, H, ρb, XH₊, XHe₊, XHe₊₊, dτ], [fb, H0, Yp]; name)
end

@named th1 = thermodynamics_saha()
@named th1_bg_conn = ODESystem([
    th1.H ~ E * th1.H0 # 1/s
    th1.ρb ~ th1.fb * bg.mat.ρ * 3*th1.H0^2 / (8*π*G) # kg/m³
], a)
@named th1_bg = compose(th1_bg_conn, th1, bg)
th1_sim = structural_simplify(th1_bg)
th1_prob = ODEProblem(th1_sim, unknowns(th1_sim) .=> NaN, (aini, atoday), parameters(th1_sim) .=> NaN; jac=true)

function solve_thermodynamics_saha(ρr0, ρm0, ρb0, H0, Yp)
    fb = ρb0 / ρm0; @assert fb <= 1
    Tini = (ρr0 * 15/π^2 * 3*H0^2/(8*π*G) * ħ^3*c^5)^(1/4) / kB / aini # TODO: relate to ρr0 once that is a parameter
    ρrini, ρmini, ρΛini = solve_background(ρr0, ρm0)(aini; idxs = [bg.rad.ρ, bg.mat.ρ, bg.de.ρ]) # integrate background from atoday back to aini # TODO: avoid when ρr0 etc. are parameters
    Xeini = 1 + Yp / (4*(1-Yp)) * 2
    prob = remake(th1_prob; u0 = [th1.Xe => Xeini, th1.T => Tini, bg.rad.ρ => ρrini, bg.mat.ρ => ρmini, bg.de.ρ => ρΛini], p = [th1.fb => fb, th1.H0 => H0, th1.Yp => Yp])

    Xeidx = SymbolicIndexingInterface.variable_index(prob, th1.Xe) # TODO: avoid
    condition(u, _, integrator) = u[Xeidx] > 0.99
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)

    return solve(prob, Rodas5P(), reltol=1e-5, callback = cb) # TODO: after switching ivar from a to b=ln(a), the integrator needs more steps. fix this? # TODO: terminate early
end

th1_sol = solve_thermodynamics_saha(ρr0, ρm0, ρb0, H0, Yp)
plot!(p[3], log10.(th1_sol[a]), @. log10(th1_sol(th1_sol[a], idxs=th1.Xe)); xlabel="lg(a)", ylabel="lg(Xe)", ylims=(-4,1)); display(p)

# Peebles recombination (2/2)
function thermodynamics_peebles(; name)
    @parameters fb H0 Yp
    @variables Xe(a) ρb(a) nb(a) ne(a) nH(a) H(a) T(a) dτ(a) α2(a) β(a) λe(a) C(a) Λα(a) Λ2γ(a) β2(a) dτ(a)
    return ODESystem([
        Da(T) ~ -T / a # T = T0 / a # TODO: more sophisticated DE for temperature evolution

        α2 ~ 9.78 * (α*ħ/me)^2/c * √(EHion/(kB*T)) * log(EHion/(kB*T)) # Dodelson (4.38) (e⁻ + p → H + γ)
        β ~ α2 / λe^3 * exp(-EHion/(kB*T)) # Dodelson (4.37)-(4.38) (γ + H → e⁻ + p)
        λe ~ h / √(2π*me*kB*T) # electron de-Broglie wavelength

        # Peebles' correction factor (Dodelson exercise 4.7)
        Λα ~ H * (3*EHion/(ħ*c))^3 / ((8*π)^2 * nH) # 1/s
        Λ2γ ~ 8.227 # 1/s
        β2 ~ α2 / λe^3 * exp(-EHion/(4*kB*T)) # 1/s (compute this instead of β2 = β * exp(3*EHion/(4*kB*T)) to avoid exp overflow)
        C ~ (Λ2γ + Λα) / (Λ2γ + Λα + β2)

        Da(Xe) ~ C * ((1-Xe)*β - Xe^2*nH*α2) / (a*H) # Dodelson (4.36)

        nb ~ ρb / mp
        nH ~ (1-Yp) * nb # TODO: correct?
        ne ~ Xe * nH

        # optical depth
        dτ ~ -ne * σT * c / (a*H) # dτ = dτ/da

        # TODO: reionization?
    ], a, [Xe, T, H, ρb, dτ], [fb, H0, Yp]; name)
end
@named th2 = thermodynamics_peebles()
@named th2_bg_conn = ODESystem([
    th2.H ~ E * th2.H0 # 1/s
    th2.ρb ~ th2.fb * bg.mat.ρ * 3*th2.H0^2 / (8*π*G) # kg/m³
], a)
@named th2_bg = compose(th2_bg_conn, th2, bg)
th2_sim = structural_simplify(th2_bg)
th2_prob = ODEProblem(th2_sim, unknowns(th2_sim) .=> NaN, (NaN, NaN), parameters(th2_sim) .=> NaN; jac=true)

function solve_thermodynamics_peebles(ρr0, ρm0, ρb0, H0, Yp)
    fb = ρb0 / ρm0; @assert fb <= 1
    aini, Xeini, Tini = solve_thermodynamics_saha(ρr0, ρm0, ρb0, H0, Yp)[[a, th1.Xe, th1.T]][end]
    ρrini, ρmini, ρΛini = solve_background(ρr0, ρm0)(aini; idxs = [bg.rad.ρ, bg.mat.ρ, bg.de.ρ]) # integrate background from atoday back to aini # TODO: avoid when ρr0 etc. are parameters
    prob = remake(th2_prob; tspan = (aini, atoday), u0 = [th2.Xe => Xeini, th2.T => Tini, bg.rad.ρ => ρrini, bg.mat.ρ => ρmini, bg.de.ρ => ρΛini], p = [th2.fb => fb, th2.H0 => H0, th2.Yp => Yp])
    return solve(prob, solver, reltol=1e-5) # TODO: after switching ivar from a to b=ln(a), the integrator needs more steps. fix this?
end

th2_sol = solve_thermodynamics_peebles(ρr0, ρm0, ρb0, H0, Yp)
plot!(p[3], log10.(th2_sol[a]), @. log10(th2_sol.(th2_sol[a], idxs=th2.Xe))); display(p)

function solve_thermodynamics(ρr0, ρm0, ρb0, H0, Yp)
    th1_sol = solve_thermodynamics_saha(ρr0, ρm0, ρb0, H0, Yp)
    th2_sol = solve_thermodynamics_peebles(ρr0, ρm0, ρb0, H0, Yp)
    as = [th1_sol[a][1:end-2]; th2_sol[a][:]] # TODO: why must I skip last 2 (and not 1) to avoid duplicate?
    dτs = [th1_sol[th1.dτ][1:end-2]; th2_sol[th2.dτ][:]]
    spl = CubicSpline(log.(-dτs), log.(as)) # update spline for dτ (e.g. to propagate derivative information through recombination, if called with dual numbers) TODO: use th_sol(a; idxs=th.dτ) directly in a type-stable way?
    return spl
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
@variables ρc(a) ρb(a) δρr(a) δρc(a) δρb(a) R(a) # TODO: get rid of
dτspl = nothing # to be set in solve_perturbations # TODO: use th_sol(a; idxs=th.dτ) directly in a type-stable way? # TODO: does autodiff work through this step? # TODO: define callable struct instead? like in https://discourse.julialang.org/t/registering-a-time-dependent-function-inside-a-function-in-modelingtoolkit-jl/100402/3
dτfunc(a) = -exp(dτspl(log(a))) # TODO: type-stable? @code_warntype dτfunc(1e0) gives warnings, but code seems fast?
@register_symbolic dτfunc(a)
@parameters fb # TODO: get rid of
@named pt_bg_conn = ODESystem([
    ρb ~ fb * bg.mat.ρ
    ρc ~ bg.mat.ρ - ρb
    δρr ~ rad.δ * bg.rad.ρ
    δρc ~ cdm.δ * ρc
    δρb ~ bar.δ * ρb
    grav.δρ ~ δρr + δρc + δρb # total energy density perturbation
    grav.ρm ~ bg.mat.ρ

    # baryon-photon interactions: Compton (Thomson) scattering # TODO: define connector type?
    R ~ 3/4 * ρb / bg.rad.ρ # Dodelson (5.74)
    rad.interaction ~ -dτfunc(a)/3 * (bar.u - 3*rad.Θ1)
    bar.interaction ~ +dτfunc(a)/R * (bar.u - 3*rad.Θ1)
], a)
pt_bg_conn = extend(pt_bg_conn, th1_bg_conn)

@named pt = compose(pt_bg_conn, rad, bar, cdm, grav, bg)
pt_sim = structural_simplify(pt)
pt_prob = ODEProblem(pt_sim, unknowns(pt_sim) .=> NaN, (aini, atoday), parameters(pt_sim) .=> NaN; jac=true) 

function solve_perturbations(kvals::AbstractArray, ρr0, ρm0, ρb0, H0, Yp)
    # TODO: fix instability for a few k-modes after splining Saha+Peebles dτ
    fb = ρb0 / ρm0; @assert fb <= 1 # TODO: avoid duplication thermo logic
    bg_sol = solve_background(ρr0, ρm0)
    global dτspl = solve_thermodynamics(ρr0, ρm0, ρb0, H0, Yp) # update spline for dτ (e.g. to propagate derivative information through recombination, if called with dual numbers) TODO: use th_sol(a; idxs=th.dτ) directly in a type-stable way?
    ρrini, ρmini, ρΛini, Eini = bg_sol(aini; idxs = [bg.rad.ρ, bg.mat.ρ, bg.de.ρ, E]) # integrate background from atoday back to aini
    function prob_func(_, i, _)
        kval = kvals[i]
        println("$i/$(length(kvals)) k = $(kval*k0) Mpc/h")
        Φini = 1.0 # arbitrary normalization (from primordial curvature power spectrum?)
        Θr0ini = Φini/2 # Dodelson (7.89)
        δcini = δbini = 3*Θr0ini # Dodelson (7.94)
        Θr1ini = -kval*Φini/(6*aini*Eini) # Dodelson (7.95) # TODO: replace aini -> a when this is fixed? https://github.com/SciML/ModelingToolkit.jl/issues/2543
        ucini = ubini = 3*Θr1ini # Dodelson (7.95)
        return remake(pt_prob; u0 = [Φ => Φini, pt_sim.rad.Θ0 => Θr0ini, pt_sim.rad.Θ1 => Θr1ini, pt_sim.bar.δ => δbini, pt_sim.bar.u => ubini, pt_sim.cdm.δ => δcini, pt_sim.cdm.u => ucini, bg.rad.ρ => ρrini, bg.mat.ρ => ρmini, bg.de.ρ => ρΛini], p = [pt_sim.fb => fb, k => kval])
    end
    probs = EnsembleProblem(pt_prob, prob_func = prob_func)
    sols = solve(probs, solver, EnsembleThreads(), trajectories = length(kvals), reltol=1e-8) # TODO: test GPU parallellization
    return sols
end

pt_sols = solve_perturbations(ks, ρr0, ρm0, ρb0, H0, Yp)
plot!(p[4], log10.(as), [pt_sol.(as; idxs=Φ) for pt_sol in pt_sols]; xlabel="lg(a)", ylabel="Φ/Φᵢ"); display(p)
plot!(p[5], log10.(as), [[log10.(abs.(pt_sol.(as; idxs=δ))) for pt_sol in pt_sols] for δ in [pt.bar.δ,pt.cdm.δ]]; color=[(1:length(ks))' (1:length(ks))'], xlabel="lg(a)", ylabel="lg(|δb|), lg(δc)"); display(p)

# power spectra
θ0 = [ρr0, ρm0, ρb0, H0, As, Yp]
P0(k, As) = @. As / k ^ 3
P(k, ρr0, ρm0, ρb0, H0, As, Yp) = P0(k, As) .* solve_perturbations(k, ρr0, ρm0, ρb0, H0, Yp)(atoday; idxs=pt.grav.Δm) .^ 2
P(k, θ) = P(k, θ...) # unpack parameters θ = [ρr0, ρm0, ρb0, H0, As]

function plot_dlgP_dθs(dlgP_dθs, name, color)
    plot!(p[7], log10.(ks*k0), dlgP_dθs[:,1]; xlabel="lg(k/(h/Mpc))", ylabel="d lg(P) / d lg(Ωr0)", color=color, label=name, ylims=(-2, 0)); display(p)
    plot!(p[8], log10.(ks*k0), dlgP_dθs[:,2]; xlabel="lg(k/(h/Mpc))", ylabel="d lg(P) / d lg(Ωm0)", color=color, label=name, ylims=(-2, 3)); display(p)
    plot!(p[9], log10.(ks*k0), dlgP_dθs[:,4]; xlabel="lg(k/(h/Mpc))", ylabel="d lg(P) / d lg(H0)", color=color, label=name); display(p)
end

# computer power spectrum and derivatives wrt. input parameters using autodiff in one go
Pres = DiffResults.JacobianResult(ks, θ0)
ForwardDiff.jacobian!(Pres, θ -> log10.(P(ks, 10 .^ θ)/k0^3), log10.(θ0))
lgPs, dlgP_dθs_ad = DiffResults.value(Pres), DiffResults.jacobian(Pres)
plot!(p[6], log10.(ks*k0), lgPs; xlabel="lg(k/(h/Mpc))", label="lg(P/(Mpc/h)³)", color=3, legend=:bottomleft); display(p)
plot_dlgP_dθs(dlgP_dθs_ad, "auto. diff.", 1)

# compute derivatives of power spectrum using finite differences
dlgP_dθs_fd = FiniteDiff.finite_difference_jacobian(θ -> log10.(P(ks, 10 .^ θ)/k0^3), log10.(θ0); relstep=1e-4) # relstep is important for finite difference accuracy!
plot_dlgP_dθs(dlgP_dθs_fd, "fin. diff.", 2)