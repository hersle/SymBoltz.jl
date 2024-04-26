using ModelingToolkit
using DifferentialEquations
using DataInterpolations
using ForwardDiff, DiffResults, FiniteDiff
using OffsetArrays
using Plots; Plots.default(label=nothing, markershape=:pixel)

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
const σT = PhysicalConstants.CODATA2018.σ_e / u"m^2"
const km  = 1u"km/m"  |> NoUnits
const Mpc = 1u"Mpc/m" |> NoUnits
const EHion = 13.59844u"eV/J" |> NoUnits
const EHe1ion = 24.58738u"eV/J" |> NoUnits
const EHe2ion = 54.41776u"eV/J" |> NoUnits

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
# TODO: GPU-parallellized EnsembleProblem
# TODO: relate parameters through parameter expressions: https://docs.sciml.ai/ModelingToolkit/stable/basics/Composition/#Variable-scope-and-parameter-expressions
# TODO: define components with @mtkmodel?

@kwdef struct Parameters
    ρr0 = 5e-5 # Ωr0
    ρm0 = 0.3 # Ωm0
    ρb0 = 0.02 # Ωb0
    H0 = 70 * km/Mpc # s^-1
    As = 2e-9
    Yp = 0.25
end
par = Parameters()
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
        Da(ρ) ~ -3/a * (ρ + P) # alternative analytical solution: ρ ~ ρ0 / a^(3*(1+w))
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

@named bg = compose(bg, rad, mat, de, grav)
bg_sim = structural_simplify(bg)
bg_prob = ODEProblem(bg_sim, unknowns(bg_sim) .=> NaN, (atoday, aini), parameters(bg_sim) .=> NaN)
function solve_background(ρr0, ρm0)
    ρΛ0 = 1 - ρr0 - ρm0 # TODO: handle with equation between parameters once ρr0 etc. are parameters?
    prob = remake(bg_prob; u0 = [bg_sim.rad.ρ => ρr0, bg_sim.mat.ρ => ρm0, bg_sim.de.ρ => ρΛ0])
    return solve(prob, Tsit5(), reltol=1e-7) # using KenCarp4 here leads to difference in AD vs FD
end
solve_background(θ::Parameters) = solve_background(θ.ρr0, θ.ρm0)

bg_sol = solve_background(par.ρr0, par.ρm0)
plot!(p[1], log10.(as), reduce(vcat, bg_sol.(as; idxs=[bg_sim.rad.Ω,bg_sim.mat.Ω,bg_sim.de.Ω])'); xlabel="lg(a)", ylabel="Ω", label=["Ωr" "Ωm" "ΩΛ"], legend=:left); display(p)
plot!(p[2], log10.(as), log10.(bg_sol.(as; idxs=E) / bg_sol(atoday; idxs=E)); xlabel="lg(a)", ylabel="lg(H/H0)"); display(p)

# background thermodynamcis / recombination
Hifelse(x, v1, v2; k=1) = 1/2 * ((v1+v2) + (v2-v1)*tanh(k*x)) # smooth transition/step function from v1 at x<0 to v2 at x>0
function recombination_helium_saha(; name, Xeconst=1e-5)
    @parameters Yp
    @variables Xe(a) XH₊(a) XHe₊(a) XHe₊₊(a) ne(a) nH(a) T(a) λ(a) R1(a) R2(a) R3(a)
    return ODESystem([
        λ ~ h / √(2π*me*kB*T) # e⁻ de-Broglie wavelength
        ne ~ Xe * nH
        R1 ~ 1 * exp(-EHion  /(kB*T)) / (λ^3 * ne) # TODO: simplify equations?
        R2 ~ 2 * exp(-EHe1ion/(kB*T)) / (λ^3 * ne)
        R3 ~ 4 * exp(-EHe2ion/(kB*T)) / (λ^3 * ne)
        XH₊ ~ R1 / (1 + R1)
        XHe₊ ~ R2 / (1 + R2 + R2*R3)
        XHe₊₊ ~ R3 * XHe₊
        Xe ~ XH₊ + Yp / (4*(1-Yp)) * (XHe₊ + 2*XHe₊₊) + Xeconst # add small constant Xe which is greater than integrator tolerance to avoid solver giving tiny negative values
    ], a, [Xe, XH₊, XHe₊, XHe₊₊, T, nH, ne], [Yp]; name)
end

function recombination_hydrogen_peebles(; name)
    @variables Xe(a) ne(a) nH(a) T(a) H(a) λ(a) α2(a) β(a) C(a) Λα(a) Λ2γ(a) β2(a)
    return ODESystem([
        ne ~ Xe * nH
        λ ~ h / √(2π*me*kB*T) # e⁻ de-Broglie wavelength
        α2 ~ 9.78 * (α*ħ/me)^2/c * √(EHion/(kB*T)) * log(EHion/(kB*T)) # Dodelson (4.38) (e⁻ + p → H + γ) # TODO: add Recfast fudge factor?
        β  ~ α2 / λ^3 * exp(-EHion/(  kB*T)) # Dodelson (4.37)-(4.38) (γ + H → e⁻ + p)
        β2 ~ α2 / λ^3 * exp(-EHion/(4*kB*T)) # 1/s (compute this instead of β2 = β * exp(3*EHion/(4*kB*T)) to avoid exp overflow)
        Λα ~ H * (3*EHion/(ħ*c))^3 / ((8*π)^2 * nH) # 1/s
        Λ2γ ~ 8.227 # 1/s
        C ~ (Λ2γ + Λα) / (Λ2γ + Λα + β2) # Peebles' correction factor (Dodelson exercise 4.7)
        Da(Xe) ~ C * ((1-Xe)*β - Xe^2*nH*α2) / (a*H) # remains ≈ 0 during Saha recombinations, so no need to manually turn off
    ], a, [Xe, H, nH, ne, T], []; name)
end

function thermodynamics_temperature(; name)
    @variables Tγ(a) Tb(a) ργ(a) ρb(a) dτ(a) fγb(a)
    return ODESystem([
        Da(Tγ) ~   -Tγ/a # Tγ = Tγ0 / a
        Da(Tb) ~ -2*Tb/a - 8/3*(mp/me)*fγb*a*dτ*(Tγ-Tb) # TODO: multiply last term by a or not?
    ], a, [Tγ, Tb, dτ, fγb], []; name)
end

function reionization_smooth_step(; name)
    y(z) = (1+z)^(3/2)
    Δy(z, Δz0) = 3/2 * (1+z)^(1/2) * Δz0
    @parameters z0 Δz0 Xe0
    @variables z(a)
    return ODESystem([
        z ~ 1/a - 1
        Xe ~ Hifelse(y(z0)-y(z), 0, Xe0; k=1/Δy(z0, Δz0)) # smooth step from 0 to Xe0
    ], a, [Xe], [z0, Δz0, Xe0]; name)
end

@parameters fb H0 Yp
@variables Xe(a) ne(a) dτ(a) H(a) ρb(a) nb(a) nH(a)
@named temp = thermodynamics_temperature()
@named saha = recombination_helium_saha()
@named peebles = recombination_hydrogen_peebles()
@named reion1 = reionization_smooth_step() # TODO: separate reionH₊, reionHe₊, reionHe₊₊
@named reion2 = reionization_smooth_step()
@named th_bg_conn = ODESystem([
    H ~ E * H0 # 1/s
    ρb ~ fb * bg.mat.ρ * 3*H0^2 / (8*π*G) # kg/m³
    nb ~ ρb / mp # 1/m³
    nH ~ (1-Yp) * nb # TODO: correct?

    temp.fγb ~ bg.rad.ρ / (fb*bg.mat.ρ) # ργ/ρb
    temp.dτ ~ dτ
    saha.T ~ temp.Tb
    saha.nH ~ nH
    peebles.T ~ temp.Tb
    peebles.nH ~ nH
    peebles.H ~ H
    
    # switch *smoothly* from Saha to Peebles when XeS ≤ 1 (see e.g. https://discourse.julialang.org/t/handling-instability-when-solving-ode-problems/9019/5) # TODO: make into a connection
    Xe ~ Hifelse(1-saha.Xe, saha.Xe, peebles.Xe; k=1e3) + reion1.Xe + reion2.Xe
    ne ~ Xe * nH
    dτ ~ -ne * σT * c / (a*H) # common optical depth dτ = dτ/da # TODO: separate in Saha/Peebles?
], a)
@named th = compose(th_bg_conn, saha, peebles, temp, reion1, reion2, bg)
th_sim = structural_simplify(th)
th_prob = ODEProblem(th_sim, unknowns(th_sim) .=> NaN, (aini, atoday), parameters(th_sim) .=> NaN; jac=true)

function solve_thermodynamics(ρr0, ρm0, ρb0, H0, Yp)
    ρrini, ρmini, ρΛini = solve_background(ρr0, ρm0)(aini; idxs=[bg.rad.ρ, bg.mat.ρ, bg.de.ρ]) # TODO: avoid duplicate logic
    fb = ρb0 / ρm0; @assert fb <= 1
    Tini = (ρr0 * 15/π^2 * 3*H0^2/(8*π*G) * ħ^3*c^5)^(1/4) / kB / aini # common initial Tb = Tγ TODO: relate to ρr0 once that is a parameter
    XeSini = 1 + Yp / (4*(1-Yp)) * 2 # TODO: avoid?
    prob = remake(th_prob; u0 = [saha.Xe => XeSini, peebles.Xe => 1.0, temp.Tγ => Tini, temp.Tb => Tini, bg.rad.ρ => ρrini, bg.mat.ρ => ρmini, bg.de.ρ => ρΛini], p = [th_sim.fb => fb, th_sim.H0 => H0, th_sim.Yp => Yp, saha.Yp => Yp, reion1.z0 => 8, reion1.Δz0 => 0.5, reion1.Xe0 => 1+Yp/(4*(1-Yp)), reion2.z0 => 3.5, reion2.Δz0 => 0.5, reion2.Xe0 => Yp/(4*(1-Yp))])
    return solve(prob, RadauIIA5(), reltol=1e-7) # CLASS uses "NDF15" (https://lesgourg.github.io/class-tour/London2014/Numerical_Methods_in_CLASS_London.pdf) TODO: after switching ivar from a to b=ln(a), the integrator needs more steps. fix this?
end
solve_thermodynamics(θ::Parameters) = solve_thermodynamics(θ.ρr0, θ.ρm0, θ.ρb0, θ.H0, θ.Yp)

th_sol = solve_thermodynamics(par.ρr0, par.ρm0, par.ρb0, par.H0, par.Yp)
plot!(p[3], log10.(th_sol[a]), stack(th_sol(th_sol[a], idxs=[saha.Xe, peebles.Xe, reion1.Xe, reion2.Xe, th_sim.Xe]))'; xlabel="lg(a)", ylabel="Xe", ylims=(0, 1.5), label=["XeS" "XeP" "XeRE1" "XeRE2" "Xe"], legend=:bottomleft); display(p)
#plot!(p[4], log10.(th_sol[a]), log10.(th_sol[th.Tγ])); display(p)
#plot!(p[4], log10.(th_sol[a]), log10.(th_sol[th.Tb])); display(p)

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

function perturbations_photon_hierarchy(lmax=6, interact=false; name)
    @variables Θ(a)[0:lmax] δ(a) interactions(a)[1:lmax-1]
    eqs = [
        Da(Θ[0]) + k/(a^2*E)*Θ[1] ~ -Da(Φ)
        Da(Θ[1]) - k/(3*a^2*E)*(Θ[0]-2*Θ[2]) ~ k/(3*a^2*E)*Ψ + interactions[1]
        [Da(Θ[l]) ~ k/(a^2*E)/(2*l+1) * (l*Θ[l-1] - (l+1)*Θ[l+1]) + interactions[l] for l in 2:lmax-1]...
        Θ[lmax] ~ 0 # TODO: integrate η(a) and use better cutoff
        δ ~ 4*Θ[0]
    ]
    if !interact
        push!(eqs, collect(interactions .~ 0)...)
    end
    return ODESystem(eqs, a; name)
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
    @variables δρ(a) Δm(a) ρm(a) Π(a)
    return ODESystem([
        Da(Φ) ~ (3/2*a^2*δρ - k^2*Φ - 3*(a*E)^2*Φ) / (3*a^3*E^2) # Dodelson (8.14) # TODO: write in more natural form?
        k^2 * (Ψ + Φ) ~ Π # anisotropic stress
        Δm ~ k^2*Φ / (3/2*a^2*ρm) # gauge-invariant overdensity (from Poisson equation) # TODO: move outside or change?
    ], a; name)
end

lmax = 6
@named rad = perturbations_photon_hierarchy(lmax, true)
@named cdm = perturbations_matter(false)
@named bar = perturbations_matter(true)
@named grav = perturbations_gravity()
@parameters fb dτspline # TODO: get rid of
@variables ρc(a) ρb(a) δρr(a) δρc(a) δρb(a) R(a) dτ(a) # TODO: get rid of
dτfunc(a, spl) = -exp(spl(log(a))) # TODO: type-stable? @code_warntype dτfunc(1e0) gives warnings, but code seems fast?
@register_symbolic dτfunc(a, spl)
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
    bar.interaction     ~ +dτ/R * (bar.u - 3*rad.Θ[1])
    rad.interactions[1] ~ -dτ/3 * (bar.u - 3*rad.Θ[1])
    rad.interactions[2] ~  dτ*rad.Θ[2]*9/10 # TODO: add polarization
    [rad.interactions[l] ~ dτ*rad.Θ[l] for l in 3:lastindex(rad.interactions)]...
    dτ ~ dτfunc(a, dτspline)

    # gravity shear stress
    grav.Π ~ -12*a^2 * bg.rad.ρ*rad.Θ[2] # TODO: add neutrinos
], a)
@named pt = compose(pt_bg_conn, rad, bar, cdm, grav, bg)
pt_sim = structural_simplify(pt)
pt_prob = ODEProblem(pt_sim, unknowns(pt_sim) .=> NaN, (aini, atoday), parameters(pt_sim) .=> NaN; jac=true) 

function solve_perturbations(kvals::AbstractArray, ρr0, ρm0, ρb0, H0, Yp)
    fb = ρb0 / ρm0; @assert fb <= 1 # TODO: avoid duplication thermo logic
    th_sol = solve_thermodynamics(ρr0, ρm0, ρb0, H0, Yp) # update spline for dτ (e.g. to propagate derivative information through recombination, if called with dual numbers) TODO: use th_sol(a; idxs=th.dτ) directly in a type-stable way?
    dτspline = CubicSpline(log.(-th_sol[th.dτ]), log.(th_sol[a])) # update spline for dτ (e.g. to propagate derivative information through recombination, if called with dual numbers) TODO: use th_sol(a; idxs=th.dτ) directly in a type-stable way?
    ρrini, ρmini, ρΛini, Eini = th_sol(aini; idxs = [bg.rad.ρ, bg.mat.ρ, bg.de.ρ, E]) # integrate background from atoday back to aini
    dτini = dτfunc(aini, dτspline)
    function prob_func(_, i, _)
        kval = kvals[i]
        println("$i/$(length(kvals)) k = $(kval*k0) Mpc/h")
        Φini = 1.0 # arbitrary normalization (from primordial curvature power spectrum?)
        lmax = lastindex(pt_sim.rad.Θ)
        Θrini = OffsetVector(Vector{Any}(undef, lmax), -1) # index from l=0 to l=lmax-1
        Θrini[0] = Φini/2 # Dodelson (7.89)
        Θrini[1] = -kval*Φini/(6*aini*Eini) # Dodelson (7.95)
        Θrini[2] = -20/45*kval/(aini^2*Eini*dτini) # TODO: change with polarization
        for l in 2:lmax-1
            Θrini[l] = -l/(2*l+1) * kval/(aini^2*Eini*dτini) * Θrini[l-1]
        end
        δcini = δbini = 3*Θrini[0] # Dodelson (7.94)
        ucini = ubini = 3*Θrini[1] # Dodelson (7.95)
        return remake(pt_prob; u0 = Dict(Φ => Φini, [pt_sim.rad.Θ[l] => Θrini[l] for l in 0:lmax-1]..., pt_sim.bar.δ => δbini, pt_sim.bar.u => ubini, pt_sim.cdm.δ => δcini, pt_sim.cdm.u => ucini, bg.rad.ρ => ρrini, bg.mat.ρ => ρmini, bg.de.ρ => ρΛini), p = [pt_sim.fb => fb, k => kval, pt_sim.dτspline => dτspline])
    end
    probs = EnsembleProblem(prob = nothing, prob_func = prob_func)
    return solve(probs, KenCarp4(), EnsembleThreads(), reltol=1e-8, trajectories = length(kvals)) # KenCarp4 and Kvaerno5 seem to work well # TODO: test GPU parallellization
end
solve_perturbations(kvals::AbstractArray, θ::Parameters) = solve_perturbations(kvals, θ.ρr0, θ.ρm0, θ.ρb0, θ.H0, θ.Yp)

pt_sols = solve_perturbations(ks, par.ρr0, par.ρm0, par.ρb0, par.H0, par.Yp)
plot!(p[4], log10.(as), [pt_sol.(as; idxs=Φ) for pt_sol in pt_sols]; xlabel="lg(a)", ylabel="Φ/Φᵢ"); display(p)
plot!(p[5], log10.(as), [[log10.(abs.(pt_sol.(as; idxs=δ))) for pt_sol in pt_sols] for δ in [pt.bar.δ,pt.cdm.δ]]; color=[(1:length(ks))' (1:length(ks))'], xlabel="lg(a)", ylabel="lg(|δb|), lg(δc)"); display(p)

# power spectra
θ0 = [par.ρr0, par.ρm0, par.ρb0, par.H0, par.As, par.Yp]
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
log10Ph3(log10θ) = log10.(P(ks, 10 .^ log10θ)/k0^3)
ForwardDiff.jacobian!(Pres, log10Ph3, log10.(θ0))
lgPs, dlgP_dθs_ad = DiffResults.value(Pres), DiffResults.jacobian(Pres)
plot!(p[6], log10.(ks*k0), lgPs; xlabel="lg(k/(h/Mpc))", label="lg(P/(Mpc/h)³)", color=3, legend=:bottomleft); display(p)
plot_dlgP_dθs(dlgP_dθs_ad, "auto. diff.", 1)

# compute derivatives of power spectrum using finite differences
dlgP_dθs_fd = FiniteDiff.finite_difference_jacobian(log10Ph3, log10.(θ0); relstep=1e-4) # relstep is important for finite difference accuracy!
plot_dlgP_dθs(dlgP_dθs_fd, "fin. diff.", 2)