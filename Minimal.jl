using ModelingToolkit
using DifferentialEquations
using DataInterpolations
using ForwardDiff, DiffResults, FiniteDiff
using Bessels: sphericalbesselj
using Trapz
using Roots
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

δkron(i, j) = (i == j ? 1 : 0)

# TODO: shooting method https://docs.sciml.ai/DiffEqDocs/stable/tutorials/bvp_example/ (not supported by ModelingToolkit: https://github.com/SciML/ModelingToolkit.jl/issues/924, https://discourse.julialang.org/t/boundary-value-problem-with-modellingtoolkit-or-diffeqoperators/57656)
# TODO: @register_symbolic from bg -> thermo -> pert instead of reintegrating: https://docs.sciml.ai/ModelingToolkit/stable/tutorials/ode_modeling/#Specifying-a-time-variable-forcing-function
# TODO: make simpler Cosmology interface
# TODO: compare runtime for finite vs. dlgP_dlgks_autodiff
# TODO: compare accuracy with class
# TODO: non-linear: higher-order perturbations vs halofit vs N-body?
# TODO: baryons: Recfast -> Recfast++ -> CosmoRec -> HyRec -> HyRec-2: call out, or integrate equations into my code to make use of my background calculation?
# TODO: composable models, generate equations
# TODO: modified gravity: (coupled) quintessence, Brans-Dicke, DGP, parametrized framework, EFT of LSS, ...
# TODO: GPU-parallellized EnsembleProblem
# TODO: relate parameters through parameter expressions: https://docs.sciml.ai/ModelingToolkit/stable/basics/Composition/#Variable-scope-and-parameter-expressions
# TODO: define components with @mtkmodel?

@kwdef struct Parameters
    Ωr0 = 5.5e-5 # Ωr0
    Ωm0 = 0.317 # Ωm0
    Ωb0 = 0.05 # Ωb0
    H0 = 67 * km/Mpc # s^-1
    As = 2.1e-9
    Yp = 0.245
end
par = Parameters()
aini, atoday = 1e-8, 1e0
as = 10 .^ range(log10(aini), log10(atoday), length=200)
k0 = 1 / 2997.92458 # h/Mpc
ks = 10 .^ range(-4, +2, length=400) / k0 # in code units of k0 = H0/c

p = plot(layout=(3,3), size=(1920, 1080), left_margin=bottom_margin=30*Plots.px); display(p) # TODO: add plot recipe!

# independent variable: scale factor
# TODO: spacetime/geometry structure?
@variables η a(η) E(η) # η is really η in units of 1/H0
a = GlobalScope(a)
Dη = Differential(η)

function background_gravity_GR(; name)
    @variables ρ(η)
    return ODESystem([
        Dη(a) ~ √(ρ)* a^2 # TODO: 8π/3 factor?
    ], η; name)
end

function background_species_constant_eos(w; name)
    @variables ρ(η) P(η) Ω(η) ρcrit(η)
    return ODESystem([
        P ~ w*ρ
        Dη(ρ) ~ -3 * Dη(a)/a * (ρ + P) # alternative analytical solution: ρ ~ ρ0 / a^(3*(1+w))
        Ω ~ ρ / ρcrit
    ], η; name)
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
    E ~ Dη(a) / a^2
], η, [a, E], [])

@named bg = compose(bg, rad, mat, de, grav)
bg_sim = structural_simplify(bg)
bg_prob = ODEProblem(bg_sim, unknowns(bg_sim) .=> NaN, (0.0, 4.0), parameters(bg_sim) .=> NaN; jac=true)
function solve_background(Ωr0, Ωm0)
    ΩΛ0 = 1 - Ωr0 - Ωm0 # TODO: handle with equation between parameters once ρr0 etc. are parameters?
    ηini = aini / √(Ωr0) # analytical radiation-dominated solution
    Ωrini = Ωr0 / aini^4 # TODO: avoid!
    Ωmini = Ωm0 / aini^3
    ΩΛini = ΩΛ0
    prob = remake(bg_prob; tspan=(ηini, 4.0), u0 = [bg_sim.a => aini, bg_sim.rad.ρ => Ωrini, bg_sim.mat.ρ => Ωmini, bg_sim.de.ρ => ΩΛini])
    # TODO: stop when a == 1, see https://github.com/hersle/SymBoltz.jl/commit/7c4beb56eb13a1bb4dc5b30f8f808ca561989e4d
    return solve(prob, KenCarp4(), reltol=1e-8) # using KenCarp4 here leads to difference in AD vs FD
end
solve_background(θ::Parameters) = solve_background(θ.Ωr0, θ.Ωm0)
ηi(bg_sol::ODESolution) = bg_sol.prob.tspan[1]
η0(bg_sol::ODESolution) = find_zero(η -> bg_sol(η, idxs=bg.a) - 1.0, bg_sol.prob.tspan)

bg_sol = solve_background(par.Ωr0, par.Ωm0)
println("ηi = $(ηi(bg_sol)), η0 = $(η0(bg_sol))")
plot!(p[1], bg_sol[η], bg_sol[a]; xlabel="η / (1/H0)", ylabel="a", ylims=(0, 1)); display(p)
plot!(p[2], log10.(bg_sol[a]), stack(bg_sol[[bg_sim.rad.Ω,bg_sim.mat.Ω,bg_sim.de.Ω]])'; xlabel="lg(a)", ylabel="Ω", label=["Ωr" "Ωm" "ΩΛ"], legend=:left); display(p)

# background thermodynamcis / recombination
Hifelse(x, v1, v2; k=1) = 1/2 * ((v1+v2) + (v2-v1)*tanh(k*x)) # smooth transition/step function from v1 at x<0 to v2 at x>0
function recombination_helium_saha(; name, Xeconst=1e-5)
    @parameters Yp
    @variables Xe(η) XH₊(η) XHe₊(η) XHe₊₊(η) ne(η) nH(η) T(η) λ(η) R1(η) R2(η) R3(η)
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
    ], η, [Xe, XH₊, XHe₊, XHe₊₊, T, nH, ne], [Yp]; name)
end

function recombination_hydrogen_peebles(; name)
    @parameters H0
    @variables Xe(η) nH(η) T(η) H(η) λ(η) α2(η) β(η) C(η) Λ2γ_Λα(η) β2_Λα(η) actiweight(η)
    return ODESystem([
        λ ~ h / √(2π*me*kB*T) # e⁻ de-Broglie wavelength
        α2 ~ 9.78 * (α*ħ/me)^2/c * √(EHion/(kB*T)) * log(EHion/(kB*T)) # Dodelson (4.38) (e⁻ + p → H + γ) # TODO: add Recfast fudge factor?
        β  ~ α2 / λ^3 * exp(-EHion/(kB*T)) # Dodelson (4.37)-(4.38) (γ + H → e⁻ + p)
        β2_Λα ~ α2 / λ^3 * exp(-EHion/(4*kB*T)) * ((8*π)^2 * (1-Xe) * nH) / (H * (3*EHion/(ħ*c))^3) # β2/Λα (compute this instead of β2 = β * exp(3*EHion/(4*kB*T)) to avoid exp overflow)
        Λ2γ_Λα ~ 8.227 * ((8*π)^2 * (1-Xe) * nH) / (H * (3*EHion/(ħ*c))^3) # Λ2γ/Λα
        C ~ Hifelse(actiweight, 1, (1 + Λ2γ_Λα) / (1 + Λ2γ_Λα + β2_Λα); k=1e3) # Peebles' correction factor (Dodelson exercise 4.7), manually activated to avoid numerical issues at early times # TODO: activate using internal quantities only! # TODO: why doesnt it work to activate from 0? activating from 1 is really unnatural
        Dη(Xe) ~ C * ((1-Xe)*β - Xe^2*nH*α2) * a / H0 # remains ≈ 0 during Saha recombinations, so no need to manually turn off (multiply by H0 on left because cide η is physical η/(1/H0))
    ], η, [Xe, H, nH, T, actiweight], [H0]; name)
end

function thermodynamics_temperature(; name)
    @variables Tγ(η) Tb(η) ργ(η) ρb(η) τ(η) fγb(η)
    return ODESystem([
        Dη(Tγ) ~ -1*Tγ * Dη(a)/a # Tγ = Tγ0 / a # TODO: introduce ℋ = Dη(a) / a?
        Dη(Tb) ~ -2*Tb * Dη(a)/a - 8/3*(mp/me)*fγb*a*Dη(τ)*(Tγ-Tb) # TODO: multiply last term by a or not?
    ], η, [Tγ, Tb, τ, fγb], []; name)
end

function reionization_smooth_step(; name)
    y(z) = (1+z)^(3/2)
    Δy(z, Δz0) = 3/2 * (1+z)^(1/2) * Δz0
    @parameters z0 Δz0 Xe0
    @variables z(η)
    return ODESystem([
        z ~ 1/a - 1
        Xe ~ Hifelse(y(z0)-y(z), 0, Xe0; k=1/Δy(z0, Δz0)) # smooth step from 0 to Xe0
    ], η, [Xe], [z0, Δz0, Xe0]; name)
end

@parameters fb H0 Yp
@variables Xe(η) ne(η) τ(η) H(η) ρb(η) nb(η) nH(η)
@named temp = thermodynamics_temperature()
@named saha = recombination_helium_saha()
@named peebles = recombination_hydrogen_peebles()
@named reion1 = reionization_smooth_step() # TODO: separate reionH₊, reionHe₊, reionHe₊₊
@named reion2 = reionization_smooth_step()
@named th_bg_conn = ODESystem([
    H ~ Dη(a) / a^2 * H0 # 1/s, H ~ ℋ / a
    ρb ~ fb * bg.mat.ρ * 3*H0^2 / (8*π*G) # kg/m³
    nb ~ ρb / mp # 1/m³
    nH ~ (1-Yp) * nb # TODO: correct?

    temp.fγb ~ bg.rad.ρ / (fb*bg.mat.ρ) # ργ/ρb
    temp.τ ~ τ
    saha.T ~ temp.Tb
    saha.nH ~ nH
    peebles.T ~ temp.Tb
    peebles.nH ~ nH
    peebles.H ~ H
    peebles.actiweight ~ 1 - saha.Xe
    
    # switch *smoothly* from Saha to Peebles when XeS ≤ 1 (see e.g. https://discourse.julialang.org/t/handling-instability-when-solving-ode-problems/9019/5) # TODO: make into a connection
    Xe ~ Hifelse(peebles.actiweight, saha.Xe, peebles.Xe; k=1e3) + reion1.Xe + reion2.Xe
    ne ~ Xe * nH
    Dη(τ) * H0 ~ -ne * σT * c * a # common optical depth τ (multiply by H0 on left because code η is physical η/(1/H0)) # TODO: separate in Saha/Peebles?
], η)
@named th = compose(th_bg_conn, saha, peebles, temp, reion1, reion2, bg)
th_sim = structural_simplify(th)
th_prob = ODEProblem(th_sim, unknowns(th_sim) .=> NaN, (0.0, 4.0), parameters(th_sim) .=> NaN; jac=true)

function solve_thermodynamics(Ωr0, Ωm0, Ωb0, H0, Yp)
    bg_sol = solve_background(Ωr0, Ωm0)
    ηini, ηtoday = ηi(bg_sol), η0(bg_sol)
    Ωrini, Ωmini, ΩΛini = bg_sol(ηini; idxs=[bg.rad.ρ, bg.mat.ρ, bg.de.ρ]) # TODO: avoid duplicate logic
    fb = Ωb0 / Ωm0; @assert fb <= 1
    Tini = (Ωrini * 15/π^2 * 3*H0^2/(8*π*G) * ħ^3*c^5)^(1/4) / kB # common initial Tb = Tγ TODO: relate to ρr0 once that is a parameter
    XeSini = 1 + Yp / (4*(1-Yp)) * 2 # TODO: avoid?
    prob = remake(th_prob; tspan = (ηini, ηtoday), u0 = [saha.Xe => XeSini, peebles.Xe => 1.0, temp.Tγ => Tini, temp.Tb => Tini, th_sim.τ => 0.0, bg.rad.ρ => Ωrini, bg.mat.ρ => Ωmini, bg.de.ρ => ΩΛini, bg.a => aini], p = [th_sim.fb => fb, th_sim.H0 => H0, th_sim.peebles.H0 => H0, th_sim.Yp => Yp, saha.Yp => Yp, reion1.z0 => 8, reion1.Δz0 => 0.5, reion1.Xe0 => 1+Yp/(4*(1-Yp)), reion2.z0 => 3.5, reion2.Δz0 => 0.5, reion2.Xe0 => Yp/(4*(1-Yp))])
    return solve(prob, RadauIIA5(), reltol=1e-8) # CLASS uses "NDF15" (https://lesgourg.github.io/class-tour/London2014/Numerical_Methods_in_CLASS_London.pdf) TODO: after switching ivar from a to b=ln(a), the integrator needs more steps. fix this?
end
solve_thermodynamics(θ::Parameters) = solve_thermodynamics(θ.Ωr0, θ.Ωm0, θ.Ωb0, θ.H0, θ.Yp)

th_sol = solve_thermodynamics(par.Ωr0, par.Ωm0, par.Ωb0, par.H0, par.Yp)
plot!(p[3], log10.(th_sol[a]), log10.(abs.(stack(th_sol[[saha.Xe, peebles.Xe, reion1.Xe, reion2.Xe, th_sim.Xe]])')); xlabel="lg(a)", ylabel="lg(Xe)", ylims=(-5, +1), label=["XeS" "XeP" "XeRE1" "XeRE2" "Xe"], legend=:bottomleft); display(p)
#plot!(p[4], log10.(th_sol[a]), log10.(stack(th_sol[[th.temp.Tγ, th.temp.Tb]])'); xlabel = "lg(a)", ylabel = "lg(T/K)", labels = ["Tγ" "Tb"]); display(p)

@variables Φ(η) Ψ(η)
@parameters k # perturbation wavenumber # TODO: associate like pt.k
k, Φ, Ψ = GlobalScope.([k, Φ, Ψ])

function perturbations_radiation(interact=false; name)
    @variables Θ0(η) Θ1(η) δ(η)
    interaction = interact ? only(@variables interaction(η)) : 0
    return ODESystem([
        Dη(Θ0) + k*Θ1 ~ -Dη(Φ) # Dodelson (5.67) or (8.10)
        Dη(Θ1) - k/3*Θ0 ~ k/3*Ψ + interaction # Dodelson (5.67) or (8.11)
        δ ~ 4*Θ0
    ], η; name)
end

function perturbations_photon_hierarchy(lmax=6, interact=false; name)
    @variables Θ0(η) Θ(η)[1:lmax] δ(η) interactions(η)[1:lmax]
    eqs = [
        Dη(Θ0) + k*Θ[1] ~ -Dη(Φ)
        Dη(Θ[1]) - k/3*(Θ0-2*Θ[2]) ~ k/3*Ψ + interactions[1]
        [Dη(Θ[l]) ~ k/(2*l+1) * (l*Θ[l-1] - (l+1)*Θ[l+1]) + interactions[l] for l in 2:lmax-1]...
        Dη(Θ[lmax]) ~ k*Θ[lmax-1] - (lmax+1) * Θ[lmax] / η + interactions[lmax]
        δ ~ 4*Θ0
    ]
    if !interact
        push!(eqs, collect(interactions .~ 0)...)
    end
    return ODESystem(eqs, η; name)
end

# TODO: merge with photon_hierarchy and have polarization flag
function perturbations_polarization_hierarchy(lmax=6; name)
    @variables Θ0(η) Θ(η)[1:lmax] dτ(η) Π(η) # TODO: index Θ[l=0] when fixed: https://github.com/SciML/ModelingToolkit.jl/pull/2671
    eqs = [
        Dη(Θ0) + k*Θ[1] ~ dτ * (Θ0 - Π/2)
        Dη(Θ[1]) - k/(2*1+1) * (1*Θ0 - (1+1)*Θ[1+1]) ~ dτ * (Θ[1] - Π/10*δkron(1,2))
        [Dη(Θ[l]) - k/(2*l+1) * (l*Θ[l-1] - (l+1)*Θ[l+1]) ~ dτ * (Θ[l] - Π/10*δkron(l,2)) for l in 2:lmax-1]...
        Dη(Θ[lmax]) ~ k*Θ[lmax-1] - (lmax+1) * Θ[lmax] / η + dτ * Θ[lmax]
    ]
    return ODESystem(eqs, η; name)
end

function perturbations_matter(interact=false; name)
    @variables δ(η) u(η)
    interaction = interact ? only(@variables interaction(η)) : 0
    return ODESystem([
        Dη(δ) + k*u ~ -3*Dη(Φ) # Dodelson (5.69) or (8.12) with i*uc -> uc
        Dη(u) + u*Dη(a)/a ~ k*Ψ + interaction # Dodelson (5.70) or (8.13) with i*uc -> uc (opposite sign convention from Hans' website)
    ], η; name)
end

function perturbations_gravity(; name)
    @variables δρ(η) Δm(η) ρm(η) Π(η)
    return ODESystem([
        Dη(Φ) ~ (3/2*a^2*δρ - k^2*Φ + 3*(Dη(a)/a)^2*Ψ) / (3*Dη(a)/a) # Dodelson (6.41) # TODO: write in more natural form?
        k^2 * (Ψ + Φ) ~ Π # anisotropic stress
        Δm ~ k^2*Φ / (3/2*a^2*ρm) # gauge-invariant overdensity (from Poisson equation) # TODO: move outside or change?
    ], η; name)
end

lmax = 6
@named rad = perturbations_photon_hierarchy(lmax, true)
@named pol = perturbations_polarization_hierarchy(lmax)
@named cdm = perturbations_matter(false)
@named bar = perturbations_matter(true)
@named grav = perturbations_gravity()
@parameters fb dτspline # TODO: get rid of
@variables ρc(η) ρb(η) δρr(η) δρc(η) δρb(η) R(η) dτ(η) # TODO: get rid of
dτfunc(η, spl) = -exp(spl(log(η))) # TODO: type-stable? @code_warntype dτfunc(1e0) gives warnings, but code seems fast?
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
    [rad.interactions[l] ~ dτ * (rad.Θ[l] - rad.Θ[2]/10*δkron(l,2)) for l in 2:lastindex(rad.interactions)]...
    dτ ~ dτfunc(η, dτspline) # TODO: spline over η

    # polarization
    pol.dτ ~ dτ
    pol.Π ~ rad.Θ[2] + pol.Θ[2] + pol.Θ0

    # gravity shear stress
    grav.Π ~ -12*a^2 * bg.rad.ρ*rad.Θ[2] # TODO: add neutrinos
], η, [Ψ, Φ], [k, fb, dτspline])
@named pt = compose(pt_bg_conn, rad, pol, bar, cdm, grav, bg)
pt_sim = structural_simplify(pt)
pt_prob = ODEProblem(pt_sim, unknowns(pt_sim) .=> NaN, (ηi(bg_sol), 4.0), parameters(pt_sim) .=> NaN; jac=true) 

function solve_perturbations(kvals::AbstractArray, Ωr0, Ωm0, Ωb0, H0, Yp)
    bg_sol = solve_background(Ωr0, Ωm0)
    ηini, ηtoday = ηi(bg_sol), η0(bg_sol)
    Ωrini, Ωmini, ΩΛini, Eini = bg_sol(ηini; idxs=[bg.rad.ρ, bg.mat.ρ, bg.de.ρ, bg.E]) # TODO: avoid duplicate logic

    fb = Ωb0 / Ωm0; @assert fb <= 1 # TODO: avoid duplication thermo logic
    th_sol = solve_thermodynamics(Ωr0, Ωm0, Ωb0, H0, Yp) # update spline for dτ (e.g. to propagate derivative information through recombination, if called with dual numbers) TODO: use th_sol(a; idxs=th.dτ) directly in a type-stable way?
    dτspline = CubicSpline(log.(-th_sol.(th_sol[η], Val{1}, idxs=th.τ)), log.(th_sol[η]); extrapolate=true) # update spline for dτ (e.g. to propagate derivative information through recombination, if called with dual numbers) # TODO: verify that extrapolate=true is ok (it is needed for autodiff CMB computation) # TODO: use th_sol(a; idxs=th.dτ) directly in a type-stable way?
    dτini = dτfunc(ηini, dτspline)

    function prob_func(_, i, _)
        kval = kvals[i]
        println("$i/$(length(kvals)) k = $(kval*k0) Mpc/h")
        Φini = 2/3 # TODO: why 2/3? # arbitrary normalization (from primordial curvature power spectrum?)
        lmax = lastindex(pt_sim.rad.Θ)
        Θrini0 = Φini/2 # Dodelson (7.89)
        Θrini = Vector{Any}(undef, lmax)
        Θrini[1] = -kval*Φini/(6*aini*Eini) # Dodelson (7.95)
        Θrini[2] = -8/15#=-20/45=#*kval/dτini * Θrini[1]# TODO: change with/without polarization; another argument for merging photons+polarization
        for l in 3:lmax
            Θrini[l] = -l/(2*l+1) * kval/dτini * Θrini[l-1]
        end
        ΘPini0 = 5/4 * Θrini[2]
        ΘPini = Vector{Any}(undef, lmax) # TODO: allow lrmax ≠ lPmax
        ΘPini[1] = -kval/(4*dτini) * Θrini[2]
        ΘPini[2] = 1/4 * Θrini[2]
        for l in 3:lmax
            ΘPini[l] = -l/(2*l+1) * kval/dτini * ΘPini[l-1]
        end
        δcini = δbini = 3*Θrini0 # Dodelson (7.94)
        ucini = ubini = 3*Θrini[1] # Dodelson (7.95)
        return remake(pt_prob; tspan = (ηini, ηtoday), u0 = Dict(Φ => Φini, pt_sim.rad.Θ0 => Θrini0, [pt_sim.rad.Θ[l] => Θrini[l] for l in 1:lmax]..., pt_sim.pol.Θ0 => ΘPini0, [pt_sim.pol.Θ[l] => ΘPini[l] for l in 1:lmax]..., pt_sim.bar.δ => δbini, pt_sim.bar.u => ubini, pt_sim.cdm.δ => δcini, pt_sim.cdm.u => ucini, bg.rad.ρ => Ωrini, bg.mat.ρ => Ωmini, bg.de.ρ => ΩΛini, bg.a => aini), p = [pt_sim.fb => fb, pt_sim.k => kval, pt_sim.dτspline => dτspline])
    end

    probs = EnsembleProblem(prob = nothing, prob_func = prob_func)
    return solve(probs, KenCarp4(), EnsembleThreads(), reltol=1e-8, trajectories = length(kvals)) # KenCarp4 and Kvaerno5 seem to work well # TODO: test GPU parallellization
end
solve_perturbations(kvals::AbstractArray, θ::Parameters) = solve_perturbations(kvals, θ.ρr0, θ.ρm0, θ.ρb0, θ.H0, θ.Yp)

#=
pt_sols = solve_perturbations(ks, par.Ωr0, par.Ωm0, par.Ωb0, par.H0, par.Yp)
for (i, pt_sol) in enumerate(pt_sols)
    plot!(p[4], log10.(pt_sol[a]), pt_sol[Φ]; xlabel="lg(a)", ylabel="Φ/Φᵢ")
    plot!(p[5], log10.(pt_sol[a]), log10.(abs.(pt_sol[pt.cdm.δ])); color=i, xlabel="lg(a)", ylabel="lg(|δb|), lg(δc)")
    plot!(p[5], log10.(pt_sol[a]), log10.(abs.(pt_sol[pt.bar.δ])); color=i, xlabel="lg(a)", ylabel="lg(|δb|), lg(δc)")
end
display(p)
=#

# power spectra
θ0 = [par.Ωr0, par.Ωm0, par.Ωb0, par.H0, par.As, par.Yp]
P0(k, As) = @. As / k ^ 3
function P(k, Ωr0, Ωm0, Ωb0, H0, As, Yp)
    bg_sol = solve_background(Ωr0, Ωm0)
    return P0(k, As) .* solve_perturbations(k, Ωr0, Ωm0, Ωb0, H0, Yp)(η0(bg_sol); idxs=pt.grav.Δm) .^ 2
end
P(k, θ) = P(k, θ...) # unpack parameters θ = [ρr0, ρm0, ρb0, H0, As]

function plot_dlgP_dθs(dlgP_dθs, name, color)
    plot!(p[7], log10.(ks*k0), dlgP_dθs[:,1]; xlabel="lg(k/(h/Mpc))", ylabel="d lg(P) / d lg(Ωr0)", color=color, label=name)
    plot!(p[8], log10.(ks*k0), dlgP_dθs[:,2]; xlabel="lg(k/(h/Mpc))", ylabel="d lg(P) / d lg(Ωm0)", color=color, label=name)
    plot!(p[9], log10.(ks*k0), dlgP_dθs[:,4]; xlabel="lg(k/(h/Mpc))", ylabel="d lg(P) / d lg(H0)", color=color, label=name)
    display(p)
end

#=
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
=#

Ωr0, Ωm0, Ωb0, H0, Yp, As = par.Ωr0, par.Ωm0, par.Ωb0, par.H0, par.Yp, par.As
ηs = exp.(range(log(ηi(bg_sol)), log(η0(bg_sol)), length=800)) # logarithmic spread to capture early-time oscillations
ls = unique([2:2:10; 10:16:2010])
ks = range(1, 1.5*maximum(ls), step=2*π/4) ./ η0(bg_sol)
Sspline_ks = range(1, 1.5*maximum(ls), step=50) ./ η0(bg_sol) # Δk = 10/η0
# TODO: only need as from a = 1e-4 till today
# TODO: write to be more memory-efficient! only need to hold ∂Θ/∂ηs(η,k) for each l!
function S(ηs::AbstractArray, ks::AbstractArray, Ωr0, Ωm0, Ωb0, H0, Yp; Sspline_ks=nothing)    
    Ss = zeros(eltype([Ωr0, Ωm0, Ωb0, H0, Yp]), (length(ηs), length(ks)))

    if !isnothing(Sspline_ks)
        Ssplinedata = S(ηs, Sspline_ks, Ωr0, Ωm0, Ωb0, H0, Yp)
        for i_η in eachindex(ηs)
            Sspline = CubicSpline(Ssplinedata[i_η,:], Sspline_ks)
            Ss[i_η,:] .= Sspline(ks)
        end
        return Ss
    end

    th_sol = solve_thermodynamics(Ωr0, Ωm0, Ωb0, H0, Yp)
    η0 = th_sol.t[end]
    τ(η) = th_sol(η, idxs=th.τ).u .- th_sol(η0, idxs=th.τ)
    τ′(η) = th_sol(η, Val{1}, idxs=th.τ).u
    τ′′(η) = th_sol(η, Val{2}, idxs=th.τ).u
    g(η) = .-τ′(η) .* exp.(.-τ(η))
    g′(η) = (τ′(η) .^ 2 - τ′′(η)) .* exp.(.-τ(η))
    
    pt_sols = solve_perturbations(ks, Ωr0, Ωm0, Ωb0, H0, Yp)
    for (i_k, (k, pt_sol)) in enumerate(zip(ks, pt_sols))
        # TODO: must be faster!! use saveat for ηs in ODESolution?
        # TODO: add source functions as observed perturbation functions? but difficult with cumulative τ(η)? must anyway wait for this to be fixed: https://github.com/SciML/ModelingToolkit.jl/issues/2697
        Θ0(η) = pt_sol(η, idxs=pt.rad.Θ0)
        Ψ(η) = pt_sol(η, idxs=pt.Ψ)
        Φ(η) = pt_sol(η, idxs=pt.Φ)
        Π(η) = pt_sol(η, idxs=pt.pol.Π)
        ub(η) = pt_sol(η, idxs=pt.bar.u)
        function Ψ′(η)
            # TODO: use pt_sol(...) when this is fixed: https://github.com/SciML/ModelingToolkit.jl/issues/2697 and https://github.com/SciML/ModelingToolkit.jl/pull/2574
            #return pt_sol(η, Val{1}, idxs=pt.Ψ).u
            # workaround: spline and take derivative
            Ψspl = CubicSpline(Ψ(η).u, η; extrapolate=true)
            return DataInterpolations.derivative.(Ref(Ψspl), η)
        end
        Φ′(η) = pt_sol(η, Val{1}, idxs=pt.Φ).u
        ub′(η) = pt_sol(η, Val{1}, idxs=pt.bar.u).u

        Ss[:,i_k] .+= g(ηs) .* (Θ0(ηs) + Ψ(ηs) + Π(ηs)/4) # SW
        Ss[:,i_k] .+= (g′(ηs).*ub(ηs) + g(ηs).*ub′(ηs)) / k # Doppler
        Ss[:,i_k] .+= exp.(.-τ(ηs)) .* (Ψ′(ηs) - Φ′(ηs)) # ISW
        # TODO: add polarization
    end

    return Ss
end

#Ss = S(ηs, ks, Ωr0, Ωm0, Ωb0, H0, Yp)
#plot(ηs, asinh.(Ss[:,[1,9]]))

function ∂Θ_∂η(ls::AbstractArray, ks::AbstractArray, ηs::AbstractArray, Ωr0, Ωm0, Ωb0, H0, Yp; kwargs...)
    # TODO: can move into ΘT() to avoid allocating for third l-dimension!
    ∂Θ_∂ηs = zeros(eltype([Ωr0, Ωm0, Ωb0, H0, Yp]), (length(ηs), length(ks), length(ls)))
    Ss = S(ηs, ks, Ωr0, Ωm0, Ωb0, H0, Yp; kwargs...)
    η0 = ηs[end] # TODO: assume!!
    for (i_l, l) in enumerate(ls)
        for (i_k, k) in enumerate(ks)
            for (i_η, η) in enumerate(ηs)
                # TODO: faster to calculate with a range for l? but only for besselj, not for sphericalbesselj. https://github.com/JuliaMath/Bessels.jl?tab=readme-ov-file#support-for-sequence-of-orders
                ∂Θ_∂ηs[i_η,i_k,i_l] = Ss[i_η,i_k] * sphericalbesselj(l, k * (η0-η))
            end
        end
    end
    return ∂Θ_∂ηs
end

#∂Θ_∂ηs = ∂Θ_∂η(ls, ks, ηs, Ωr0, Ωm0, Ωb0, H0, Yp; Sspline_ks)

function ΘT(ls::AbstractArray, ks::AbstractArray, ηs::AbstractArray, Ωr0, Ωm0, Ωb0, H0, Yp; kwargs...)
    # TODO: just integrate the spline! https://discourse.julialang.org/t/how-to-speed-up-the-numerical-integration-with-interpolation/96223/5
    ∂Θ_∂ηs = ∂Θ_∂η(ls, ks, ηs, Ωr0, Ωm0, Ωb0, H0, Yp; kwargs...)
    return [trapz(ηs, ∂Θ_∂ηs[:,i_k,i_l]) for i_l in eachindex(ls), i_k in eachindex(ks)] # TODO: return all Θls in shape (size(ls), size(ks))
end

#Θs = ΘT(ls, ks, ηs, Ωr0, Ωm0, Ωb0, H0, Yp)

function dCl_dk(ls::AbstractArray, ks::AbstractArray, ηs::AbstractArray, Ωr0, Ωm0, Ωb0, H0, As, Yp; kwargs...)
    Θls = ΘT(ls, ks, ηs, Ωr0, Ωm0, Ωb0, H0, Yp; kwargs...)
    # TODO: just integrate the spline! https://discourse.julialang.org/t/how-to-speed-up-the-numerical-integration-with-interpolation/96223/5
    return stack([@. 2/π * ks^2 * P0(ks, As) * Θls[i_l,:]^2 for i_l in eachindex(ls)])
end

# TODO: integrate over log(a) instead of a!
function Cl(ls::AbstractArray, ks::AbstractArray, ηs::AbstractArray, Ωr0, Ωm0, Ωb0, H0, As, Yp; kwargs...)
    dCl_dks = dCl_dk(ls, ks, ηs, Ωr0, Ωm0, Ωb0, H0, As, Yp; kwargs...)
    # TODO: just integrate the spline! https://discourse.julialang.org/t/how-to-speed-up-the-numerical-integration-with-interpolation/96223/5
    Cls = trapz(ks, dCl_dks, Val(1)) # integrate over k
    Cls .+= @. (ks[1] - 0.0) * (0.0 + dCl_dks[1,:]) # add extra trapz point at (k, dCl_dk) = (0.0, 0.0)
    return Cls
end

function Dl(ls::AbstractArray, ks::AbstractArray, ηs::AbstractArray, Ωr0, Ωm0, Ωb0, H0, As, Yp; kwargs...)
    return Cl(ls, ks, ηs, Ωr0, Ωm0, Ωb0, H0, As, Yp; kwargs...) .* ls .* (ls .+ 1) / (2*π)
end
Dl(ls::AbstractArray, ks::AbstractArray, ηs::AbstractArray, θ; kwargs...) = Dl(ls, ks, ηs, θ...; kwargs...) # unpack parameters θ = [ρr0, ρm0, ρb0, H0, As]

#Dls = Dl(ls, ks, ηs, Ωr0, Ωm0, Ωb0, H0, As, Yp; Sspline_ks)
#plot(ls, Dls; xlabel="l", ylabel="Dl = l (l+1) Cl / 2π")

# differentiated CMB power spectrum
lgDlres = DiffResults.JacobianResult(Float64.(ls), θ0)
lgDl(lgθ) = log10.(Dl(ls, ks, ηs, 10 .^ lgθ; Sspline_ks))
ForwardDiff.jacobian!(lgDlres, lgDl, log10.(θ0))
lgDls, dlgDl_dθs_ad = DiffResults.value(lgDlres), DiffResults.jacobian(lgDlres)

p = plot(layout=(2,1), size=(800, 1000), left_margin=bottom_margin=30*Plots.px); display(p)
plot!(p[1], ls, 10 .^ lgDls / 1e-12; xlabel = "l", ylabel = "Dₗ=l (l+1) Cₗ / 2π / 10⁻¹²", title = "CMB power spectrum"); display(p)
plot!(p[2], ls, dlgDl_dθs_ad; xlabel = "l", ylabel = "d lg(Dₗ) / d lg(θᵢ)", labels = "θᵢ=" .* ["Ωr0" "Ωm0" "Ωb0" "H0" "As" "Yp"]); display(p)