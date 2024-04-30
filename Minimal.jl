using ModelingToolkit
using DifferentialEquations
using DataInterpolations
using ForwardDiff, DiffResults, FiniteDiff
using OffsetArrays
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
    Ωr0 = 5e-5 # Ωr0
    Ωm0 = 0.3 # Ωm0
    Ωb0 = 0.02 # Ωb0
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
@variables η a(η) E(η) # η is really η in units of 1/H0
a = GlobalScope(a)
Dη = Differential(η)

function background_gravity_GR(; name)
    @variables ρ(η)
    return ODESystem([
        Dη(a) ~ √(ρ * a^4) # TODO: 8π/3 factor?
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
solve_background(θ::Parameters) = solve_background(θ.ρr0, θ.ρm0)
ηi(bg_sol::ODESolution) = bg_sol.prob.tspan[1]
η0(bg_sol::ODESolution) = find_zero(η -> bg_sol(η, idxs=a) - 1.0, bg_sol.prob.tspan)

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
    @variables Xe(η) ne(η) nH(η) T(η) H(η) λ(η) α2(η) β(η) C(η) Λα(η) Λ2γ(η) β2(η)
    return ODESystem([
        ne ~ Xe * nH
        λ ~ h / √(2π*me*kB*T) # e⁻ de-Broglie wavelength
        α2 ~ 9.78 * (α*ħ/me)^2/c * √(EHion/(kB*T)) * log(EHion/(kB*T)) # Dodelson (4.38) (e⁻ + p → H + γ) # TODO: add Recfast fudge factor?
        β  ~ α2 / λ^3 * exp(-EHion/(  kB*T)) # Dodelson (4.37)-(4.38) (γ + H → e⁻ + p)
        β2 ~ α2 / λ^3 * exp(-EHion/(4*kB*T)) # 1/s (compute this instead of β2 = β * exp(3*EHion/(4*kB*T)) to avoid exp overflow)
        Λα ~ H * (3*EHion/(ħ*c))^3 / ((8*π)^2 * nH) # 1/s
        Λ2γ ~ 8.227 # 1/s
        C ~ (Λ2γ + Λα) / (Λ2γ + Λα + β2) # Peebles' correction factor (Dodelson exercise 4.7)
        Dη(Xe) * H0 ~ C * ((1-Xe)*β - Xe^2*nH*α2) * a # remains ≈ 0 during Saha recombinations, so no need to manually turn off (multiply by H0 on left because cide η is physical η/(1/H0))
    ], η, [Xe, H, nH, ne, T], [H0]; name)
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
    
    # switch *smoothly* from Saha to Peebles when XeS ≤ 1 (see e.g. https://discourse.julialang.org/t/handling-instability-when-solving-ode-problems/9019/5) # TODO: make into a connection
    Xe ~ Hifelse(1-saha.Xe, saha.Xe, peebles.Xe; k=1e3) + reion1.Xe + reion2.Xe
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
    return solve(prob, RadauIIA5(), reltol=1e-7) # CLASS uses "NDF15" (https://lesgourg.github.io/class-tour/London2014/Numerical_Methods_in_CLASS_London.pdf) TODO: after switching ivar from a to b=ln(a), the integrator needs more steps. fix this?
end
solve_thermodynamics(θ::Parameters) = solve_thermodynamics(θ.Ωr0, θ.Ωm0, θ.Ωb0, θ.H0, θ.Yp)

th_sol = solve_thermodynamics(par.Ωr0, par.Ωm0, par.Ωb0, par.H0, par.Yp)
plot!(p[3], log10.(th_sol[a]), stack(th_sol[[saha.Xe, peebles.Xe, reion1.Xe, reion2.Xe, th_sim.Xe]])'; xlabel="lg(a)", ylabel="Xe", ylims=(0, 1.5), label=["XeS" "XeP" "XeRE1" "XeRE2" "Xe"], legend=:bottomleft); display(p)
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
    @variables Θ(η)[0:lmax] δ(η) interactions(η)[1:lmax-1]
    eqs = [
        Dη(Θ[0]) + k*Θ[1] ~ -Dη(Φ)
        Dη(Θ[1]) - k/3*(Θ[0]-2*Θ[2]) ~ k/3*Ψ + interactions[1]
        [Dη(Θ[l]) ~ k/(2*l+1) * (l*Θ[l-1] - (l+1)*Θ[l+1]) + interactions[l] for l in 2:lmax-1]...
        Θ[lmax] ~ 0 # TODO: integrate η(a) and use better cutoff
        δ ~ 4*Θ[0]
    ]
    if !interact
        push!(eqs, collect(interactions .~ 0)...)
    end
    return ODESystem(eqs, η; name)
end

function perturbations_polarization_hierarchy(lmax=6; name)
    @variables Θ(η)[0:lmax] dτ(η) Π(η)
    eqs = [
        Dη(Θ[0]) + k*Θ[1] ~ dτ * (Θ[0] - Π/2)
        [Dη(Θ[l]) - k/(2*l+1) * (l*Θ[l-1] - (l+1)*Θ[l+1]) ~ dτ * (Θ[l] - Π/10*δkron(l,2)) for l in 1:lmax-1]...
        Θ[lmax] ~ 0 # TODO: integrate η(a) and use better cutoff
    ]
    return ODESystem(eqs, η; name)
end

function perturbations_matter(interact=false; name)
    @variables δ(η) u(η)
    interaction = interact ? only(@variables interaction(η)) : 0
    return ODESystem([
        Dη(δ) + k*u ~ -3*Dη(Φ) # Dodelson (5.69) or (8.12) with i*uc -> uc
        Dη(u) + u*Dη(a)/a ~ k*Ψ + interaction # Dodelson (5.70) or (8.13) with i*uc -> uc
    ], η; name)
end

function perturbations_gravity(; name)
    @variables δρ(η) Δm(η) ρm(η) Π(η)
    return ODESystem([
        Dη(Φ) ~ (3/2*a^2*δρ - k^2*Φ - 3*(Dη(a)/a)^2*Φ) / (3*Dη(a)/a) # Dodelson (8.14) # TODO: write in more natural form?
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
    [rad.interactions[l] ~ dτ * (rad.Θ[l] - pol.Π/10*δkron(l,2)) for l in 2:lastindex(rad.interactions)]...
    dτ ~ dτfunc(η, dτspline) # TODO: spline over η

    # polarization
    pol.dτ ~ dτ
    pol.Π ~ rad.Θ[2] + pol.Θ[2] + pol.Θ[0]

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
    dτspline = CubicSpline(log.(-th_sol.(th_sol[η], Val{1}, idxs=th.τ)), log.(th_sol[η])) # update spline for dτ (e.g. to propagate derivative information through recombination, if called with dual numbers) TODO: use th_sol(a; idxs=th.dτ) directly in a type-stable way?
    dτini = dτfunc(ηini, dτspline)

    function prob_func(_, i, _)
        kval = kvals[i]
        println("$i/$(length(kvals)) k = $(kval*k0) Mpc/h")
        Φini = 1.0 # arbitrary normalization (from primordial curvature power spectrum?)
        lmax = lastindex(pt_sim.rad.Θ)
        Θrini = OffsetVector(Vector{Any}(undef, lmax), -1) # index from l=0 to l=lmax-1
        Θrini[0] = Φini/2 # Dodelson (7.89)
        Θrini[1] = -kval*Φini/(6*aini*Eini) # Dodelson (7.95)
        Θrini[2] = -8/15*kval/(aini^2*Eini*dτini) # TODO: change with/without polarization
        for l in 3:lmax-1
            Θrini[l] = -l/(2*l+1) * kval/(aini^2*Eini*dτini) * Θrini[l-1]
        end
        ΘPini = OffsetVector(Vector{Any}(undef, lmax), -1) # index from l=0 to l=lmax-1 # TODO: allow lrmax ≠ lPmax
        ΘPini[0] = 5/4 * Θrini[2]
        ΘPini[1] = -kval/(4*aini^2*Eini*dτini) * Θrini[2]
        ΘPini[2] = 1/4 * Θrini[2]
        for l in 3:lmax-1
            ΘPini[l] = -l/(2*l+1) * kval/(aini^2*Eini*dτini) * ΘPini[l-1]
        end
        δcini = δbini = 3*Θrini[0] # Dodelson (7.94)
        ucini = ubini = 3*Θrini[1] # Dodelson (7.95)
        return remake(pt_prob; tspan = (ηini, ηtoday), u0 = Dict(Φ => Φini, [pt_sim.rad.Θ[l] => Θrini[l] for l in 0:lmax-1]..., [pt_sim.pol.Θ[l] => ΘPini[l] for l in 0:lmax-1]..., pt_sim.bar.δ => δbini, pt_sim.bar.u => ubini, pt_sim.cdm.δ => δcini, pt_sim.cdm.u => ucini, bg.rad.ρ => Ωrini, bg.mat.ρ => Ωmini, bg.de.ρ => ΩΛini, bg.a => aini), p = [pt_sim.fb => fb, pt_sim.k => kval, pt_sim.dτspline => dτspline])
    end

    probs = EnsembleProblem(prob = nothing, prob_func = prob_func)
    return solve(probs, KenCarp4(), EnsembleThreads(), reltol=1e-8, trajectories = length(kvals)) # KenCarp4 and Kvaerno5 seem to work well # TODO: test GPU parallellization
end
solve_perturbations(kvals::AbstractArray, θ::Parameters) = solve_perturbations(kvals, θ.ρr0, θ.ρm0, θ.ρb0, θ.H0, θ.Yp)

pt_sols = solve_perturbations(ks, par.Ωr0, par.Ωm0, par.Ωb0, par.H0, par.Yp)
for (i, pt_sol) in enumerate(pt_sols)
    plot!(p[4], log10.(pt_sol[a]), pt_sol[Φ]; xlabel="lg(a)", ylabel="Φ/Φᵢ")
    plot!(p[5], log10.(pt_sol[a]), log10.(abs.(pt_sol[pt.cdm.δ])); color=i, xlabel="lg(a)", ylabel="lg(|δb|), lg(δc)")
    plot!(p[5], log10.(pt_sol[a]), log10.(abs.(pt_sol[pt.bar.δ])); color=i, xlabel="lg(a)", ylabel="lg(|δb|), lg(δc)")
end
display(p)

# power spectra
θ0 = [par.Ωr0, par.Ωm0, par.Ωb0, par.H0, par.As, par.Yp]
P0(k, As) = @. As / k ^ 3
P(k, Ωr0, Ωm0, Ωb0, H0, As, Yp) = P0(k, As) .* solve_perturbations(k, Ωr0, Ωm0, Ωb0, H0, Yp)(atoday; idxs=pt.grav.Δm) .^ 2
P(k, θ) = P(k, θ...) # unpack parameters θ = [ρr0, ρm0, ρb0, H0, As]

#=
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
=#

# TODO: CMB power spectrum
#=
ρr0, ρm0, ρb0, H0, Yp, As = par.ρr0, par.ρm0, par.ρb0, par.H0, par.Yp, par.As
# TODO: only need as from a = 1e-4 till today
# TODO: spline_first logic for each k!
function S(xs::AbstractArray, ks::AbstractArray, ρr0, ρm0, ρb0, H0, Yp)
    Ss = zeros(length(xs), length(ks))
    D = ForwardDiff.derivative

    th_sol = solve_thermodynamics(ρr0, ρm0, ρb0, H0, Yp)
    a(x) = exp.(x) # TODO: move dot to callers
    E(x) = th_sol(a(x), idxs=bg.E)
    τ(x) = th_sol(a(x), idxs=th.τ) .- th_sol(atoday, idxs=th.τ)
    g(x) = D(x -> exp(-τ(x)), x)

    pt_sols = solve_perturbations(ks, ρr0, ρm0, ρb0, H0, Yp)
    for (i_k, (k, pt_sol)) in enumerate(zip(ks, pt_sols))
        Θ0(x) = pt_sol(a(x), idxs=pt.rad.Θ[0])
        Ψ(x) = pt_sol(a(x), idxs=pt.Ψ)
        Φ(x) = pt_sol(a(x), idxs=pt.Φ)
        Π(x) = pt_sol(a(x), idxs=pt.grav.Π)
        ub(x) = pt_sol(a(x), idxs=pt.bar.u)

        S_SW = g.(xs) .* (Θ0(xs) + Ψ(xs) + Π(xs)/4) # TODO: avoid g elementwise .?
        S_ISW = exp.(.-τ(xs)) .* D.(x -> Ψ(x) - Φ(x), xs) # TODO: oscillates unless low tolerance in ODE solver? but removed after adding Doppler?
        S_Doppler = -1/k * D.(x -> a(x) * E(x) * g(x) * ub(x), xs)
        Ss[:,i_k] = S_SW + S_ISW + S_Doppler
    end

    return Ss
end

function ΘT(ls::AbstractArray, ks::AbstractArray, xs::AbstractArray, ρr0, ρm0, ρb0, H0, Yp)
    Ss = S(xs, ks, ρr0, ρm0, ρb0, H0, Yp)

    bg_sol = solve_background(ρr0, ρm0)
    a(x) = exp(x)
    Δηs = bg_sol(atoday; idxs=bg.η) .- bg_sol(a.(xs); idxs=bg.η)
    ys = ks' .* Δηs # argument to Bessel function

    # TODO: transform integral to log(a)
    # TODO: just integrate the spline! https://discourse.julialang.org/t/how-to-speed-up-the-numerical-integration-with-interpolation/96223/5
    ∂Θ_∂xs = Ss .* stack(sphericalbesselj.(l, ys) for l in ls)
    println(size(∂Θ_∂xs)) # (as, ks, ls)
    return [trapz(xs, ∂Θ_∂xs[:,i_k,i_l]) for i_l in eachindex(ls), i_k in eachindex(ks)] # TODO: return all Θls in shape (size(ls), size(ks))
end

# TODO: integrate over log(a) instead of a!
function Cl(ls::AbstractArray, ks::AbstractArray, as::AbstractArray, ρr0, ρm0, ρb0, H0, As, Yp)
    Θls = ΘT(ls, ks, log.(as), ρr0, ρm0, ρb0, H0, Yp)
    # TODO: just integrate the spline! https://discourse.julialang.org/t/how-to-speed-up-the-numerical-integration-with-interpolation/96223/5
    return [2/π .* trapz(ks, @. ks^2 * P0(ks, As) * Θls[i_l,:]^2) for i_l in eachindex(ls)]
end

function Dl(ls::AbstractArray, ks::AbstractArray, as::AbstractArray, ρr0, ρm0, ρb0, H0, As, Yp)
    return Cl(ls, ks, as, ρr0, ρm0, ρb0, H0, As, Yp) .* ls .* (ls .+ 1) / (2*π)
end

#=
ks = 10 .^ range(-4, +2, length=10) / k0 # in code units of k0 = H0/c
as = 10 .^ range(-4, 0, length=1000)
Ss = S(as, ks, ρr0, ρm0, ρb0, H0, Yp)
plot()
#plot!(log10.(as), Ss[:,1])
plot!(log10.(as), Ss[:,9])
=#

lmax = 1000
η0 = -bg_sol(aini, idxs=bg.η)
kη0s = range(1, 2*lmax, step=2*π/8) # TODO: stop should be *higher* than lmax
ks = kη0s / η0
as = 10 .^ range(-4, 0, length=600)
ls = range(1, lmax, step=10)

Dls = Dl(ls, ks, as, ρr0, ρm0, ρb0, H0, As, Yp)
plot(log10.(ls), Dls)
=#