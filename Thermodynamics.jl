struct ThermodynamicsSystem
    sys::ODESystem
    ssys::ODESystem
    prob::ODEProblem
    bg::BackgroundSystem
end

# TODO: separate components more meaningfully
# TODO: Hydrogen, Helium, baryons, photons, ...

# background thermodynamics / recombination
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

function recombination_hydrogen_peebles(g; name)
    @parameters H0
    @variables Xe(η) nH(η) T(η) H(η) λ(η) α2(η) β(η) C(η) Λ2γ_Λα(η) β2_Λα(η) actiweight(η)
    return ODESystem([
        λ ~ h / √(2π*me*kB*T) # e⁻ de-Broglie wavelength
        α2 ~ 9.78 * (α*ħ/me)^2/c * √(EHion/(kB*T)) * log(EHion/(kB*T)) # Dodelson (4.38) (e⁻ + p → H + γ) # TODO: add Recfast fudge factor?
        β  ~ α2 / λ^3 * exp(-EHion/(kB*T)) # Dodelson (4.37)-(4.38) (γ + H → e⁻ + p)
        β2_Λα ~ α2 / λ^3 * exp(-EHion/(4*kB*T)) * ((8*π)^2 * (1-Xe) * nH) / (H * (3*EHion/(ħ*c))^3) # β2/Λα (compute this instead of β2 = β * exp(3*EHion/(4*kB*T)) to avoid exp overflow)
        Λ2γ_Λα ~ 8.227 * ((8*π)^2 * (1-Xe) * nH) / (H * (3*EHion/(ħ*c))^3) # Λ2γ/Λα
        C ~ Hifelse(actiweight, 1, (1 + Λ2γ_Λα) / (1 + Λ2γ_Λα + β2_Λα); k=1e3) # Peebles' correction factor (Dodelson exercise 4.7), manually activated to avoid numerical issues at early times # TODO: activate using internal quantities only! # TODO: why doesnt it work to activate from 0? activating from 1 is really unnatural
        Dη(Xe) ~ C * ((1-Xe)*β - Xe^2*nH*α2) * g.a / H0 # remains ≈ 0 during Saha recombinations, so no need to manually turn off (multiply by H0 on left because cide η is physical η/(1/H0))
    ], η, [Xe, H, nH, T, actiweight], [H0]; name)
end

function thermodynamics_temperature(g; name)
    @variables Tγ(η) Tb(η) ργ(η) ρb(η) τ(η) fγb(η)
    return ODESystem([
        Dη(Tγ) ~ -1*Tγ * g.ℰ # Tγ = Tγ0 / a # TODO: introduce ℋ = Dη(a) / a?
        Dη(Tb) ~ -2*Tb * g.ℰ - 8/3*(mp/me)*fγb*g.a*Dη(τ)*(Tγ-Tb) # TODO: multiply last term by a or not?
    ], η, [Tγ, Tb, τ, fγb], []; name)
end

function reionization_smooth_step(g, z0, Δz0, Xe0; name)
    y(z) = (1+z)^(3/2)
    Δy(z, Δz0) = 3/2 * (1+z)^(1/2) * Δz0
    @variables Xe(η) z(η)
    return ODESystem([
        z ~ 1/g.a - 1
        Xe ~ Hifelse(y(z0)-y(z), 0, Xe0; k=1/Δy(z0, Δz0)) # smooth step from 0 to Xe0
    ], η; name)
end

function ThermodynamicsSystem(bg::BackgroundSystem, Herec::ODESystem, Hrec::ODESystem, temp::ODESystem, reions::AbstractArray{ODESystem}; name)
    @parameters fb H0
    @variables Xe(η) ne(η) τ(η) H(η) ρb(η) nb(η) nH(η)
    connections = ODESystem([
        H ~ bg.sys.g.E * H0 # 1/s # TODO: avoid duplicate name with background H
        ρb ~ fb * bg.sys.mat.ρ * H0^2/G # kg/m³
        nb ~ ρb / mp # 1/m³
        nH ~ (1-Herec.Yp) * nb # TODO: correct?

        temp.fγb ~ bg.sys.rad.ρ / (fb*bg.sys.mat.ρ) # ργ/ρb
        temp.τ ~ τ
        Herec.T ~ temp.Tb
        Herec.nH ~ nH
        Hrec.T ~ temp.Tb
        Hrec.nH ~ nH
        Hrec.H ~ H
        Hrec.actiweight ~ 1 - Herec.Xe
        
        # switch *smoothly* from Saha to Peebles when XeS ≤ 1 (see e.g. https://discourse.julialang.org/t/handling-instability-when-solving-ode-problems/9019/5) # TODO: make into a connection
        Xe ~ Hifelse(Hrec.actiweight, Herec.Xe, Hrec.Xe; k=1e3) + sum(reion.Xe for reion in reions)
        ne ~ Xe * nH
        Dη(τ) * H0 ~ -ne * σT * c * bg.sys.g.a # common optical depth τ (multiply by H0 on left because code η is physical η/(1/H0)) # TODO: separate in Saha/Peebles?
    ], η; name)
    sys = compose(connections, [Herec; Hrec; temp; reions; bg.sys])
    ssys = structural_simplify(sys)
    prob = ODEProblem(ssys, unknowns(ssys) .=> NaN, (0.0, 4.0), parameters(ssys) .=> NaN; jac=true)
    return ThermodynamicsSystem(sys, ssys, prob, bg)
end

function solve(th::ThermodynamicsSystem, Ωr0, Ωm0, Ωb0, h, Yp; aini=1e-8, aend=1.0, solver=RadauIIA5(), reltol=1e-8, kwargs...)
    H0 = h * 100 * km/Mpc
    bg_sol = solve(th.bg, Ωr0, Ωm0)
    ηini, ηtoday = bg_sol[η][begin], bg_sol[η][end]
    ρrini, ρmini, ρΛini = bg_sol(ηini; idxs=[th.bg.sys.rad.ρ, th.bg.sys.mat.ρ, th.bg.sys.de.ρ]) # TODO: avoid duplicate logic
    fb = Ωb0 / Ωm0; @assert fb <= 1
    Tini = (ρrini * 15/π^2 * H0^2/G * ħ^3*c^5)^(1/4) / kB # common initial Tb = Tγ TODO: relate to ρr0 once that is a parameter
    XeSini = 1 + Yp / (4*(1-Yp)) * 2 # TODO: avoid?
    prob = remake(th.prob; tspan = (ηini, ηtoday), u0 = [th.ssys.Herec.Xe => XeSini, th.ssys.Hrec.Xe => 1.0, th.ssys.temp.Tγ => Tini, th.ssys.temp.Tb => Tini, th.ssys.τ => 0.0, th.bg.sys.rad.ρ => ρrini, th.bg.sys.mat.ρ => ρmini, th.bg.sys.de.ρ => ρΛini, th.bg.ssys.g.a => aini], p = [th.ssys.fb => fb, th.ssys.H0 => H0, th.ssys.Hrec.H0 => H0, th.ssys.Herec.Yp => Yp, th.ssys.reion1.Herec₊Yp => Yp, th.ssys.reion2.Herec₊Yp => Yp]) # TODO: avoid 2xH0, 2xYp, ...
    return solve(prob, solver; reltol, kwargs...) # CLASS uses "NDF15" (https://lesgourg.github.io/class-tour/London2014/Numerical_Methods_in_CLASS_London.pdf) TODO: after switching ivar from a to b=ln(a), the integrator needs more steps. fix this?
end