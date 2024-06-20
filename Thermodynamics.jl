import SpecialFunctions: zeta as ζ

struct ThermodynamicsSystem
    sys::ODESystem
    ssys::ODESystem
    bg::BackgroundSystem
end

function thermodynamics_hydrogen_recombination_saha(Eion = NoUnits(13.59844u"eV/J"); kwargs...)
    @parameters Y
    @variables X(η) Xe(η) T(η) λ(η) R(η) ne(η) n(η) nb(η) nγ(η) n(η)
    return ODESystem([
        n ~ Y * nb
        λ ~ h / √(2π*me*kB*T)
        R ~ exp(-Eion/(kB*T)) / (λ^3 * ne)
        X ~ R / (1 + R) # = nH⁺/nH
        Xe ~ Y * X
    ], η; defaults = [X => 1.0, Xe => Y], kwargs...)
end

function thermodynamics_helium_recombination_saha(; Eion1 = NoUnits(24.58738u"eV/J"), Eion2 = NoUnits(54.41776u"eV/J"), kwargs...)
    pars = @parameters Y
    vars = @variables Xe(η) X1(η) X2(η) T(η) λ(η) R1(η) R2(η) ne(η) n(η) nb(η) nγ(η)
    return ODESystem([
        λ ~ h / √(2π*me*kB*T)
        R1 ~ 2 * exp(-Eion1/(kB*T)) / (λ^3 * ne)
        R2 ~ 4 * exp(-Eion2/(kB*T)) / (λ^3 * ne)
        X1 ~ R1 / (1 + R1 + R1*R2)
        X2 ~ R2 * X1
        Xe ~ Y/4 * (X1 + 2*X2)
    ], η, vars, pars; defaults = [X2 => 1.0, X1 => 0.0, Xe => Y/2], kwargs...)
end

# background thermodynamics / recombination
smoothifelse(x, v1, v2; k=1) = 1/2 * ((v1+v2) + (v2-v1)*tanh(k*x)) # smooth transition/step function from v1 at x<0 to v2 at x>0
smoothmin(v1, v2; kwargs...) = smoothifelse(v1 - v2, v1, v2; kwargs...)
smoothmax(v1, v2; kwargs...) = smoothifelse(v1 - v2, v2, v1; kwargs...)

function thermodynamics_hydrogen_recombination_peebles(g; E1 = -NoUnits(13.59844u"eV/J"), kwargs...)
    pars = @parameters Y
    vars = @variables X(η) Xe(η) T(η) λe(η) ne(η) n(η) n1s(η) λ(η) λα(η) α2(η) β(η) C(η) Λα⁻¹(η) Λ2γ_Λα(η) β2_Λα(η) nb(η) nγ(η)
    Eion = -E1
    E(n) = E1 / n^2
    return ODESystem([
        n ~ Y * nb
        λe ~ h / √(2π*me*kB*T) # h / √(2π*me*kB*T) # e⁻ de-Broglie wavelength
        α2 ~ 9.78 * (α*ħ/me)^2/c * √(Eion/kB/T) * log(Eion/kB/T) # Dodelson (4.38) (e⁻ + p → H + γ) # TODO: add Recfast fudge factor?
        β ~ α2 * exp(E(1)/kB/T) / λe^3 # Dodelson (4.37)-(4.38) (γ + H → e⁻ + p)

        # TODO: fix Peebles' correction factor! it is unstable
        #n1s ~ smoothifelse(0.8 - X, 0, 1-X; k=1e2) * n # (1-X) * n # TODO: n1s sometimes goes negative!
        #λα ~ 8π*ħ*c / (3*Eion)
        #Λα⁻¹ ~ λα^3 * n1s / (8π * g.H) # TODO: this becomes negative, but should always be positive
        #β2_Λα ~ α2 / λe^3 * exp(E(2)/kB/T) * Λα⁻¹ # β2/Λα (compute this instead of β2 = β * exp(3*Eion/(4*kB*T)) to avoid exp overflow)
        #Λ2γ_Λα ~ 8.227 * Λα⁻¹ # Λ2γ/Λα
        C ~ 1 # (1 + Λ2γ_Λα) / (1 + Λ2γ_Λα + β2_Λα) # Peebles' correction factor (Dodelson exercise 4.7)

        Dη(X) * g.H0 ~ g.a * C * ((1-X) * β - α2*X^2*n) # TODO: do min(0, ...) to avoid increasing? # remains ≈ 0 during Saha recombinations, so no need to manually turn off (multiply by H0 on left because cide η is physical η/(1/H0))
        Xe ~ X * Y
    ], η, vars, pars; defaults = [X => 1.0, Xe => Y], kwargs...)
end
# testing:
# plot(log.(th_sol.t), th_sol[th.sys.H.Λα⁻¹])

# TODO: does this work correctly together with He?
function thermodynamics_hydrogen_recombination_baumann(g; Eion = NoUnits(13.59844u"eV/J"), ϵ=1e-20, kwargs...)
    @parameters Y λ
    @variables X(η) Xeq(η) T(η) ne(η) n(η) R⁻¹(η) x(η) nb(η) nγ(η)
    return ODESystem([
        # Baumann (3.3.108) (ϵ ≪ 1 maintains numerical stability)
        R⁻¹ ~ exp(-Eion/kB/T) / (2*ζ(3)/π^2 * nγ/nb * (2π*kB*T/(me*c^2))^(3/2))
        Xeq ~ -R⁻¹/2 + √(R⁻¹ + R⁻¹^2/4 + ϵ)

        # Baumann (3.3.123) (λ = λ(x=1) set as parameter)
        x ~ Eion / (kB*T)
        Dη(X) ~ -Dη(x) * λ/x^2 * (X^2 - Xeq^2)
    ], η, [X, Xeq, T, ne, n, R⁻¹, nb, nγ], [Y, λ]; defaults = [X => 1.0], kwargs...)
end

function reionization_tanh(g, z0, Δz0, Xe0; kwargs...)
    y(z) = (1+z)^(3/2)
    Δy(z, Δz0) = 3/2 * (1+z)^(1/2) * Δz0
    @variables Xe(η) z(η) T(η) ne(η) nb(η) nγ(η)
    @parameters Y
    @named base = ODESystem(Equation[], η, [T, ne, nb, nγ], [Y])
    return extend(ODESystem([
        z ~ 1/g.a - 1
        Xe ~ smoothifelse(y(z0)-y(z), 0, Xe0; k=1/Δy(z0, Δz0)) # smooth step from 0 to Xe0
    ], η; kwargs...), base)
end

function thermodynamics_ΛCDM(bg::BackgroundSystem; kwargs...)
    defaults = Dict() # TODO: merge with kwargs
    #@named H = thermodynamics_hydrogen_recombination_saha()
    @named H = thermodynamics_hydrogen_recombination_peebles(bg.sys.g)
    #@named H = thermodynamics_hydrogen_recombination_baumann(bg.sys.g)
    @named He = thermodynamics_helium_recombination_saha()
    @named Hre = reionization_tanh(bg.sys.g, 8.0, 0.5, H.Y); push!(defaults, Hre.H₊Y => H.Y) # TODO: avoid extra default?
    @named Here1 = reionization_tanh(bg.sys.g, 8.0, 0.5, He.Y/4); push!(defaults, Here1.He₊Y => He.Y)
    @named Here2 = reionization_tanh(bg.sys.g, 3.5, 0.5, He.Y/4); push!(defaults, Here2.He₊Y => He.Y)
    return ThermodynamicsSystem(bg, ODESystem[#=H, He, Hre, Here1, Here2=#]; defaults, kwargs...)
end

# TODO: make BaryonSystem or something, then merge into a background_baryon component?
# TODO: integrate using E/kB*T as independent variable?
# TODO: make e⁻ and γ species
function ThermodynamicsSystem(bg::BackgroundSystem, atoms::AbstractArray{ODESystem}; Xeϵ=0.0, defaults = Dict(), kwargs...)
    mHe_mH = 4.002603 / 1.00784 # ⪅ 4 ("not4") # TODO: move to Constants
    @parameters Tγ0 Yp fHe = Yp / (mHe_mH*(1-Yp)) # fHe = nHe/nH
    @variables Xe(η) ne(η) τ(η) = 0.0 dτ(η) ρb(η) Tγ(η) Tb(η) βb(η) cs²(η) λe(η)
    @variables XH⁺(η) nH(η) αH(η) βH(η) KH(η) CH(η) # H <-> H⁺
    @variables XHe⁺(η) nHe(η) αHe(η) βHe(η) KHe(η) KHe0⁻¹(η) KHe1⁻¹(η) KHe2⁻¹(η) γ2Ps(η) CHe(η) # He <-> He⁺
    @variables XHe⁺⁺(η) RHe⁺(η) AHe⁺(η) BHe⁺(η) CHe⁺(η) # He⁺ <-> He⁺⁺
    @variables αHe3(η) βHe3(η) τHe3(η) pHe3(η) CHe3(η) γ2Pt(η) # Helium triplet correction
    @variables DXHe⁺_Dη_singlet(η) DXHe⁺_Dη_triplet(η)
    push!(defaults, Tγ0 => (bg.sys.ph.ρ0 * 15/π^2 * bg.sys.g.H0^2/G * ħ^3*c^5)^(1/4) / kB) # TODO: make part of background species?

    # RECFAST implementation of Hydrogen and Helium recombination (https://arxiv.org/pdf/astro-ph/9909275 + https://arxiv.org/abs/astro-ph/9912182))
    ΛH = 8.2245809 # s⁻¹
    ΛHe = 51.3 # s⁻¹
    λ_H_∞_1s   =  91.17534471485433e-9 # n=∞ -> n=1 (ionization energy)
    λ_H_2s_1s  = 121.56700176979052e-9 # n=2 -> n=1
    λ_He_∞_1s  =  50.425904246895875e-9 # n=∞ -> n=1 (first ionization energy)
    λ_He_2s_1s =  60.14045177050301e-9 # n=2, l=1 -> n=1, l=1
    λ_He_2p_1s =  58.433437749406686e-9 # n=2, l=2 -> n=1, l=1
    E_H_∞_1s   = h*c / λ_H_∞_1s # E_∞  - E_1s
    E_H_2s_1s  = h*c / λ_H_2s_1s # E_2s - E_1s
    E_H_∞_2s   = E_H_∞_1s - E_H_2s_1s # E_∞ - E_2s
    E_He_∞_1s  = h*c / λ_He_∞_1s # E_∞  - E_1s
    E_He_2s_1s = h*c / λ_He_2s_1s # E_2s - E_1s
    E_He_2p_1s = h*c / λ_He_2p_1s # E_2p - E_1s
    E_He_2p_2s = E_He_2p_1s - E_He_2s_1s # E_2p - E_2s
    E_He_∞_2s  = E_He_∞_1s - E_He_2s_1s # E_∞ - E_2s
    E_He⁺_∞_1s = 54.4178 * eV # E_∞ - E_1s (second ionization energy)
    A2Ps = 1.798287e9
    A2Pt = 177.58e0
    λ_He_2p_1s_tri = 59.1411e-9 # TODO: rename s,t to singlet,triplet?
    λ_He_2s_1s_tri = 62.5563e-9
    E_He_2p_1s_tri = h*c / λ_He_2p_1s_tri
    E_He_2s_1s_tri = h*c / λ_He_2s_1s_tri
    E_He_2p_2s_tri = E_He_2p_1s_tri - E_He_2s_1s_tri
    λwtf = 260.0463e-9 # TODO: wtf is this?
    αH_fit(T; F=1.14, a=4.309, b=-0.6166, c=0.6703, d=0.5300, T0=1e4) = F * 1e-19 * a * (T/T0)^b / (1 + c * (T/T0)^d) # fitting formula to Hummer's table (fudge factor 1.14 here is equivalent to way RECFAST does it)
    αHe_fit(T; q=10^(-16.744), p=0.711, T1=10^5.114, T2=3.0) = q / (√(T/T2) * (1+√(T/T2))^(1-p) * (1+√(T/T1))^(1+p)) # fitting formula
    αHe3_fit(T; q=10^(-16.306), p=0.761, T1=10^5.114, T2=3.0) = q / (√(T/T2) * (1+√(T/T2))^(1-p) * (1+√(T/T1))^(1+p)) # fitting formula for Helium triplet correction
    KH_KH0_fit(a, A, z, w) = A*exp(-((log(a)+z)/w)^2)
    KH_KH0_fit(a) = KH_KH0_fit(a, -0.14, 7.28, 0.18) + KH_KH0_fit(a, 0.079, 6.73, 0.33)
    initialization_eqs = [
        XHe⁺ ~ 1, # TODO: add first order correction?
        XH⁺ ~ 1 - αH/βH, # + O((α/β)²); from solving β*(1-X) = α*X*Xe*n with Xe=X
        Tb ~ Tγ
    ]
    g = bg.sys.g
    connections = ODESystem([
        ρb ~ bg.sys.bar.ρ * g.H0^2/G # kg/m³ (convert from H0=1 units to SI units)
        nH ~ (1-Yp) * ρb/mp # 1/m³
        nHe ~ fHe * nH # 1/m³

        Tγ ~ Tγ0 / g.a # alternative derivative: Dη(Tγ) ~ -1*Tγ * g.ℰ
        Dη(Tb) ~ -2*Tb*g.ℰ - g.a/g.H0 * 8/3*σT*aR*Tγ^4 / (me*c) * Xe / (1+fHe+Xe) * (Tb-Tγ) # baryon temperature
        βb ~ 1 / (kB*Tb) # inverse temperature ("coldness")
        λe ~ h / √(2π*me/βb) # e⁻ de-Broglie wavelength
        cs² ~ kB/(mp*c^2) * (Tb - Dη(Tb)/g.ℰ) # https://arxiv.org/pdf/astro-ph/9506072 eq. (69) # TODO: proper mean molecular weight

        # H⁺ + e⁻ recombination
        αH ~ αH_fit(Tb)
        βH ~ αH / λe^3 * exp(-βb*E_H_∞_2s)
        KH ~ λ_H_2s_1s^3 / (8π*g.H) * (1 + KH_KH0_fit(g.a))
        CH ~ (1 + KH*ΛH*nH*(1-XH⁺)) /
             (1 + KH*(ΛH+βH)*nH*(1-XH⁺))
        Dη(XH⁺) ~ -g.a/g.H0 * CH * (αH*XH⁺*ne - βH*(1-XH⁺)*exp(-βb*E_H_2s_1s)) # XH⁺ = nH⁺ / nH; multiplied by H0 on left because cide η is physical η/(1/H0)

        # He⁺ + e⁻ singlet recombination
        αHe ~ αHe_fit(Tb)
        βHe ~ 4 * αHe / λe^3 * exp(-βb*E_He_∞_2s)
        KHe0⁻¹ ~ (8π*g.H) / λ_He_2p_1s^3 # RECFAST He flag 0
        KHe1⁻¹ ~ -exp(-3*A2Ps*nHe*(1-XHe⁺+1e-10)/KHe0⁻¹) * KHe0⁻¹ # RECFAST He flag 1 (close to zero modification?) # TODO: not that good, reliability depends on ϵ to avoid division by 0; try to use proper Saha ICs with XHe⁺ ≠ 1.0 and remove it
        γ2Ps ~ 3*A2Ps*fHe*(1-XHe⁺+1e-10)*c^2 / (1.436289e-22*8π*√(2π/(βb*mp*mHe_mH*c^2))*(1-XH⁺+1e-10)*(c/λ_He_2p_1s)^3) # TODO: introduce mHe and ν_He_2p_1s?
        KHe2⁻¹ ~ A2Ps/(1+0.36*γ2Ps^0.86)*3*nHe*(1-XHe⁺) # RECFAST He flag 2 (Doppler correction) # TODO: increase reliability, particularly at initial time
        KHe ~ 1 / (KHe0⁻¹ + KHe1⁻¹ + KHe2⁻¹) # corrections to inverse KHe are additive
        CHe ~ (exp(-βb*E_He_2p_2s) + KHe*ΛHe*nHe*(1-XHe⁺)) /
              (exp(-βb*E_He_2p_2s) + KHe*(ΛHe+βHe)*nHe*(1-XHe⁺))
        DXHe⁺_Dη_singlet ~ -g.a/g.H0 * CHe * (αHe*XHe⁺*ne - βHe*(1-XHe⁺)*exp(-βb*E_He_2s_1s))

        # He⁺ + e⁻ triplet recombination
        αHe3 ~ αHe3_fit(Tb)
        βHe3 ~ 4/3 * αHe3 / λe^3 * exp(-βb*h*c/λwtf)
        τHe3 ~ A2Pt*nHe*(1-XHe⁺+1e-10)*3 * λ_He_2p_1s_tri^3/(8π*g.H)
        pHe3 ~ (1 - exp(-τHe3)) / τHe3
        γ2Pt ~ 3*A2Pt*fHe*(1-XHe⁺+1e-10)*c^2 / (8π*1.484872e-22*c/λ_He_2p_1s_tri*√(2π/(βb*mp*mHe_mH*c^2))*(1-XH⁺+1e-10)) / (c/λ_He_2p_1s_tri)^2
        CHe3 ~ (1e-10 + A2Pt*(pHe3+1/(1+0.66*γ2Pt^0.9)/3)*exp(-βb*E_He_2p_2s_tri)) /
               (1e-10 + A2Pt*(pHe3+1/(1+0.66*γ2Pt^0.9)/3)*exp(-βb*E_He_2p_2s_tri) + βHe3) # added 1e-10 to avoid NaN at late times (does not change early behavior) # TODO: is sign in p-s exponentials wrong/different to what it is in just CHe?
        DXHe⁺_Dη_triplet ~ -g.a/g.H0 * CHe3 * (αHe3*XHe⁺*ne - βHe3*(1-XHe⁺)*3*exp(-βb*E_He_2s_1s_tri))

        # He⁺ + e⁻ total recombination
        Dη(XHe⁺) ~ DXHe⁺_Dη_singlet + DXHe⁺_Dη_triplet

        # He⁺⁺ + e⁻ recombination
        RHe⁺ ~ 1 * exp(-βb*E_He⁺_∞_1s) / (nH * λe^3) # right side of equation (6) in https://arxiv.org/pdf/astro-ph/9909275
        AHe⁺ ~ 1
        BHe⁺ ~ RHe⁺ - 1 - fHe
        CHe⁺ ~ -RHe⁺ * (1 + 2*fHe)
        XHe⁺⁺ ~ (-BHe⁺ + √(BHe⁺^2 - 4*AHe⁺*CHe⁺)) / (2*AHe⁺) - 1 - fHe # solve quadratic equation A*(XHe⁺⁺)² + B*XHe⁺⁺ + C = 0 from Saha equation (6) in https://arxiv.org/pdf/astro-ph/9909275; XHe⁺⁺ = Xe - 1*XH⁺ - fHe*XHe⁺ (with XH⁺ = XHe⁺ = 1)

        # electrons
        Xe ~ 1*XH⁺ + fHe*XHe⁺ + XHe⁺⁺ # TODO: redefine XHe⁺⁺ so it is also 1 at early times!
        ne ~ Xe * nH # TODO: redefine Xe = ne/nb ≠ ne/nH

        Dη(τ) ~ -g.a/g.H0 * ne * σT * c # common optical depth τ
        dτ ~ Dη(τ)
    ], η; defaults, initialization_eqs, kwargs...)
    sys = compose(connections, [atoms; bg.sys])
    ssys = structural_simplify(sys) # alternatively, disable simplifcation and construct "manually" to get helium Xe in the system
    return ThermodynamicsSystem(sys, ssys, bg)
end

function solve(th::ThermodynamicsSystem, Ωγ0, Ων0, Ωc0, Ωb0, h, Yp; aini=1e-4, aend=1.0, solver=Rodas5P(), reltol=1e-8, kwargs...)
    bg = th.bg
    bg_sol = solve(bg, Ωγ0, Ων0, Ωc0, Ωb0; aini, aend)
    ηini, ηtoday = bg_sol[η][begin], bg_sol[η][end]
    ΩΛ0 = bg_sol.ps[bg.ssys.de.Ω0]

    # TODO: use defaults for th.ssys.Xe => 1 + Yp/2
    prob = ODEProblem(th.ssys, [bg.ssys.g.a => aini], (ηini, ηtoday), [bg.sys.ph.Ω0 => Ωγ0, bg.sys.neu.Ω0 => Ων0, bg.sys.cdm.Ω0 => Ωc0, bg.sys.bar.Ω0 => Ωb0, bg.sys.de.Ω0 => ΩΛ0, bg.sys.g.H0 => H100 * h, th.ssys.Yp => Yp])

    # make solver take smaller steps when some quantity goes out of bounds: https://docs.sciml.ai/DiffEqDocs/stable/basics/faq/#My-ODE-goes-negative-but-should-stay-positive,-what-tools-can-help?
    #XHindex = variable_index(th.ssys, th.ssys.H.X)
    #isoutofdomain = (u, p, t) -> (XH = u[XHindex]; XH > 1.0)

    return solve(prob, solver; reltol, #=isoutofdomain,=# kwargs...) # CLASS uses "NDF15" (https://lesgourg.github.io/class-tour/London2014/Numerical_Methods_in_CLASS_London.pdf)
end