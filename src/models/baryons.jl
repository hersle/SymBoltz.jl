"""
    recombination_recfast(g; reionization = true, Hswitch = 1, Heswitch = 6, kwargs...)

Recombination physics for Hydrogen and Helium (including fudge factors) based on RECFAST 1.5.2.

References
==========
- https://www.astro.ubc.ca/people/scott/recfast.html
- https://www.astro.ubc.ca/people/scott/recfast.for
- https://arxiv.org/abs/astro-ph/9909275
- https://arxiv.org/abs/astro-ph/9912182).
- https://arxiv.org/abs/1110.0247
"""
function recombination_recfast(g, YHe, fHe; reionization = true, Hswitch = 1, Heswitch = 6, kwargs...)
    pars = @parameters begin
        XlimC = 0.99, [description = "Set CHe = CH = 1 for larger XH⁺ and XHe⁺ to avoid instabilities"]
        FH, [description = "Hydrogen fudge factor"] # to emulate more accurate and expensive multi-level calculation (https://arxiv.org/pdf/astro-ph/9909275)
    end
    vars = @variables begin
        Xe(τ), [description = "Free electron fraction contribution"]
        ne(τ), [description = "Total free electron number density"]
        λe(τ), [description = "Electron de-Broglie wavelength"]
        H(τ), [description = "Cosmic Hubble function in SI units"]

        T(τ), [description = "Temperature"]
        β(τ), [description = "Inverse temperature (coldness)"]

        XH⁺(τ), [description = "H ionization fraction n(H⁺)/nH"]
        nH(τ), [description = "Total H number density"]
        αH(τ), βH(τ), KH(τ), KHfitfactor(τ), CH(τ)

        nHe(τ), [description = "Total He number density"]
        XHe⁺(τ), [description = "Singly ionized He fraction"]
        XHe⁺⁺(τ), [description = "Doubly ionized He fraction"]
        αHe(τ), βHe(τ), RHe⁺(τ), τHe(τ), KHe(τ), invKHe0(τ), invKHe1(τ), invKHe2(τ), CHe(τ), DXHe⁺(τ), DXHet⁺(τ) # invK = 1 / K
    end

    ΛH = 8.2245809 # s⁻¹
    ΛHe = 51.3 # s⁻¹

    defaults = [
        XHe⁺ => 1.0 # TODO: add first order correction?
        XH⁺ => 1.0 # - αH/βH # + O((α/β)²); from solving β*(1-X) = α*X*Xe*n with Xe=X
    ]

    αHfit(T; F=FH, a=4.309, b=-0.6166, c=0.6703, d=0.5300, T₀=1e4) = F * 1e-19 * a * (T/T₀)^b / (1 + c * (T/T₀)^d) # fitting formula to Hummer's table (fudge factor here is equivalent to the way RECFAST does it)
    αHefit(T; q=NaN, p=NaN, T1=10^5.114, T2=3.0) = q / (√(T/T2) * (1+√(T/T2))^(1-p) * (1+√(T/T1))^(1+p)) # fitting formula

    eqs = [
        β ~ 1 / (kB*T)
        λe ~ h / √(2π*me/β) # e⁻ de-Broglie wavelength
        H ~ H100 * g.h * g.H

        # H⁺ + e⁻ recombination
        αH ~ αHfit(T)
        βH ~ αH / λe^3 * exp(-β*EH∞2s)
        KH ~ KHfitfactor/8π * λH2s1s^3 / H # KHfitfactor ≈ 1; see above
        CH ~ smoothifelse(XH⁺ - XlimC, (1 + KH*ΛH*nH*(1-XH⁺)) / (1 + KH*(ΛH+βH)*nH*(1-XH⁺)), 1; k = 1e3) # CLASS has FH in denominator; SymBoltz has it in αH (similar to Rdown in CLASS)
        D(XH⁺) ~ -g.a/(H100*g.h) * CH * (αH*XH⁺*ne - βH*(1-XH⁺)*exp(-β*EH2s1s)) # XH⁺ = nH⁺ / nH; multiplied by H₀ on left because side τ is physical τ/(1/H₀)

        # He⁺ + e⁻ singlet recombination
        αHe ~ αHefit(T; q=10^(-16.744), p=0.711)
        βHe ~ 4 * αHe / λe^3 * exp(-β*EHe∞2s)
        KHe ~ 1 / (invKHe0 + invKHe1 + invKHe2) # corrections are additive in inverse KHe
        invKHe0 ~ 8π*H / λHe2p1s^3
        CHe ~ smoothifelse(XHe⁺ - XlimC, (exp(-β*EHe2p2s) + KHe*ΛHe*nHe*(1-XHe⁺)) / (exp(-β*EHe2p2s) + KHe*(ΛHe+βHe)*nHe*(1-XHe⁺)), 1; k = 1e3) # TODO: normal ifelse()? https://github.com/SciML/ModelingToolkit.jl/issues/3897
        DXHe⁺ ~ -g.a/(H100*g.h) * CHe * (αHe*XHe⁺*ne - βHe*(1-XHe⁺)*exp(-β*EHe2s1s))

        # He⁺ + e⁻ total recombination
        D(XHe⁺) ~ DXHe⁺ + DXHet⁺ # singlet + triplet

        # He⁺⁺ + e⁻ recombination
        RHe⁺ ~ 1 * exp(-β*EHe⁺∞1s) / (nH * λe^3) # right side of equation (6) in https://arxiv.org/pdf/astro-ph/9909275
        XHe⁺⁺ ~ 2*RHe⁺*fHe / (1+fHe+RHe⁺) / (1 + √(1 + 4*RHe⁺*fHe/(1+fHe+RHe⁺)^2)) # solve quadratic Saha equation (6) in https://arxiv.org/pdf/astro-ph/9909275 with the method of https://arxiv.org/pdf/1011.3758#equation.6.96

        Xe ~ 1*XH⁺ + fHe*XHe⁺ + XHe⁺⁺ # TODO: redefine XHe⁺⁺ so it is also 1 at early times?
    ]

    if Hswitch == 0
        push!(defaults, FH => 1.14) # original fudge factor
        push!(eqs, KHfitfactor ~ 1)
    elseif Hswitch == 1
        push!(defaults, FH => 1.125) # fudged fudge factor in RECFAST 1.5.2 to match new He physics in https://arxiv.org/abs/1110.0247
        KHfitfactorfunc(a, A, z, w) = A*exp(-((log(a)+z)/w)^2) # Gaussian fit in log(a)-space
        push!(eqs, KHfitfactor ~ 1 + KHfitfactorfunc(g.a, -0.14, 7.28, 0.18) + KHfitfactorfunc(g.a, 0.079, 6.73, 0.33))
    else
        error("Supported H switches are 0 and 1. Got $Hswitch.")
    end

    if Heswitch == 0 # no corrections
        append!(eqs, [DXHet⁺ ~ 0, invKHe1 ~ 0, invKHe2 ~ 0])
    elseif Heswitch == 6 # all corrections (Doppler, triplet etc.)
        ϵ = 1e-9 # original RECFAST switches off He corrections when XHe⁺ ≈ 1.0; I add a tiny ϵ to avoid numerical instabilities
        A2ps = 1.798287e9 # A 2p singlet
        A2pt = 177.58e0 # A 2p triplet
        γHe(; A=NaN, σ=NaN, f=NaN) = 3*A*fHe*(1-XHe⁺+ϵ)*c^2 / (8π*σ*√(2π/(β*mHe*c^2))*(1-XH⁺+ϵ)*f^3)
        append!(vars, @variables γ2ps(τ) αHet(τ) βHet(τ) τHet(τ) pHet(τ) CHet(τ) CHetnum(τ) γ2pt(τ))
        append!(eqs, [
            τHe ~ 3*A2ps*nHe*(1-XHe⁺+ϵ) / invKHe0
            invKHe1 ~ -exp(-τHe) * invKHe0 # RECFAST He flag 1

            γ2ps ~ γHe(A = A2ps, σ = 1.436289e-22, f = fHe2p1s)
            invKHe2 ~ A2ps/(1+0.36*γ2ps^0.86)*3*nHe*(1-XHe⁺) # RECFAST He flag 2 (Doppler correction)

            # He⁺ + e⁻ triplet recombination
            αHet ~ αHefit(T; q=10^(-16.306), p=0.761)
            βHet ~ 4/3 * αHet / λe^3 * exp(-β*EHet∞2s)
            τHet ~ A2pt*nHe*(1-XHe⁺+ϵ)*3 * λHet2p1s^3/(8π*H)
            pHet ~ (1 - exp(-τHet)) / τHet
            γ2pt ~ γHe(A = A2pt, σ = 1.484872e-22, f = fHet2p1s)
            CHetnum ~ A2pt*(pHet+1/(1+0.66*γ2pt^0.9)/3)*exp(-β*EHet2p2s) # numerator of CHet
            CHet ~ (ϵ + CHetnum) / (ϵ + CHetnum + βHet) # TODO: is sign in p-s exponentials wrong/different to what it is in just CHe?
            DXHet⁺ ~ -g.a/(H100*g.h) * CHet * (αHet*XHe⁺*ne - βHet*(1-XHe⁺)*3*exp(-β*EHet2s1s))
        ])
    else
        error("Supported He switches are 0 and 6. Got $Heswitch.") # TODO support more granular switches 1-5?
    end
    description = "Baryon-photon recombination thermodynamics (RECFAST)"
    return System(eqs, τ, vars, pars; defaults, description, kwargs...)
end

"""
    reionization_tanh(g; kwargs...)

Reionization physics with a free electron function that is activated around a given redshift `z` by a tanh function.
Equivalent to the reionization model in CAMB.

References
==========
- https://cosmologist.info/notes/CAMB.pdf#section*.10
"""
function reionization_tanh(g, z, Δz, n, Xemax; kwargs...)
    vars = @variables begin
        Xe(τ), [description = "free electron fraction contribution"]
    end
    f(_z) = n % 1 == 1//2 ? √(1+_z) * (1+_z)^Int(n-1//2) : (1+_z)^n
    eqs = [
        Xe ~ smoothifelse(f(z) - f(g.z), 0, Xemax; k = 1/(n*(1+z)^(n-1)*Δz))
    ]
    description = "Reionization with tanh-like (activation function) contribution to the free electron fraction"
    return System(eqs, τ, vars, []; description, kwargs...)
end

"""
    baryons(g; recombination = true, reionization = true, Hswitch = 1, Heswitch = 6, name = :b, kwargs...)

Create a particle species for baryons in the spacetime with metric `g`.
"""
function baryons(g; recombination = true, reionization = true, Hswitch = 1, Heswitch = 6, name = :b, kwargs...)
    description = "Baryonic matter"
    b = matter(g; adiabatic = false, θinteract = true, name, description, kwargs...) |> complete

    pars = @parameters begin
        YHe, [description = "Primordial He abundance or mass fraction ρ(He)/(ρ(H)+ρ(He))"]
        fHe, [description = "Primordial He/H nucleon ratio n(He)/n(H)"]
        κ0 = NaN, [description = "Optical depth today (set retrospectively)"] # to make the real κ = 0 today
    end
    vars = @variables begin
        κ(τ), [description = "Optical depth normalized to 0 today"]
        _κ(τ), [description = "Optical depth normalized to 0 initially"]
        κ̇(τ), [description = "Optical depth derivative"]
        I(τ), [description = "Optical depth exponential exp(-κ)"]
        v(τ), [description = "Visibility function"]
        v̇(τ), [description = "Visibility function derivative"]
        cₛ²(τ), [description = "Thermal speed of sound squared"]
        T(τ), [description = "Baryon temperature"]
        Tγ(τ), [description = "Photon temperature"]
        ΔT(τ) = 0, [description = "Baryon-photon temperature difference"] # Tb ≈ Tγ at early times
        DTγ(τ), [description = "Photon temperature derivative"]
        DT(τ), [description = "Baryon temperature derivative"]
        μc²(τ), [description = "Mean molecular weight multiplied by speed of light squared"]
        Xe(τ), [description = "Total free electron fraction"]
        nH(τ), [description = "Total H number density"]
        nHe(τ), [description = "Total He number density"]
        ne(τ), [description = "Total free electron number density"]
    end

    comps = []
    eqs = [
        # parameter equations
        fHe ~ YHe / (mHe/mH*(1-YHe)) # fHe = nHe/nH

        D(_κ) ~ -g.a/(H100*g.h) * ne * σT * c # optical depth derivative
        κ̇ ~ D(_κ) # optical depth derivative
        κ ~ _κ - κ0 # optical depth offset such that κ = 0 today (non-NaN only after integration)
        I ~ exp(-κ)
        v ~ D(exp(-κ)) |> expand_derivatives # visibility function
        v̇ ~ D(v)
        cₛ² ~ kB/μc² * (T - D(T)/3g.ℋ) # https://arxiv.org/pdf/astro-ph/9506072 eq. (68)
        μc² ~ mH*c^2 / (1 + (mH/mHe-1)*YHe + Xe*(1-YHe))

        DT ~ -2*T*g.ℋ - g.a/g.h * 8/3*σT*aR/H100*Tγ^4 / (me*c) * Xe / (1+fHe+Xe) * ΔT # baryon temperature
        DTγ ~ D(Tγ) # or -1*Tγ*g.ℋ
        D(ΔT) ~ DT - DTγ # solve ODE for D(T-Tγ), since solving it for D(T) instead is extremely sensitive to T-Tγ≈0 at early times
        T ~ ΔT + Tγ

        nH ~ (1-YHe) * b.ρ*(H100*g.h)^2/GN / mH # 1/m³; convert b.ρ from H₀=1 units to SI units
        nHe ~ fHe * nH # 1/m³
        ne ~ Xe * nH # TODO: redefine Xe = ne/nb ≠ ne/nH?
    ]
    defaults = [
        _κ => 0.0
    ]

    if recombination
        @named rec = recombination_recfast(g, ParentScope(YHe), ParentScope(fHe); Hswitch, Heswitch)
        push!(eqs, rec.nH ~ nH, rec.nHe ~ nHe, rec.ne ~ ne, rec.T ~ T)
        push!(comps, rec)
    end

    if reionization
        @named rei1 = reionization_tanh(g, 7.6711, 0.5, 3//2, 1 + ParentScope(fHe))
        @named rei2 = reionization_tanh(g, 3.5, 0.5, 1, ParentScope(fHe))
        push!(comps, rei1, rei2)
    end

    push!(eqs, Xe ~ sum(comp.Xe for comp in comps; init = 0))

    b = extend(b, System(eqs, τ, vars, pars; defaults, name); description)
    b = compose(b, comps)
    return b
end
