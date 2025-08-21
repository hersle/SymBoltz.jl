"""
    thermodynamics_recombination_recfast(g; reionization = true, Hswitch = 1, Heswitch = 6, kwargs...)

Recombination physics for Hydrogen and Helium (including fudge factors) based on RECFAST 1.5.2.

References
==========
- https://www.astro.ubc.ca/people/scott/recfast.html
- https://www.astro.ubc.ca/people/scott/recfast.for
- https://arxiv.org/abs/astro-ph/9909275
- https://arxiv.org/abs/astro-ph/9912182).
- https://arxiv.org/abs/1110.0247
"""
function thermodynamics_recombination_recfast(g; reionization = true, Hswitch = 1, Heswitch = 6, kwargs...)
    pars = @parameters begin
        YHe, [description = "Primordial He abundance or mass fraction ρ(He)/(ρ(H)+ρ(He))"]
        fHe, [description = "Primordial He/H nucleon ratio n(He)/n(H)"]
        κ0 = NaN, [description = "Optical depth today (set retrospectively)"] # to make the real κ = 0 today
        XlimC = 0.99, [description = "Set CHe = CH = 1 for larger XH⁺ and XHe⁺ to avoid instabilities"]
        FH, [description = "Hydrogen fudge factor"] # to emulate more accurate and expensive multi-level calculation (https://arxiv.org/pdf/astro-ph/9909275)
    end
    vars = @variables begin
        Xe(τ), [description = "Total free electron fraction"]
        ne(τ), [description = "Free electron number density"]
        λe(τ), [description = "Electron de-Broglie wavelength"]

        κ(τ), [description = "Optical depth normalized to 0 today"]
        _κ(τ), [description = "Optical depth normalized to 0 initially"]
        κ̇(τ), [description = "Optical depth derivative"]
        I(τ), [description = "Optical depth exponential exp(-κ)"]
        v(τ), [description = "Visibility function"]
        v̇(τ), [description = "Visibility function derivative"]

        Tγ(τ), [description = "Photon temperature"]
        Tb(τ), [description = "Baryon temperature"]
        ΔT(τ) = 0, [description = "Baryon-photon temperature difference"] # Tb ≈ Tγ at early times
        DTγ(τ), [description = "Photon temperature derivative"]
        DTb(τ), [description = "Baryon temperature derivative"]
        βb(τ), [description = "Baryon inverse temperature (coldness)"]
        cₛ²(τ), [description = "Thermal speed of sound squared"]

        ρb(τ), [description = "Baryon mass density"]
        μc²(τ), [description = "Mean molecular weight multiplied by speed of light squared"]

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
    me = PhysicalConstants.CODATA2018.m_e / u"kg"
    mH = elements[:H].atomic_mass / u"kg" |> NoUnits
    mHe = elements[:He].atomic_mass / u"kg" |> NoUnits

    # Hydrogen singlet transitions
    λH∞1s   =  91.17534e-9; fH∞1s  = c/λH∞1s;  EH∞1s  = h*fH∞1s # ∞ - 1s (read: wavelength Hydrogen ∞ to 1s)
    λH2s1s  = 121.56700e-9; fH2s1s = c/λH2s1s; EH2s1s = h*fH2s1s # 2s - 1s
    EH∞2s  = EH∞1s - EH2s1s # E(∞) - E(2s)

    # Helium singlet transitions
    λHe∞1s  =  50.42590e-9; fHe∞1s  = c/λHe∞1s;  EHe∞1s  = h*fHe∞1s
    λHe2s1s =  60.14045e-9; fHe2s1s = c/λHe2s1s; EHe2s1s = h*fHe2s1s
    λHe2p1s =  58.43344e-9; fHe2p1s = c/λHe2p1s; EHe2p1s = h*fHe2p1s
    EHe2p2s = EHe2p1s - EHe2s1s
    EHe∞2s  = EHe∞1s - EHe2s1s
    EHe⁺∞1s = 54.4178 * eV

    defaults = [
        XHe⁺ => 1.0 # TODO: add first order correction?
        XH⁺ => 1.0 # - αH/βH # + O((α/β)²); from solving β*(1-X) = α*X*Xe*n with Xe=X
        _κ => 0.0
    ]

    αHfit(T; F=FH, a=4.309, b=-0.6166, c=0.6703, d=0.5300, T₀=1e4) = F * 1e-19 * a * (T/T₀)^b / (1 + c * (T/T₀)^d) # fitting formula to Hummer's table (fudge factor here is equivalent to the way RECFAST does it)
    αHefit(T; q=NaN, p=NaN, T1=10^5.114, T2=3.0) = q / (√(T/T2) * (1+√(T/T2))^(1-p) * (1+√(T/T1))^(1+p)) # fitting formula

    eqs = [
        # parameter equations
        fHe ~ YHe / (mHe/mH*(1-YHe)) # fHe = nHe/nH

        nH ~ (1-YHe) * ρb/mH # 1/m³
        nHe ~ fHe * nH # 1/m³

        DTb ~ -2*Tb*g.ℰ - g.a/g.h * 8/3*σT*aR/H100*Tγ^4 / (me*c) * Xe / (1+fHe+Xe) * ΔT # baryon temperature
        DTγ ~ D(Tγ) # or -1*Tγ*g.ℰ
        D(ΔT) ~ DTb - DTγ # solve ODE for D(Tb-Tγ), since solving it for D(Tb) instead is extremely sensitive to Tb-Tγ≈0 at early times
        Tb ~ ΔT + Tγ
        βb ~ 1 / (kB*Tb)
        λe ~ h / √(2π*me/βb) # e⁻ de-Broglie wavelength
        μc² ~ mH*c^2 / ((1 + (mH/mHe-1)*YHe + Xe*(1-YHe)))
        cₛ² ~ kB/μc² * (Tb - D(Tb)/3g.ℰ) # https://arxiv.org/pdf/astro-ph/9506072 eq. (68)

        # H⁺ + e⁻ recombination
        αH ~ αHfit(Tb)
        βH ~ αH / λe^3 * exp(-βb*EH∞2s)
        KH ~ KHfitfactor/8π * λH2s1s^3 / g.H # KHfitfactor ≈ 1; see above
        CH ~ smoothifelse(XH⁺ - XlimC, (1 + KH*ΛH*nH*(1-XH⁺)) / (1 + KH*(ΛH+βH)*nH*(1-XH⁺)), 1; k = 1e3) # CLASS has FH in denominator; SymBoltz has it in αH (similar to Rdown in CLASS)
        D(XH⁺) ~ -g.a/(H100*g.h) * CH * (αH*XH⁺*ne - βH*(1-XH⁺)*exp(-βb*EH2s1s)) # XH⁺ = nH⁺ / nH; multiplied by H₀ on left because side τ is physical τ/(1/H₀)

        # He⁺ + e⁻ singlet recombination
        αHe ~ αHefit(Tb; q=10^(-16.744), p=0.711)
        βHe ~ 4 * αHe / λe^3 * exp(-βb*EHe∞2s)
        KHe ~ 1 / (invKHe0 + invKHe1 + invKHe2) # corrections are additive in inverse KHe
        invKHe0 ~ 8π*g.H / λHe2p1s^3
        CHe ~ smoothifelse(XHe⁺ - XlimC, (exp(-βb*EHe2p2s) + KHe*ΛHe*nHe*(1-XHe⁺)) / (exp(-βb*EHe2p2s) + KHe*(ΛHe+βHe)*nHe*(1-XHe⁺)), 1; k = 1e3) # TODO: normal ifelse()? https://github.com/SciML/ModelingToolkit.jl/issues/3897
        DXHe⁺ ~ -g.a/(H100*g.h) * CHe * (αHe*XHe⁺*ne - βHe*(1-XHe⁺)*exp(-βb*EHe2s1s))

        # He⁺ + e⁻ total recombination
        D(XHe⁺) ~ DXHe⁺ + DXHet⁺ # singlet + triplet

        # He⁺⁺ + e⁻ recombination
        RHe⁺ ~ 1 * exp(-βb*EHe⁺∞1s) / (nH * λe^3) # right side of equation (6) in https://arxiv.org/pdf/astro-ph/9909275
        XHe⁺⁺ ~ 2*RHe⁺*fHe / (1+fHe+RHe⁺) / (1 + √(1 + 4*RHe⁺*fHe/(1+fHe+RHe⁺)^2)) # solve quadratic Saha equation (6) in https://arxiv.org/pdf/astro-ph/9909275 with the method of https://arxiv.org/pdf/1011.3758#equation.6.96

        D(_κ) ~ -g.a/(H100*g.h) * ne * σT * c # optical depth derivative
        κ̇ ~ D(_κ) # optical depth derivative
        κ ~ _κ - κ0 # optical depth offset such that κ = 0 today (non-NaN only after integration)
        I ~ exp(-κ)

        v ~ D(exp(-κ)) |> expand_derivatives # visibility function
        v̇ ~ D(v)
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
        λHet∞2s = 260.0463e-9; fHet∞2s = c/λHet∞2s; EHet∞2s = h*fHet∞2s # ∞ - 2³s; ionization of lowest triplet state (4.77 or 4.8 eV) (read: "wavelength Helium triplet ∞ to 2s")
        λHet2p1s = 59.1411e-9; fHet2p1s = c/λHet2p1s; EHet2p1s = h*fHet2p1s
        λHet2s1s = 62.5563e-9; fHet2s1s = c/λHet2s1s; EHet2s1s = h*fHet2s1s
        EHet2p2s = EHet2p1s - EHet2s1s
        A2ps = 1.798287e9 # A 2p singlet
        A2pt = 177.58e0 # A 2p triplet
        γHe(; A=NaN, σ=NaN, f=NaN) = 3*A*fHe*(1-XHe⁺+ϵ)*c^2 / (8π*σ*√(2π/(βb*mHe*c^2))*(1-XH⁺+ϵ)*f^3)
        append!(vars, @variables γ2ps(τ) αHet(τ) βHet(τ) τHet(τ) pHet(τ) CHet(τ) CHetnum(τ) γ2pt(τ))
        append!(eqs, [
            τHe ~ 3*A2ps*nHe*(1-XHe⁺+ϵ) / invKHe0
            invKHe1 ~ -exp(-τHe) * invKHe0 # RECFAST He flag 1

            γ2ps ~ γHe(A = A2ps, σ = 1.436289e-22, f = fHe2p1s)
            invKHe2 ~ A2ps/(1+0.36*γ2ps^0.86)*3*nHe*(1-XHe⁺) # RECFAST He flag 2 (Doppler correction)

            # He⁺ + e⁻ triplet recombination
            αHet ~ αHefit(Tb; q=10^(-16.306), p=0.761)
            βHet ~ 4/3 * αHet / λe^3 * exp(-βb*EHet∞2s)
            τHet ~ A2pt*nHe*(1-XHe⁺+ϵ)*3 * λHet2p1s^3/(8π*g.H)
            pHet ~ (1 - exp(-τHet)) / τHet
            γ2pt ~ γHe(A = A2pt, σ = 1.484872e-22, f = fHet2p1s)
            CHetnum ~ A2pt*(pHet+1/(1+0.66*γ2pt^0.9)/3)*exp(-βb*EHet2p2s) # numerator of CHet
            CHet ~ (ϵ + CHetnum) / (ϵ + CHetnum + βHet) # TODO: is sign in p-s exponentials wrong/different to what it is in just CHe?
            DXHet⁺ ~ -g.a/(H100*g.h) * CHet * (αHet*XHe⁺*ne - βHet*(1-XHe⁺)*3*exp(-βb*EHet2s1s))
        ])
    else
        error("Supported He switches are 0 and 6. Got $Heswitch.") # TODO support more granular switches 1-5?
    end

    if reionization
        @named rei1 = reionization_tanh(g)
        @named rei2 = reionization_tanh(g)
        append!(defaults, [
            rei1.z => 7.6711, rei1.Δz => 0.5, rei1.n => 3/2,
            rei2.z => 3.5, rei2.Δz => 0.5, rei2.n => 1,
        ])
        append!(eqs, [
            rei1.Xemax ~ 1 + fHe
            rei2.Xemax ~ fHe
        ])
        reis = [rei1, rei2]
    else
        reis = []
    end

    append!(eqs, [
        Xe ~ 1*XH⁺ + fHe*XHe⁺ + XHe⁺⁺ + sum(rei.Xe for rei in reis; init = 0) # total free electron fraction # TODO: redefine XHe⁺⁺ so it is also 1 at early times?
        ne ~ Xe * nH # TODO: redefine Xe = ne/nb ≠ ne/nH?
    ])

    description = "Baryon-photon recombination thermodynamics (RECFAST)"
    rec = System(eqs, τ, vars, pars; defaults, description, kwargs...)
    rec = compose(rec, reis)
    return rec
end

"""
    reionization_tanh(g; kwargs...)

Reionization physics with a free electron function that is activated around a given redshift `z` by a tanh function.
Equivalent to the reionization model in CAMB.

References
==========
- https://cosmologist.info/notes/CAMB.pdf#section*.10
"""
function reionization_tanh(g; kwargs...)
    pars = @parameters begin
        z, [description = "reionization redshift"]
        Δz, [description = "reionization redshift width"]
        n, [description = "reionization exponent"]
        Xemax, [description = "fully ionized contribution to the free electron fraction"]
    end
    vars = @variables begin
        Xe(τ), [description = "free electron fraction contribution"]
    end
    eqs = [
        Xe ~ smoothifelse((1+z)^n - (1+g.z)^n, 0, Xemax; k = 1/(n*(1+z)^(n-1)*Δz))
    ]
    description = "Reionization with tanh-like (activation function) contribution to the free electron fraction"
    return System(eqs, τ, vars, pars; description, kwargs...)
end
