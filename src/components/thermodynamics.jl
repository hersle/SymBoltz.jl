# TODO: make BaryonSystem or something, then merge into a background_baryon component?
# TODO: make e⁻ and γ species
function thermodynamics_recombination_recfast(g; reionization = true, kwargs...)
    pars = @parameters begin
        Yp, [description = "Primordial He mass fraction ρ(He)/(ρ(H)+ρ(He))"] # TODO: correct?
        fHe, [description = "Primordial He-to-H nucleon ratio n(He)/n(H)"] # TODO: correct?
        re1z, [description = "1st reionization redshift"]
        re2z, [description = "2nd reionization redshift"]
        κ0, [description = "Optical depth today"] # to make the real κ = 0 today
    end
    vars = @variables begin
        Xe(τ), [description = "Total free electron fraction"]
        ne(τ), [description = "Free electron number density"]
        λe(τ), [description = "Electron de-Broglie wavelength"]

        κ(τ), [description = "Optical depth normalized to 0 today"]
        _κ(τ), [description = "Optical depth normalized to 0 initially"]
        κ̇(τ), [description = "Optical depth derivative"]
        v(τ), [description = "Visibility function"]
        v̇(τ), [description = "Visibility function derivative"]

        Tγ(τ), [description = "Photon temperature"]
        Tb(τ), [description = "Baryon temperature"]
        ΔT(τ), [description = "Difference between baryon and photon temperature"]
        DTγ(τ), [description = "Photon temperature derivative"]
        DTb(τ), [description = "Baryon temperature derivative"]
        βb(τ), [description = "inverse temperature (coldness)"]
        cₛ²(τ), [description = "Speed of sound squared"]

        ρb(τ), [description = "Baryon density"]
        μ(τ), [description = "Mean molecular weight"]

        re1Xe(τ), [description = "1st reionization free electron fraction contribution"]
        re2Xe(τ), [description = "2nd reionization free electron fraction contribution"]

        XH⁺(τ), [description = "Hydrogen ionization factor"] # H <-> H⁺
        nH(τ), [description = "Total Hydrogen number density"]
        αH(τ), βH(τ), KH(τ), KH0(τ), KH1(τ), CH(τ)

        nHe(τ), [description = "Total Helium number density"]
        XHe⁺(τ), [description = "Singly ionized Helium fraction"] # He <-> He⁺
        XHe⁺⁺(τ), [description = "Doubly ionized Helium fraction"] # He⁺ <-> He⁺⁺
        αHe(τ), βHe(τ), KHe(τ), KHe0⁻¹(τ), KHe1⁻¹(τ), KHe2⁻¹(τ), γ2Ps(τ), CHe(τ) # K(...)⁻¹ = 1 / K(...)
        RHe⁺(τ), αHe3(τ), βHe3(τ), τHe3(τ), pHe3(τ), CHe3(τ), γ2Pt(τ), DXHe⁺_singlet(τ), DXHe⁺_triplet(τ) # He triplet correction variables
    end

    # RECFAST implementation of Hydrogen and Helium recombination (https://arxiv.org/pdf/astro-ph/9909275 + https://arxiv.org/abs/astro-ph/9912182))
    ΛH = 8.2245809 # s⁻¹
    ΛHe = 51.3 # s⁻¹
    A2Ps = 1.798287e9
    A2Pt = 177.58e0
    αH_fit(T; F=1.14, a=4.309, b=-0.6166, c=0.6703, d=0.5300, T₀=1e4) = F * 1e-19 * a * (T/T₀)^b / (1 + c * (T/T₀)^d) # fitting formula to Hummer's table (fudge factor 1.14 here is equivalent to way RECFAST does it)
    αHe_fit(T, q, p, T1, T2) = q / (√(T/T2) * (1+√(T/T2))^(1-p) * (1+√(T/T1))^(1+p)) # fitting formula
    αHe_fit(T) = αHe_fit(T, 10^(-16.744), 0.711, 10^5.114, 3.0)
    αHe3_fit(T) = αHe_fit(T, 10^(-16.306), 0.761, 10^5.114, 3.0)
    KH_KH0_fit(a, A, z, w) = A*exp(-((log(a)+z)/w)^2)
    KH_KH0_fit(a) = KH_KH0_fit(a, -0.14, 7.28, 0.18) + KH_KH0_fit(a, 0.079, 6.73, 0.33)
    defaults = [
        XHe⁺ => 1.0 # TODO: add first order correction?
        XH⁺ => 1.0 - αH/βH # + O((α/β)²); from solving β*(1-X) = α*X*Xe*n with Xe=X
        fHe => Yp / (mHe/mH*(1-Yp)) # fHe = nHe/nH # TODO: factor mHe/mH?
        _κ => 0.0
        κ0 => NaN
        re1z => 7.6711
        re2z => 3.5
        ΔT => 0.0 # i.e. Tb ~ Tγ at early times
    ]
    description = "Baryon-photon recombination thermodynamics (RECFAST)"

    # reionization utility functions
    y(z) = (1+z)^(3/2)
    Δy(z, Δz0) = 3/2 * (1+z)^(1/2) * Δz0
    smoothifelse(x, v1, v2; k=1) = 1/2 * ((v1+v2) + (v2-v1)*tanh(k*x)) # smooth transition/step function from v1 at x<0 to v2 at x>0

    return ODESystem([
        nH ~ (1-Yp) * ρb/mH # 1/m³
        nHe ~ fHe * nH # 1/m³

        DTb ~ -2*Tb*g.ℰ - g.a/g.h * 8/3*σT*aR/H100*Tγ^4 / (me*c) * Xe / (1+fHe+Xe) * ΔT # baryon temperature
        DTγ ~ D(Tγ) # or -1*Tγ*g.ℰ
        D(ΔT) ~ DTb - DTγ # solve ODE for D(Tb-Tγ), since solving it for D(Tb) instead is extremely sensitive to Tb-Tγ≈0 at early times
        Tb ~ ΔT + Tγ
        βb ~ 1 / (kB*Tb)
        λe ~ h / √(2π*me/βb) # e⁻ de-Broglie wavelength
        μ ~ mH / ((1 + (mH/mHe-1)*Yp + Xe*(1-Yp)))
        cₛ² ~ kB/(μ*c^2) * (Tb - D(Tb)/3g.ℰ) # https://arxiv.org/pdf/astro-ph/9506072 eq. (68)

        # H⁺ + e⁻ recombination
        αH ~ αH_fit(Tb)
        βH ~ αH / λe^3 * exp(-βb*E_H_∞_2s)
        KH0 ~ λ_H_2s_1s^3 / (8π*g.H)
        KH1 ~ KH0 * KH_KH0_fit(g.a)
        KH ~ KH0 + KH1
        CH ~ (1 + KH*ΛH*nH*(1-XH⁺+1e-10)) /
             (1 + KH*(ΛH+βH)*nH*(1-XH⁺+1e-10))
        D(XH⁺) ~ -g.a/(H100*g.h) * CH * (αH*XH⁺*ne - βH*(1-XH⁺)*exp(-βb*E_H_2s_1s)) # XH⁺ = nH⁺ / nH; multiplied by H₀ on left because side τ is physical τ/(1/H₀)

        # He⁺ + e⁻ singlet recombination
        αHe ~ αHe_fit(Tb)
        βHe ~ 4 * αHe / λe^3 * exp(-βb*E_He_∞_2s)
        KHe0⁻¹ ~ (8π*g.H) / λ_He_2p_1s^3 # RECFAST He flag 0
        KHe1⁻¹ ~ -exp(-3*A2Ps*nHe*(1-XHe⁺+1e-10)/KHe0⁻¹) * KHe0⁻¹ # RECFAST He flag 1 (close to zero modification?) # TODO: not that good, reliability depends on ϵ to avoid division by 0; try to use proper Saha ICs with XHe⁺ ≠ 1.0 and remove it
        γ2Ps ~ abs(3*A2Ps*fHe*(1-XHe⁺+1e-10)*c^2 / (1.436289e-22*8π*√(2π/(βb*mHe*c^2))*(1-XH⁺+1e-10)*(f_He_2p_1s)^3)) # abs to reduce chance for early-time crash # TODO: address properly # TODO: introduce ν_He_2p_1s?
        KHe2⁻¹ ~ A2Ps/(1+0.36*γ2Ps^0.86)*3*nHe*(1-XHe⁺) # RECFAST He flag 2 (Doppler correction) # TODO: increase reliability, particularly at initial time
        KHe ~ 1 / (KHe0⁻¹ + KHe1⁻¹ + KHe2⁻¹) # corrections to inverse KHe are additive
        CHe ~ (exp(-βb*E_He_2p_2s) + KHe*ΛHe*nHe*(1-XHe⁺)) /
              (exp(-βb*E_He_2p_2s) + KHe*(ΛHe+βHe)*nHe*(1-XHe⁺))
        DXHe⁺_singlet ~ -g.a/(H100*g.h) * CHe * (αHe*XHe⁺*ne - βHe*(1-XHe⁺)*exp(-βb*E_He_2s_1s))

        # He⁺ + e⁻ triplet recombination
        αHe3 ~ αHe3_fit(Tb)
        βHe3 ~ 4/3 * αHe3 / λe^3 * exp(-βb*E_He_∞_2s_tri)
        τHe3 ~ A2Pt*nHe*(1-XHe⁺+1e-10)*3 * λ_He_2p_1s_tri^3/(8π*g.H)
        pHe3 ~ (1 - exp(-τHe3)) / τHe3
        γ2Pt ~ abs(3*A2Pt*fHe*(1-XHe⁺+1e-10)*c^2 / (8π*1.484872e-22*f_He_2p_1s_tri*√(2π/(βb*mHe*c^2))*(1-XH⁺+1e-10)) / (f_He_2p_1s_tri)^2) # abs to reduce chance for early-time crash # TODO: address properly
        CHe3 ~ (1e-10 + A2Pt*(pHe3+1/(1+0.66*γ2Pt^0.9)/3)*exp(-βb*E_He_2p_2s_tri)) /
               (1e-10 + A2Pt*(pHe3+1/(1+0.66*γ2Pt^0.9)/3)*exp(-βb*E_He_2p_2s_tri) + βHe3) # added 1e-10 to avoid NaN at late times (does not change early behavior) # TODO: is sign in p-s exponentials wrong/different to what it is in just CHe?
        DXHe⁺_triplet ~ -g.a/(H100*g.h) * CHe3 * (αHe3*XHe⁺*ne - βHe3*(1-XHe⁺)*3*exp(-βb*E_He_2s_1s_tri))

        # He⁺ + e⁻ total recombination
        D(XHe⁺) ~ DXHe⁺_singlet + DXHe⁺_triplet

        # He⁺⁺ + e⁻ recombination
        RHe⁺ ~ 1 * exp(-βb*E_He⁺_∞_1s) / (nH * λe^3) # right side of equation (6) in https://arxiv.org/pdf/astro-ph/9909275
        XHe⁺⁺ ~ 2*RHe⁺*fHe / (1+fHe+RHe⁺) / (1 + √(1 + 4*RHe⁺*fHe/(1+fHe+RHe⁺)^2)) # solve quadratic Saha equation (6) in https://arxiv.org/pdf/astro-ph/9909275 with the method of https://arxiv.org/pdf/1011.3758#equation.6.96

        # reionization
        re1Xe ~ reionization ? smoothifelse(y(re1z)-y(g.z), 0, 1+fHe; k=1/Δy(re1z, 0.5)) : 0 # 1st reionization: H⁺ and He⁺ simultaneously
        re2Xe ~ reionization ? smoothifelse(y(re2z)-y(g.z), 0, fHe; k=1/Δy(re2z, 0.5)) : 0 # 2nd reionization: He⁺⁺

        # electrons
        Xe ~ 1*XH⁺ + fHe*XHe⁺ + XHe⁺⁺ + re1Xe + re2Xe # TODO: redefine XHe⁺⁺ so it is also 1 at early times!
        ne ~ Xe * nH # TODO: redefine Xe = ne/nb ≠ ne/nH

        D(_κ) ~ -g.a/(H100*g.h) * ne * σT * c # optical depth derivative
        κ̇ ~ D(_κ) # optical depth derivative
        κ ~ _κ - κ0 # optical depth offset such that κ = 0 today (non-NaN only after integration)

        v ~ D(exp(-κ)) |> expand_derivatives # visibility function
        v̇ ~ D(v)
    ], τ, vars, pars; defaults, description, kwargs...)
end
