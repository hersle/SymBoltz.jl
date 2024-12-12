# TODO: add again
function reionization_tanh(g, z0, Δz0, Xe0; kwargs...)
    y(z) = (1+z)^(3/2)
    Δy(z, Δz0) = 3/2 * (1+z)^(1/2) * Δz0
    @variables Xe(t) z(t) T(t) ne(t) nb(t) nγ(t)
    @parameters Y
    @named base = ODESystem(Equation[], t, [T, ne, nb, nγ], [Y])
    return extend(ODESystem([
        z ~ 1/g.a - 1
        Xe ~ smoothifelse(y(z0)-y(z), 0, Xe0; k=1/Δy(z0, Δz0)) # smooth step from 0 to Xe0
    ], t; kwargs...), base)
end

# TODO: make BaryonSystem or something, then merge into a background_baryon component?
# TODO: make e⁻ and γ species
function thermodynamics_recombination_recfast(g; kwargs...)
    @parameters Yp fHe # fHe = nHe/nH
    @variables Xe(t) ne(t) τ(t) τ̇(t) τ̈(t) τ⃛(t) v(t) v̇(t) v̈(t) ρb(t) Tγ(t) Tb(t) DTb(t) βb(t) μ(t) cₛ²(t) λe(t)
    @variables XH⁺(t) nH(t) αH(t) βH(t) KH(t) KH0(t) KH1(t) CH(t) # H <-> H⁺
    @variables XHe⁺(t) nHe(t) αHe(t) βHe(t) KHe(t) KHe0⁻¹(t) KHe1⁻¹(t) KHe2⁻¹(t) γ2Ps(t) CHe(t) # He <-> He⁺
    @variables XHe⁺⁺(t) RHe⁺(t) # He⁺ <-> He⁺⁺
    @variables αHe3(t) βHe3(t) τHe3(t) pHe3(t) CHe3(t) γ2Pt(t) # Helium triplet correction
    @variables DXHe⁺_singlet(t) DXHe⁺_triplet(t)

    # RECFAST implementation of Hydrogen and Helium recombination (https://arxiv.org/pdf/astro-ph/9909275 + https://arxiv.org/abs/astro-ph/9912182))
    ΛH = 8.2245809 # s⁻¹
    ΛHe = 51.3 # s⁻¹
    A2Ps = 1.798287e9
    A2Pt = 177.58e0
    λwtf = 260.0463e-9; fwtf = c/λwtf; Ewtf = h*fwtf # TODO: wtf is this?
    αH_fit(T; F=1.14, a=4.309, b=-0.6166, c=0.6703, d=0.5300, T₀=1e4) = F * 1e-19 * a * (T/T₀)^b / (1 + c * (T/T₀)^d) # fitting formula to Hummer's table (fudge factor 1.14 here is equivalent to way RECFAST does it)
    αHe_fit(T, q, p, T1, T2) = q / (√(T/T2) * (1+√(T/T2))^(1-p) * (1+√(T/T1))^(1+p)) # fitting formula
    αHe_fit(T) = αHe_fit(T, 10^(-16.744), 0.711, 10^5.114, 3.0)
    αHe3_fit(T) = αHe_fit(T, 10^(-16.306), 0.761, 10^5.114, 3.0)
    KH_KH0_fit(a, A, z, w) = A*exp(-((log(a)+z)/w)^2)
    KH_KH0_fit(a) = KH_KH0_fit(a, -0.14, 7.28, 0.18) + KH_KH0_fit(a, 0.079, 6.73, 0.33)
    initialization_eqs = [
        XHe⁺ ~ 1 # TODO: add first order correction?
        XH⁺ ~ 1 - αH/βH # + O((α/β)²); from solving β*(1-X) = α*X*Xe*n with Xe=X
        Tb ~ Tγ
    ]
    defaults = [
        fHe => Yp / (mHe/mH*(1-Yp))
        τ => 0.0
    ]
    description = "Baryon-photon recombination thermodynamics (RECFAST)"
    return ODESystem([
        nH ~ (1-Yp) * ρb/mH # 1/m³
        nHe ~ fHe * nH # 1/m³

        D(Tb) ~ -2*Tb*g.ℰ - g.a/g.H₀ * 8/3*σT*aR*Tγ^4 / (me*c) * Xe / (1+fHe+Xe) * (Tb-Tγ) # baryon temperature
        DTb ~ D(Tb)
        βb ~ 1 / (kB*Tb) # inverse temperature ("coldness")
        λe ~ h / √(2π*me/βb) # e⁻ de-Broglie wavelength
        μ ~ mH / ((1 + (mH/mHe-1)*Yp + Xe*(1-Yp))) # mean molecular weight
        cₛ² ~ kB/(μ*c^2) * (Tb - D(Tb)/3g.ℰ) # https://arxiv.org/pdf/astro-ph/9506072 eq. (68)

        # H⁺ + e⁻ recombination
        αH ~ αH_fit(Tb)
        βH ~ αH / λe^3 * exp(-βb*E_H_∞_2s)
        KH0 ~ λ_H_2s_1s^3 / (8π*g.H)
        KH1 ~ KH0 * KH_KH0_fit(g.a)
        KH ~ KH0 + KH1
        CH ~ (1 + KH*ΛH*nH*(1-XH⁺+1e-10)) /
             (1 + KH*(ΛH+βH)*nH*(1-XH⁺+1e-10))
        D(XH⁺) ~ -g.a/g.H₀ * CH * (αH*XH⁺*ne - βH*(1-XH⁺)*exp(-βb*E_H_2s_1s)) # XH⁺ = nH⁺ / nH; multiplied by H₀ on left because side t is physical t/(1/H₀)

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
        DXHe⁺_singlet ~ -g.a/g.H₀ * CHe * (αHe*XHe⁺*ne - βHe*(1-XHe⁺)*exp(-βb*E_He_2s_1s))

        # He⁺ + e⁻ triplet recombination
        αHe3 ~ αHe3_fit(Tb)
        βHe3 ~ 4/3 * αHe3 / λe^3 * exp(-βb*Ewtf)
        τHe3 ~ A2Pt*nHe*(1-XHe⁺+1e-10)*3 * λ_He_2p_1s_tri^3/(8π*g.H)
        pHe3 ~ (1 - exp(-τHe3)) / τHe3
        γ2Pt ~ abs(3*A2Pt*fHe*(1-XHe⁺+1e-10)*c^2 / (8π*1.484872e-22*f_He_2p_1s_tri*√(2π/(βb*mHe*c^2))*(1-XH⁺+1e-10)) / (f_He_2p_1s_tri)^2) # abs to reduce chance for early-time crash # TODO: address properly
        CHe3 ~ (1e-10 + A2Pt*(pHe3+1/(1+0.66*γ2Pt^0.9)/3)*exp(-βb*E_He_2p_2s_tri)) /
               (1e-10 + A2Pt*(pHe3+1/(1+0.66*γ2Pt^0.9)/3)*exp(-βb*E_He_2p_2s_tri) + βHe3) # added 1e-10 to avoid NaN at late times (does not change early behavior) # TODO: is sign in p-s exponentials wrong/different to what it is in just CHe?
        DXHe⁺_triplet ~ -g.a/g.H₀ * CHe3 * (αHe3*XHe⁺*ne - βHe3*(1-XHe⁺)*3*exp(-βb*E_He_2s_1s_tri))

        # He⁺ + e⁻ total recombination
        D(XHe⁺) ~ DXHe⁺_singlet + DXHe⁺_triplet

        # He⁺⁺ + e⁻ recombination
        RHe⁺ ~ 1 * exp(-βb*E_He⁺_∞_1s) / (nH * λe^3) # right side of equation (6) in https://arxiv.org/pdf/astro-ph/9909275
        XHe⁺⁺ ~ 2*RHe⁺*fHe / (1+fHe+RHe⁺) / (1 + √(1 + 4*RHe⁺*fHe/(1+fHe+RHe⁺)^2)) # solve quadratic Saha equation (6) in https://arxiv.org/pdf/astro-ph/9909275 with the method of https://arxiv.org/pdf/1011.3758#equation.6.96

        # electrons
        Xe ~ 1*XH⁺ + fHe*XHe⁺ + XHe⁺⁺ # TODO: redefine XHe⁺⁺ so it is also 1 at early times!
        ne ~ Xe * nH # TODO: redefine Xe = ne/nb ≠ ne/nH

        τ̇ ~ -g.a/g.H₀ * ne * σT * c # common optical depth τ
        D(τ) ~ τ̇
        v ~ D(exp(-τ)) # visibility function
        v̇ ~ D(v)
        v̈ ~ D(v̇)
        τ̈ ~ D(τ̇) # TODO: unstable at thermodynamics stage
        τ⃛ ~ D(τ̈) # TODO: unstable at thermodynamics stage
    ], t, [ρb, Xe, XH⁺, XHe⁺, XHe⁺⁺, τ, τ̇, τ̈, τ⃛, v, v̇, v̈, Tb, Tγ, μ, cₛ²], [Yp, fHe]; initialization_eqs, defaults, description, kwargs...)
end

function thermodynamics_ΛCDM(bg::ODESystem; spline=false, kwargs...)
    if spline
        @named rec = thermodynamics_recombination_splined(bg.g)
        eqs = []
    else
        @named rec = thermodynamics_recombination_recfast(bg.g)
        eqs = [rec.ρb ~ bg.bar.ρ * bg.g.H₀^2/GN, rec.Tγ ~ bg.ph.T] # kg/m³ (convert from H₀=1 units to SI units)
    end
    th = ODESystem(eqs, t; kwargs...)
    return compose(th, rec, bg)
end

function thermodynamics_recombination_splined(; kwargs...)
    dummyspline = CubicSpline([NaN, NaN, NaN], 0.0:1.0:2.0)
    vars = @variables τ(t) τ̇(t) τ̈(t) τ⃛(t) cₛ²(t) v(t) v̇(t) v̈(t) #Tb(t)
    pars = @parameters τspline::CubicSpline τ̇spline::CubicSpline cₛ²spline::CubicSpline #Tbspline::CubicSpline
    defaults = [τspline => dummyspline, τ̇spline => dummyspline, cₛ²spline => dummyspline]
    return ODESystem([
        τ ~ value(τspline, t)
        τ̇ ~ derivative(τspline, t)
        τ̈ ~ derivative(τspline, t, 2) # TODO: unstable?
        τ⃛ ~ derivative(τ̇spline, t, 2) # DataInterpolations doesn't support 3rd order derivatives, so spline τ̇ and take its 2nd order derivative # TODO: remove workaround
        v ~ -τ̇ * exp(-τ) # TODO: do automatically
        v̇ ~ (-τ̈ + τ̇^2) * exp(-τ) # TODO: automatic
        v̈ ~ (τ⃛ + 3*τ̇*τ̈ - τ̇^3) * exp(-τ) # TODO: automatic
        cₛ² ~ value(cₛ²spline, t)
        #Tb ~ exp(Tbspline(log(t)))
    ], t, vars, pars; defaults, kwargs...) # connect perturbation τ with spline evaluation
end
