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

        XH⁺(τ) = 1.0, [description = "H ionization fraction n(H⁺)/nH"] # TODO: add first order correction?
        nH(τ), [description = "Total H number density"]
        αH(τ), βH(τ), KH(τ), KHfitfactor(τ), CH(τ)

        nHe(τ), [description = "Total He number density"]
        XHe⁺(τ) = 1.0, [description = "Singly ionized He fraction"] # - αH/βH # + O((α/β)²); from solving β*(1-X) = α*X*Xe*n with Xe=X
        XHe⁺⁺(τ), [description = "Doubly ionized He fraction"]
        αHe(τ), βHe(τ), RHe⁺(τ), τHe(τ), KHe(τ), invKHe0(τ), invKHe1(τ), invKHe2(τ), CHe(τ), DXHe⁺(τ), DXHet⁺(τ) # invK = 1 / K
    end

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
        CH ~ smoothifelse(XH⁺ - XlimC, (1 + KH*ΛH2s1s*nH*(1-XH⁺)) / (1 + KH*(ΛH2s1s+βH)*nH*(1-XH⁺)), 1; k = 1e3) # CLASS has FH in denominator; SymBoltz has it in αH (similar to Rdown in CLASS)
        D(XH⁺) ~ -g.a/(H100*g.h) * CH * (αH*XH⁺*ne - βH*(1-XH⁺)*exp(-β*EH2s1s)) # XH⁺ = nH⁺ / nH; multiplied by H₀ on left because side τ is physical τ/(1/H₀)

        # He⁺ + e⁻ singlet recombination
        αHe ~ αHefit(T; q=10^(-16.744), p=0.711)
        βHe ~ 4 * αHe / λe^3 * exp(-β*EHe∞2s)
        KHe ~ 1 / (invKHe0 + invKHe1 + invKHe2) # corrections are additive in inverse KHe
        invKHe0 ~ 8π*H / λHe2p1s^3
        CHe ~ smoothifelse(XHe⁺ - XlimC, (exp(-β*EHe2p2s) + KHe*ΛHe2s1s*nHe*(1-XHe⁺)) / (exp(-β*EHe2p2s) + KHe*(ΛHe2s1s+βHe)*nHe*(1-XHe⁺)), 1; k = 1e3) # TODO: normal ifelse()? https://github.com/SciML/ModelingToolkit.jl/issues/3897
        DXHe⁺ ~ -g.a/(H100*g.h) * CHe * (αHe*XHe⁺*ne - βHe*(1-XHe⁺)*exp(-β*EHe2s1s))

        # He⁺ + e⁻ total recombination
        D(XHe⁺) ~ DXHe⁺ + DXHet⁺ # singlet + triplet

        # He⁺⁺ + e⁻ recombination
        RHe⁺ ~ 1 * exp(-β*EHe⁺∞1s) / (nH * λe^3) # right side of equation (6) in https://arxiv.org/pdf/astro-ph/9909275
        XHe⁺⁺ ~ 2*RHe⁺*fHe / (1+fHe+RHe⁺) / (1 + √(1 + 4*RHe⁺*fHe/(1+fHe+RHe⁺)^2)) # solve quadratic Saha equation (6) in https://arxiv.org/pdf/astro-ph/9909275 with the method of https://arxiv.org/pdf/1011.3758#equation.6.96

        Xe ~ 1*XH⁺ + fHe*XHe⁺ + XHe⁺⁺ # TODO: redefine XHe⁺⁺ so it is also 1 at early times?
    ]

    ics = Dict()
    if Hswitch == 0
        push!(ics, FH => 1.14) # original fudge factor
        push!(eqs, KHfitfactor ~ 1)
    elseif Hswitch == 1
        push!(ics, FH => 1.125) # fudged fudge factor in RECFAST 1.5.2 to match new He physics in https://arxiv.org/abs/1110.0247
        KHfitfactorfunc(a, A, z, w) = A*exp(-((log(a)+z)/w)^2) # Gaussian fit in log(a)-space
        push!(eqs, KHfitfactor ~ 1 + KHfitfactorfunc(g.a, -0.14, 7.28, 0.18) + KHfitfactorfunc(g.a, 0.079, 6.73, 0.33))
    else
        error("Supported H switches are 0 and 1. Got $Hswitch.")
    end

    if Heswitch == 0 # no corrections
        append!(eqs, [DXHet⁺ ~ 0, invKHe1 ~ 0, invKHe2 ~ 0])
    elseif Heswitch == 6 # all corrections (Doppler, triplet etc.)
        # RECFAST switches off He corrections when XH⁺ ≈ XHe⁺ ≈ 1, but we use a smooth+symmetric regularization of
        # (1-X) that makes it a small positive number, even if numerical errors causes X to drift slightly above 1
        reg(x; ϵ = 1e-9) = √(x^2 + ϵ^2) # regularize x so it remains small and positive even if x→0 or drifts to x<0
        γHe(; A=NaN, σ=NaN, f=NaN) = 3*A*fHe*reg(1-XHe⁺)*c^2 / (8π*σ*√(2π/(β*mHe*c^2))*reg(1-XH⁺)*f^3)
        append!(vars, @variables γ2ps(τ) αHet(τ) βHet(τ) τHet(τ) pHet(τ) CHet(τ) CHetnum(τ) γ2pt(τ))
        append!(eqs, [
            τHe ~ 3*AHe2p1s*nHe*reg(1-XHe⁺) / invKHe0
            invKHe1 ~ -exp(-τHe) * invKHe0 # RECFAST He flag 1

            γ2ps ~ γHe(A = AHe2p1s, σ = σHe2p1s, f = fHe2p1s)
            invKHe2 ~ AHe2p1s/(1+0.36*γ2ps^0.86)*3*nHe*(1-XHe⁺) # RECFAST He flag 2 (Doppler correction)

            # He⁺ + e⁻ triplet recombination
            αHet ~ αHefit(T; q=10^(-16.306), p=0.761)
            βHet ~ 4/3 * αHet / λe^3 * exp(-β*EHet∞2s)
            τHet ~ AHet2p1s*nHe*reg(1-XHe⁺)*3 * λHet2p1s^3/(8π*H)
            pHet ~ (1 - exp(-τHet)) / τHet
            γ2pt ~ γHe(A = AHet2p1s, σ = σHet2p1s, f = fHet2p1s)
            CHetnum ~ AHet2p1s*(pHet+1/(1+0.66*γ2pt^0.9)/3)*exp(-β*EHet2p2s) # numerator of CHet
            CHet ~ reg(CHetnum) / (reg(CHetnum) + βHet) # TODO: is sign in p-s exponentials wrong/different to what it is in just CHe?
            DXHet⁺ ~ -g.a/(H100*g.h) * CHet * (αHet*XHe⁺*ne - βHet*(1-XHe⁺)*3*exp(-β*EHet2s1s))
        ])
    else
        error("Supported He switches are 0 and 6. Got $Heswitch.") # TODO support more granular switches 1-5?
    end
    description = "Baryon-photon recombination thermodynamics (RECFAST)"
    return System(eqs, τ, vars, pars; initial_conditions = ics, description, kwargs...)
end

# HyRec2 rate tables (https://github.com/nanoomlee/HYREC-2), interpolated with cubic B-splines
# Extrapolate flatly above Tr = 0.4 eV (H is in Saha equilibrium for any rates) and linearly below Tr = 0.004 eV (z ≲ 16)
const hyrec_dir = joinpath(@__DIR__, "..", "..", "data", "hyrec")
hyrec_read(file) = parse.(Float64, split(read(joinpath(hyrec_dir, file), String)))
hyrec_spline(y, extrap, x...) = Interpolations.extrapolate(Interpolations.scale(Interpolations.interpolate(y, Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid())))), x...), extrap)
const hyrec_lnTr = range(log(0.004), log(0.4), length = 100) # ln(Tr/eV)
const hyrec_TmTr = range(0.1, 1.0, length = 40) # Tm/Tr
const hyrec_lnα_tables = let α = reshape(hyrec_read("Alpha_inf.dat"), 4, 40, 100) # α2s, α2p (Tm<Tr), α2s, α2p (Tm>Tr) × Tm/Tr × Tr
    [hyrec_spline(log.(α[i,:,:]'), ((Interpolations.Line(), Interpolations.Flat()), Interpolations.Line()), hyrec_lnTr, hyrec_TmTr) for i in 1:2]
end
const hyrec_lnR_table = hyrec_spline(log.(hyrec_read("R_inf.dat")), ((Interpolations.Line(), Interpolations.Flat()),), hyrec_lnTr)
const hyrec_Δ_tables = let Δ = reshape(hyrec_read("fit_swift.dat"), 5, :) # Tr/K, Δ(fid), ∂Δ/∂ωcb, ∂Δ/∂ωH, ∂Δ/∂Neff
    [hyrec_spline(Δ[i,:], Interpolations.Flat(), range(Δ[1,1], Δ[1,end], length = size(Δ, 2))) for i in 2:5]
end

# ln(αᵢ/(cm³/s)) for i = 1 (2s) and 2 (2p) as function of ln(Tr/eV) and Tm/Tr.
hyrec_lnα(i, lnTr, TmTr) = hyrec_lnα_tables[i](lnTr, TmTr)
hyrec_lnα_grad(i, j, lnTr, TmTr) = j == 1 ? ForwardDiff.derivative(x -> hyrec_lnα(i, x, TmTr), lnTr) : ForwardDiff.derivative(x -> hyrec_lnα(i, lnTr, x), TmTr) # Interpolations.gradient fails with mixed extrapolation
# ln(R2p2s/(1/s)) as function of ln(Tr/eV).
hyrec_lnR(lnTr) = hyrec_lnR_table(lnTr)
hyrec_lnR_grad(lnTr) = ForwardDiff.derivative(hyrec_lnR, lnTr)
# SWIFT correction function (i = 1) and its derivatives wrt. ωcb, ωH, Neff (i = 2, 3, 4) as function of Tr/K.
hyrec_Δ(i, Tr) = hyrec_Δ_tables[i](Tr)
hyrec_Δ_grad(i, Tr) = ForwardDiff.derivative(x -> hyrec_Δ(i, x), Tr)
@register_symbolic hyrec_lnα(i, lnTr, TmTr)
@register_symbolic hyrec_lnα_grad(i, j, lnTr, TmTr)
@register_symbolic hyrec_lnR(lnTr)
@register_symbolic hyrec_lnR_grad(lnTr)
@register_symbolic hyrec_Δ(i, Tr)
@register_symbolic hyrec_Δ_grad(i, Tr)
@register_derivative hyrec_lnα(i, lnTr, TmTr) 2 hyrec_lnα_grad(i, 1, lnTr, TmTr)
@register_derivative hyrec_lnα(i, lnTr, TmTr) 3 hyrec_lnα_grad(i, 2, lnTr, TmTr)
@register_derivative hyrec_lnR(lnTr) 1 hyrec_lnR_grad(lnTr)
@register_derivative hyrec_Δ(i, Tr) 2 hyrec_Δ_grad(i, Tr)

"""
    recombination_hyrec(g, YHe, fHe; kwargs...)

Recombination physics for Hydrogen and Helium based on HyRec2 in its default SWIFT mode.
Hydrogen follows the effective multi-level atom with 2s and 2p states (EMLA2s2p) and tabulated effective rates,
with the Lyman-α escape rate corrected by the SWIFT function Δ fitted to full radiative transfer calculations.
Helium follows HyRec2's He II → He I equation (no fudge factors) and Saha equilibrium for He III → He II.
Internal quantities use HyRec2's units (cm, s, K or eV).

References
==========
- https://arxiv.org/abs/2007.14114 (HyRec2)
- https://arxiv.org/abs/1011.3758 (HyRec)
- https://github.com/nanoomlee/HYREC-2
"""
function recombination_hyrec(g, YHe, fHe; kwargs...)
    pars = @parameters begin
        ωcb, [description = "Reduced baryon and cold dark matter density parameter (for SWIFT correction only)"]
        Neff, [description = "Effective number of neutrinos (for SWIFT correction only)"]
    end
    vars = @variables begin
        Xe(τ), [description = "Free electron fraction contribution"]
        ne(τ), [description = "Total free electron number density"] # 1/m³
        nH(τ), [description = "Total H number density"] # 1/m³
        nHe(τ), [description = "Total He number density"] # 1/m³
        T(τ), [description = "Baryon temperature"] # K
        Tγ(τ), [description = "Photon temperature"] # K
        H(τ), [description = "Cosmic Hubble function in 1/s"]
        n(τ), [description = "Total H number density in 1/cm³"]
        xe(τ), [description = "Total free electron fraction ne/nH"]

        XH⁺(τ) = 1.0, [description = "H ionization fraction n(H⁺)/nH"]
        x1s(τ), [description = "Neutral H fraction n(H 1s)/nH"]
        x1sreg(τ), [description = "Neutral H fraction regularized to stay positive"]
        Tr(τ), [description = "Photon temperature in eV"]
        Tm(τ), [description = "Baryon temperature in eV"]
        α2s(τ), α2p(τ), [description = "Effective recombination coefficients to 2s and 2p in cm³/s"]
        Dα2s(τ), Dα2p(τ), [description = "Effective recombination coefficients relative to Tm = Tr in cm³/s"]
        β2s(τ), β2p(τ), [description = "Effective photoionization rates from 2s and 2p in 1/s"]
        R2p2s(τ), [description = "Effective 2p → 2s transfer rate in 1/s"]
        Δ(τ), [description = "SWIFT correction to the Lyman-α escape rate"]
        RLyax1s(τ), [description = "Lyman-α escape rate multiplied by x1s in 1/s"]
        Γ2s(τ), Γ2px1s(τ), [description = "Inverse lifetimes of 2s and 2p (latter multiplied by x1s) in 1/s"]
        C2s(τ), C2p(τ), [description = "Generalized Peebles C factors for 2s and 2p"]
        sH(τ), [description = "H Saha ratio xe*xHII/x1s"]

        XHe⁺(τ) = 1.0, [description = "Singly ionized He fraction"]
        XHe⁺⁺(τ), [description = "Doubly ionized He fraction"]
        xHeI(τ), [description = "Neutral He fraction n(He I)/nH"]
        xHeII(τ), [description = "Singly ionized He fraction n(He II)/nH"]
        s0(τ), sHe(τ), sHe⁺(τ), [description = "He Saha ratios"]
        ηc(τ), [description = "H continuum opacity parameter in s (1/etacinv in HyRec2)"]
        Γ2pinc(τ), [description = "Incoherent width of He 2¹P in 1/s"]
        τ2p(τ), [description = "He 2¹P → 1¹S Sobolev optical depth"]
        Δνline(τ), τc(τ), enh(τ), [description = "He continuum opacity escape enhancement"]
        pesc(τ), [description = "He 2¹P escape probability multiplied by exp(-6989 K/Tr)"]
        ydown(τ), [description = "He I recombination rate coefficient in 1/s"]
    end

    kBeV = 8.617343e-5 # eV/K
    EI = 13.598286071938324 # H ionization energy (eV)
    SAHA = 3.016103031869581e21 # (2π μe/h²)^(3/2) in eV^(-3/2) cm⁻³
    LYA = 4.662899067555897e15 # 8π/(3λLyα³) in cm⁻³
    L2s1s = 8.2206 # 2s → 1s two-photon decay rate (1/s)
    reg(x; ϵ = 1e-5) = √(x^2 + ϵ^2) # regularize x so it remains positive; ϵ ≫ solver tolerance so the solver sees smooth rates in Saha equilibrium (x → 0)
    P(τ) = (1 - exp(-τ)) / τ # Sobolev escape probability (→ 1/τ for τ ≫ 1)
    nB(E) = exp(-E/Tγ) / (1 - exp(-E/Tγ)) # photon occupation number 1/(exp(E/Tγ)-1) for energy E in K (written to avoid overflow)
    αHeB(T; q=10^(-10.744), p=0.711, T1=10^5.114, T2=3.0) = q / (√(T/T2) * (1+√(T/T2))^(1-p) * (1+√(T/T1))^(1+p)) # He I case B recombination coefficient in cm³/s
    q3 = (2.7255 / (Tγ * g.a))^3 # SWIFT fiducial (T₀fid / T₀)³
    ωH = nH * mH * g.a^3 * 8π*GN / (3*H100^2) # reduced H density parameter

    eqs = [
        H ~ H100 * g.h * g.H
        n ~ nH / 1e6
        xe ~ ne / nH
        Tr ~ kBeV * Tγ
        Tm ~ kBeV * T

        # H⁺ + e⁻ recombination (rec_swift_hyrec_dxHIIdlna in HyRec2, multiplied through by x1s to avoid 1/x1s)
        x1s ~ 1 - XH⁺
        x1sreg ~ reg(x1s) # use in rate coefficients (optically thick Lyα escape ∝ 1/x1s makes rates nearly singular as x1s → 0), but keep exact x1s in net rates
        α2s ~ exp(hyrec_lnα(1, log(Tr), Tm/Tr))
        α2p ~ exp(hyrec_lnα(2, log(Tr), Tm/Tr))
        Dα2s ~ α2s - exp(hyrec_lnα(1, log(Tr), 1.0))
        Dα2p ~ α2p - exp(hyrec_lnα(2, log(Tr), 1.0))
        β2s ~ (α2s - Dα2s) * SAHA * Tr^(3/2) * exp(-EI/4Tr) # detailed balance with α(Tm = Tr)
        β2p ~ (α2p - Dα2p) * SAHA * Tr^(3/2) * exp(-EI/4Tr) / 3
        R2p2s ~ exp(hyrec_lnR(log(Tr)))
        Δ ~ hyrec_Δ(1, Tγ) + (ωcb - 0.14175)*q3 * hyrec_Δ(2, Tγ) + (ωH - 0.02242*(1-0.246738546372))*q3 * hyrec_Δ(3, Tγ) + (Neff - 3.046) * hyrec_Δ(4, Tγ)
        RLyax1s ~ LYA * H / n / (1 + Δ)
        Γ2s ~ β2s + 3*R2p2s + L2s1s
        Γ2px1s ~ (β2p + R2p2s) * x1sreg + RLyax1s
        C2s ~ (L2s1s + 3*R2p2s*RLyax1s/Γ2px1s) / (Γ2s - 3*R2p2s^2*x1sreg/Γ2px1s)
        C2p ~ (RLyax1s + x1sreg*R2p2s*L2s1s/Γ2s) / (Γ2px1s - 3*R2p2s^2*x1sreg/Γ2s)
        sH ~ SAHA * Tr^(3/2) * exp(-EI/Tr) / n
        D(XH⁺) ~ -g.ℋ * n/H * ((sH*x1s*Dα2s + α2s*(xe*XH⁺ - sH*x1s))*C2s + (sH*x1s*Dα2p + α2p*(xe*XH⁺ - sH*x1s))*C2p) # dlna/dτ = ℋ

        # He⁺ + e⁻ recombination (rec_helium_dxHeIIdlna in HyRec2)
        xHeII ~ fHe * XHe⁺
        xHeI ~ fHe * (1 - XHe⁺)
        s0 ~ 2.414194e15 * Tγ^(3/2) / n * 4
        sHe ~ s0 * exp(-285325/Tγ)
        ηc ~ n * reg(x1s; ϵ = 1e-12) / (9.15776e22 * H) # H continuum opacity is sensitive to small x1s
        Γ2pinc ~ 1.976e6 * (1 + nB(6989)) + 6.03e6 * nB(19754) + 1.06e8 * nB(21539) + 2.18e6 * nB(28496) + 3.37e7 * nB(29224) + 1.04e6 * nB(32414) + 1.51e7 * nB(32781)
        τ2p ~ 4.277e-8 * n/H * fHe * reg(1 - XHe⁺) # optically thick escape ∝ 1/τ2p makes rate nearly singular as XHe⁺ → 1
        Δνline ~ Γ2pinc * τ2p / (4π^2)
        τc ~ Δνline * ηc
        enh ~ √(1 + π^2*τc) + 7.74*τc/(1 + 70*τc)
        pesc ~ enh * P(τ2p) * exp(-6989/Tγ) + (1 - exp(-1.023e-7*τ2p)) / τ2p * (0.964525*exp(-4042/Tγ) - enh*exp(-6.14e13*ηc - 6989/Tγ)) # exp(-6989/Tr) absorbed to avoid overflow at late times
        ydown ~ 1 / (s0 * exp(-46090/Tγ) / (50.94 + 3*1.7989e9*pesc) + 1 / (1e3 * n * αHeB(T))) # HyRec2's rate, capped at 10³ × case B recombination to avoid overflow at late times (HyRec2 stops He then)
        D(XHe⁺) ~ g.ℋ / fHe * ydown * (xHeI*sHe - xHeII*xe) / H

        # He⁺⁺ + e⁻ recombination (Saha equilibrium; rec_xesaha_HeII_III in HyRec2)
        sHe⁺ ~ 2.414194e15 * Tγ^(3/2) * exp(-631462.7/Tγ) / n
        XHe⁺⁺ ~ 2*sHe⁺*fHe / (1+sHe⁺+fHe) / (1 + √(1 + 4*sHe⁺*fHe/(1+sHe⁺+fHe)^2))

        Xe ~ XH⁺ + fHe*XHe⁺ + XHe⁺⁺
    ]
    description = "Baryon-photon recombination thermodynamics (HyRec2)"
    return System(eqs, τ, vars, pars; description, kwargs...)
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
The `recombination` model is `:recfast` (or `true`; with options `Hswitch` and `Heswitch`), `:hyrec` or `false` (none).
"""
function baryons(g; recombination = true, reionization = true, Hswitch = 1, Heswitch = 6, name = :b, kwargs...)
    description = "Baryonic matter"
    b = matter(g; adiabatic = false, interact = true, continuity_pressure = false, name, description, kwargs...) |> complete

    pars = @parameters begin
        YHe, [description = "Primordial He abundance or mass fraction ρ(He)/(ρ(H)+ρ(He))"]
        fHe = YHe / (mHe/mH*(1-YHe)), [description = "Primordial He/H nucleon ratio n(He)/n(H)"] # fHe = nHe/nH
    end
    vars = @variables begin
        κ(τ) = 0.0, [backwards = true, description = "Optical depth (0 today, so integrate it backwards)"]
        κ̇(τ), [description = "Optical depth derivative"]
        I(τ), [description = "Optical depth exponential exp(-κ)"]
        v(τ), [description = "Visibility function"]
        v̇(τ), [description = "Visibility function derivative"]
        cₛ²(τ), [description = "Thermal speed of sound squared (from gas temperature; different from other species)"]
        T(τ), [description = "Baryon temperature"]
        Tγ(τ), [description = "Photon temperature"]
        ΔT(τ) = 0.0, [description = "Baryon-photon temperature difference"] # Tb ≈ Tγ at early times
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
        κ̇ ~ -g.a/(H100*g.h) * ne * σT * c # optical depth derivative
        D(κ) ~ κ̇
        I ~ exp(-κ)
        v ~ D(exp(-κ)) |> expand_derivatives # visibility function
        v̇ ~ D(v)
        cₛ² ~ kB/μc² * (T - D(T)/3g.ℋ) # thermal (adiabatic) speed of sound Ṗ/ρ̇ of the gas from https://arxiv.org/pdf/astro-ph/9506072 eq. 68; different from other species' cₛ²
        b.δP ~ cₛ² * b.ρ * b.δ # pressure perturbation in Euler equation (not continuity equation) and sourcing the total in gravity (see https://arxiv.org/pdf/astro-ph/9506072 eq. 67) # TODO: consistent thermal w = kB*T/μc² instead (needs time-dependent w and background+thermodynamics coupling)
        μc² ~ mH*c^2 / (1 + (mH/mHe-1)*YHe + Xe*(1-YHe))

        DT ~ -2*T*g.ℋ - g.a/g.h * 8/3*σT*aR/H100*Tγ^4 / (me*c) * Xe / (1+fHe+Xe) * ΔT # baryon temperature
        DTγ ~ D(Tγ) # or -1*Tγ*g.ℋ
        D(ΔT) ~ DT - DTγ # solve ODE for D(T-Tγ), since solving it for D(T) instead is extremely sensitive to T-Tγ≈0 at early times
        T ~ ΔT + Tγ

        nH ~ (1-YHe) * b.ρ*(H100*g.h)^2/GN / mH # 1/m³; convert b.ρ from H₀=1 units to SI units
        nHe ~ fHe * nH # 1/m³
        ne ~ Xe * nH # TODO: redefine Xe = ne/nb ≠ ne/nH?
    ]

    if recombination == true || recombination == :recfast
        @named rec = recombination_recfast(g, ParentScope(YHe), ParentScope(fHe); Hswitch, Heswitch)
        push!(eqs, rec.nH ~ nH, rec.nHe ~ nHe, rec.ne ~ ne, rec.T ~ T)
        push!(comps, rec)
    elseif recombination == :hyrec
        @named rec = recombination_hyrec(g, ParentScope(YHe), ParentScope(fHe))
        push!(eqs, rec.nH ~ nH, rec.nHe ~ nHe, rec.ne ~ ne, rec.T ~ T, rec.Tγ ~ Tγ)
        push!(comps, rec)
    elseif recombination != false
        error("Unknown recombination model $recombination")
    end

    if reionization
        @named rei1 = reionization_tanh(g, 7.6711, 0.5, 3//2, 1 + ParentScope(fHe))
        @named rei2 = reionization_tanh(g, 3.5, 0.5, 1, ParentScope(fHe))
        push!(comps, rei1, rei2)
    end

    push!(eqs, Xe ~ sum(comp.Xe for comp in comps; init = 0))

    b = extend(b, System(eqs, τ, vars, pars; name); description)
    b = compose(b, comps)
    return b
end
