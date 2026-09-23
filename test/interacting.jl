# General interacting baryon - cold dark matter - dark energy model
# Pairwise interactions ∇_μ T^μν_c = +Qbc^ν + QcΛ^ν, ∇_μ T^μν_b = -Qbc^ν + QbΛ^ν, ∇_μ T^μν_Λ = -QcΛ^ν - QbΛ^ν,
# given by symbolic expressions for the energy transfers Qij, their perturbations δQij and momentum transfers fQij
# Background densities are integrated backwards from their values today, instead of integrating forwards and using the shooting method.
# The tested αcΛ ranges avoid physical instabilities (for QcΛ = αcΛ*H*F: doom factor αcΛ*∂F/∂ρΛ/(3(1+w)) ≳ 1, ρc, ρΛ < 0 in the past, or ρc → 0 at early times).

using SymBoltz, Test

@independent_variables χ # lookback conformal time χ = τ₀ - τ (χ = 0 today, χ > 0 in the past)
Dχ = Differential(χ) # χ-derivative operator
Dτ(X) = -Dχ(X) # τ-derivative operator

function interacting_model(; Qbc = 0, QbΛ = 0, QcΛ = 0, δQbc = 0, δQbΛ = 0, δQcΛ = 0, fQbc = 0, fQbΛ = 0, fQcΛ = 0, name = :QΛCDM)
# Constants, some functions and atomic energy levels defined in internal files
@unpack kB, ħ, c, GN, H100, eV, me, mH, mHe, σT, aR, δkron, smoothifelse, λH2s1s, EH2s1s, EH∞2s, EHe2s1s, λHe2p1s, fHe2p1s, EHe2p2s, EHe∞2s, EHe⁺∞1s, EHet∞2s, λHet2p1s, fHet2p1s, EHet2s1s, EHet2p2s = SymBoltz
lγmax = 10
lνmax = 10
lhmax = 10
reg(x; ϵ = 1e-9) = √(x^2 + ϵ^2)
ΛH = 8.2245809
ΛHe = 51.3
A2ps = 1.798287e9
A2pt = 177.58e0
αHfit(T; F=1.125, a=4.309, b=-0.6166, c=0.6703, d=0.5300, T₀=1e4) = F * 1e-19 * a * (T/T₀)^b / (1 + c * (T/T₀)^d)
αHefit(T; q=NaN, p=NaN, T1=10^5.114, T2=3.0) = q / (√(T/T2) * (1+√(T/T2))^(1-p) * (1+√(T/T1))^(1+p))
KHfitfactorfunc(a, A, z, w) = A*exp(-((log(a)+z)/w)^2)
γHe(; A=NaN, σ=NaN, f=NaN) = 3*A*fHe*reg(1-XHe⁺)*c^2 / (8π*σ*√(2π/(β*mHe*c^2))*reg(1-XH⁺)*f^3)

# Massive neutrino distribution function and quadrature momenta
nx = 4 # number of momenta
f₀(x) = 1 / (exp(x) + 1)
dlnf₀_dlnx(x) = -x / (1 + exp(-x))
x, W = SymBoltz.momentum_quadrature(f₀, nx)
x² = x .^ 2
∫dx_x²_f₀(f) = sum(collect(f .* W))

# 1) Parameters (add your own)
pars = @parameters begin
    k, # wavenumber
    h, H0SI, # Hubble parameter in SI units (most equations have units where H0=1 and do not need these)
    Ωc0, # cold dark matter
    Ωb0, YHe, fHe, # baryons and recombination
    Tγ0, Ωγ0, # photons
    Ων0, Tν0, Nν, # massless neutrinos
    mh, mh_eV, Nh, Th0, Ωh0, yh0, Iρh0, # massive neutrinos
    ΩΛ0, w0, wa, cΛs2, # dark energy (cosmological constant or w0wa)
    zre1, Δzre1, nre1, # 1st reionization
    zre2, Δzre2, nre2, # 2nd reionization
    C, # integration constant in initial conditions
    As, ns, # primordial power spectrum
    αbc, αcΛ, αbΛ # interactions
end

# 2) Background (χ) and perturbation (χ,k) variables (add your own)
vars = @variables begin
    a(χ), [backwards = true], z(χ), ℋ(χ), H(χ), Ψ(χ,k), Φ(χ,k), τ(χ), # metric
    ρ(χ), P(χ), δρ(χ,k), Π(χ,k), # gravity
    ρb(χ), [backwards = true], Pb(χ), wb(χ), Tb(χ), δb(χ,k), Δb(χ,k), θb(χ,k), # baryons
    κ(χ), [backwards = true], v(χ), csb2(χ), β(χ), ΔT(χ), DTb(χ), μc²(χ), Xe(χ), ne(χ), λe(χ), HSI(χ), # recombination
    XH⁺(χ), nH(χ), αH(χ), βH(χ), KH(χ), KHfitfactor(χ), CH(χ), # Hydrogen recombination
    nHe(χ), XHe⁺(χ), XHe⁺⁺(χ), αHe(χ), βHe(χ), RHe⁺(χ), τHe(χ), KHe(χ), invKHe0(χ), invKHe1(χ), invKHe2(χ), CHe(χ), DXHe⁺(χ), DXHet⁺(χ), γ2ps(χ), αHet(χ), βHet(χ), τHet(χ), pHet(χ), CHet(χ), CHetnum(χ), γ2pt(χ), # Helium recombination
    Xre1(χ), Xre2(χ), # reionization
    ργ(χ), Pγ(χ), wγ(χ), Tγ(χ), Fγ0(χ,k), Fγ(χ,k)[1:lγmax], Gγ0(χ,k), Gγ(χ,k)[1:lγmax], δγ(χ,k), θγ(χ,k), σγ(χ,k), Πγ(χ,k), # photons
    ρc(χ), [backwards = true], Pc(χ), wc(χ), δc(χ,k), Δc(χ,k), θc(χ,k), # cold dark matter
    ρν(χ), Pν(χ), wν(χ), Tν(χ), Fν0(χ,k), Fν(χ,k)[1:lνmax], δν(χ,k), θν(χ,k), σν(χ,k), # massless neutrinos
    ρh(χ), Ph(χ), wh(χ), Ωh(χ), Th(χ), yh(χ), csh2(χ,k), δh(χ,k), Δh(χ,k), σh(χ,k), uh(χ,k), θh(χ,k), Eh(χ)[1:nx], ψh0(χ,k)[1:nx], ψh(χ,k)[1:nx,1:lhmax], Iρh(χ), IPh(χ), Iδρh(χ,k), # massive neutrinos
    ρΛ(χ), [backwards = true], PΛ(χ), wΛ(χ), cΛa2(χ), δΛ(χ,k), θΛ(χ,k), ΔΛ(χ,k), # dark energy (cosmological constant or w0wa)
    Qb(χ), Qc(χ), QΛ(χ), δQb(χ,k), δQc(χ,k), δQΛ(χ,k), fQb(χ,k), fQc(χ,k), fQΛ(χ,k), θ(χ,k), # interactions
    fν(χ), # misc
    ρm(χ,k), Δm(χ,k), # matter source functions
    ST_SW(χ,k), ST_ISW(χ,k), ST_Doppler(χ,k), ST_polarization(χ,k), ST(χ,k), SE(χ,k), Sψ(χ,k) # CMB source functions
end

# 3) Equations for time evolution (modify or add your own)
eqs = [
    # metric equations
    Dτ(a) ~ a * ℋ
    z ~ 1/a - 1
    ℋ ~ a * H
    Dτ(τ) ~ 1 # conformal time τ = τ₀ - χ since the big bang

    # gravity equations
    H ~ √(8π/3 * ρ) # 1st Friedmann equation
    Dτ(Φ) ~ -4π/3*a^2/ℋ*δρ - k^2/(3ℋ)*Φ - ℋ*Ψ
    k^2 * (Φ - Ψ) ~ 12π * a^2 * Π
    ρ ~ ρc + ρb + ργ + ρν + ρh + ρΛ
    P ~ Pγ + Pν + Ph + PΛ
    δρ ~ δc*ρc + δb*ρb + δγ*ργ + δν*ρν + δh*ρh + δΛ*ρΛ
    Π ~ (1+wγ)*ργ*σγ + (1+wν)*ρν*σν + (1+wh)*ρh*σh

    # baryon recombination
    β ~ 1 / (kB*Tb)
    λe ~ 2π*ħ / √(2π*me/β)
    HSI ~ H0SI * H
    Dτ(κ) ~ -a/H0SI * ne * σT * c
    v ~ expand_derivatives(Dτ(exp(-κ)))
    csb2 ~ kB/μc² * (Tb - DTb/(3ℋ))
    μc² ~ mH*c^2 / (1 + (mH/mHe-1)*YHe + Xe*(1-YHe))
    DTb ~ -2Tb*ℋ - a/h * 8/3*σT*aR/H100*Tγ^4 / (me*c) * Xe / (1+fHe+Xe) * ΔT
    Dτ(ΔT) ~ DTb + ℋ*Tγ # Tγ ∝ 1/a
    Tb ~ ΔT + Tγ
    nH ~ (1-YHe) * ρb*H0SI^2/GN / mH
    nHe ~ fHe * nH
    ne ~ Xe * nH
    Xe ~ XH⁺ + fHe*XHe⁺ + XHe⁺⁺ + Xre1 + Xre2

    # baryon H⁺ + e⁻ recombination
    αH ~ αHfit(Tb)
    βH ~ αH / λe^3 * exp(-β*EH∞2s)
    KHfitfactor ~ 1 + KHfitfactorfunc(a, -0.14, 7.28, 0.18) + KHfitfactorfunc(a, 0.079, 6.73, 0.33)
    KH ~ KHfitfactor/8π * λH2s1s^3 / HSI
    CH ~ smoothifelse(XH⁺ - 0.99, (1 + KH*ΛH*nH*(1-XH⁺)) / (1 + KH*(ΛH+βH)*nH*(1-XH⁺)), 1; k = 1e3)
    Dτ(XH⁺) ~ -a/H0SI * CH * (αH*XH⁺*ne - βH*(1-XH⁺)*exp(-β*EH2s1s))

    # baryon He⁺ + e⁻ singlet recombination
    αHe ~ αHefit(Tb; q=10^(-16.744), p=0.711)
    βHe ~ 4 * αHe / λe^3 * exp(-β*EHe∞2s)
    KHe ~ 1 / (invKHe0 + invKHe1 + invKHe2)
    invKHe0 ~ 8π*HSI / λHe2p1s^3
    τHe ~ 3*A2ps*nHe*reg(1-XHe⁺) / invKHe0
    invKHe1 ~ -exp(-τHe) * invKHe0
    γ2ps ~ γHe(A = A2ps, σ = 1.436289e-22, f = fHe2p1s)
    invKHe2 ~ A2ps/(1+0.36*γ2ps^0.86)*3*nHe*(1-XHe⁺)
    CHe ~ smoothifelse(XHe⁺ - 0.99, (exp(-β*EHe2p2s) + KHe*ΛHe*nHe*(1-XHe⁺)) / (exp(-β*EHe2p2s) + KHe*(ΛHe+βHe)*nHe*(1-XHe⁺)), 1; k = 1e3)
    DXHe⁺ ~ -a/H0SI * CHe * (αHe*XHe⁺*ne - βHe*(1-XHe⁺)*exp(-β*EHe2s1s))

    # baryon He⁺ + e⁻ triplet recombination
    αHet ~ αHefit(Tb; q=10^(-16.306), p=0.761)
    βHet ~ 4/3 * αHet / λe^3 * exp(-β*EHet∞2s)
    τHet ~ 3*A2pt*nHe*reg(1-XHe⁺) * λHet2p1s^3/(8π*HSI)
    pHet ~ (1 - exp(-τHet)) / τHet
    γ2pt ~ γHe(A = A2pt, σ = 1.484872e-22, f = fHet2p1s)
    CHetnum ~ A2pt*(pHet+1/(1+0.66*γ2pt^0.9)/3)*exp(-β*EHet2p2s)
    CHet ~ reg(CHetnum) / (reg(CHetnum) + βHet)
    DXHet⁺ ~ -a/H0SI * CHet * (αHet*XHe⁺*ne - βHet*(1-XHe⁺)*3*exp(-β*EHet2s1s))

    # baryon He⁺ + e⁻ total recombination
    Dτ(XHe⁺) ~ DXHe⁺ + DXHet⁺

    # baryon He⁺⁺ + e⁻ recombination
    RHe⁺ ~ exp(-β*EHe⁺∞1s) / (nH * λe^3)
    XHe⁺⁺ ~ 2RHe⁺*fHe / (1+fHe+RHe⁺) / (1 + √(1 + 4RHe⁺*fHe/(1+fHe+RHe⁺)^2))

    # reionization
    Xre1 ~ smoothifelse((1+zre1)^nre1 - (1+z)^nre1, 0, 1 + fHe; k = 1/(nre1*(1+zre1)^(nre1-1)*Δzre1))
    Xre2 ~ smoothifelse((1+zre2)^nre2 - (1+z)^nre2, 0, 0 + fHe; k = 1/(nre2*(1+zre2)^(nre2-1)*Δzre2))

    # baryons
    wb ~ 0
    Pb ~ wb*ρb
    Δb ~ δb + (3ℋ*(1+wb) - a*Qb/ρb)*θb/k^2 # gauge-independent with ρb′ = -3ℋ(1+wb)ρb + aQb

    # photons
    Tγ ~ Tγ0 / a
    ργ ~ 3/8π * Ωγ0 / a^4
    wγ ~ 1/3
    Pγ ~ wγ * ργ
    Dτ(Fγ0) ~ -k*Fγ[1] + 4*Dτ(Φ)
    Dτ(Fγ[1]) ~ k/3*(Fγ0-2Fγ[2]+4Ψ) - 4/3 * Dτ(κ)/k * (θb - θγ)
    [Dτ(Fγ[l]) ~ k/(2l+1) * (l*Fγ[l-1] - (l+1)*Fγ[l+1]) + Dτ(κ) * (Fγ[l] - δkron(l,2)/10*Πγ) for l in 2:lγmax-1]...
    Dτ(Fγ[lγmax]) ~ k*Fγ[lγmax-1] - (lγmax+1) / τ * Fγ[lγmax] + Dτ(κ) * Fγ[lγmax]
    δγ ~ Fγ0
    θγ ~ 3k*Fγ[1]/4
    σγ ~ Fγ[2]/2
    Πγ ~ Fγ[2] + Gγ0 + Gγ[2]
    Dτ(Gγ0) ~ -k * Gγ[1] + Dτ(κ) * (Gγ0 - Πγ/2)
    Dτ(Gγ[1]) ~ k/3 * (1*Gγ0 - 2*Gγ[2]) + Dτ(κ) * Gγ[1]
    [Dτ(Gγ[l]) ~ k/(2l+1) * (l*Gγ[l-1] - (l+1)*Gγ[l+1]) + Dτ(κ) * (Gγ[l] - δkron(l,2)/10*Πγ) for l in 2:lγmax-1]...
    Dτ(Gγ[lγmax]) ~ k*Gγ[lγmax-1] - (lγmax+1) / τ * Gγ[lγmax] + Dτ(κ) * Gγ[lγmax]

    # cold dark matter
    wc ~ 0
    Pc ~ wc*ρc
    Δc ~ δc + (3ℋ*(1+wc) - a*Qc/ρc)*θc/k^2

    # massless neutrinos
    ρν ~ 3/8π * Ων0 / a^4
    wν ~ 1/3
    Pν ~ wν * ρν
    Tν ~ Tν0 / a
    Dτ(Fν0) ~ -k*Fν[1] + 4*Dτ(Φ)
    Dτ(Fν[1]) ~ k/3*(Fν0-2Fν[2]+4Ψ)
    [Dτ(Fν[l]) ~ k/(2l+1) * (l*Fν[l-1] - (l+1)*Fν[l+1]) for l in 2:lνmax-1]...
    Dτ(Fν[lνmax]) ~ k*Fν[lνmax-1] - (lνmax+1) / τ * Fν[lνmax]
    δν ~ Fν0
    θν ~ 3k*Fν[1]/4
    σν ~ Fν[2]/2

    # massive neutrinos
    Th ~ Th0 / a
    yh ~ yh0 * a
    Iρh ~ ∫dx_x²_f₀(Eh)
    IPh ~ ∫dx_x²_f₀(x² ./ Eh)
    ρh ~ 2Nh/(2π^2) * (kB*Th)^4/(ħ*c)^3 * Iρh / ((H0SI*c)^2/GN)
    Ph ~ 2Nh/(6π^2) * (kB*Th)^4/(ħ*c)^3 * IPh / ((H0SI*c)^2/GN)
    wh ~ Ph / ρh
    Iδρh ~ ∫dx_x²_f₀(Eh .* ψh0)
    δh ~ Iδρh / Iρh
    Δh ~ δh + 3ℋ*(1+wh)*θh/k^2
    uh ~ ∫dx_x²_f₀(x .* ψh[:,1]) / (Iρh + IPh/3)
    θh ~ k * uh
    σh ~ 2/3 * ∫dx_x²_f₀(x² ./ Eh .* ψh[:,2]) / (Iρh + IPh/3)
    csh2 ~ ∫dx_x²_f₀(x² ./ Eh .* ψh0) / Iδρh
    [Eh[i] ~ √(x[i]^2 + yh^2) for i in 1:nx]...
    [Dτ(ψh0[i]) ~ -k * x[i]/Eh[i] * ψh[i,1] - Dτ(Φ) * dlnf₀_dlnx(x[i]) for i in 1:nx]...
    [Dτ(ψh[i,1]) ~ k/3 * x[i]/Eh[i] * (ψh0[i] - 2ψh[i,2]) - k/3 * Eh[i]/x[i] * Ψ * dlnf₀_dlnx(x[i]) for i in 1:nx]...
    [Dτ(ψh[i,l]) ~ k/(2l+1) * x[i]/Eh[i] * (l*ψh[i,l-1] - (l+1) * ψh[i,l+1]) for i in 1:nx, l in 2:lhmax-1]...
    [Dτ(ψh[i,lhmax]) ~ k/(2lhmax+1) * x[i]/Eh[i] * (lhmax*ψh[i,lhmax-1] - (lhmax+1) * ((2lhmax+1) * Eh[i]/x[i] * ψh[i,lhmax] / (k*τ) - ψh[i,lhmax-1])) for i in 1:nx]...

    # dark energy (cosmological constant or w0wa)
    wΛ ~ w0 + wa*(1-a)
    PΛ ~ wΛ*ρΛ
    cΛa2 ~ wΛ + ρΛ * Dτ(wΛ) / Dτ(ρΛ) # completely general
    ΔΛ ~ δΛ + (3ℋ*(1+wΛ) - a*QΛ/ρΛ)*θΛ/k^2

    # neutrino-to-radiation fraction
    fν ~ (ρν + ρh) / (ρν + ρh + ργ)

    # matter source functions
    ρm ~ ρb + ρc + ρh
    Δm ~ (ρb*Δb + ρc*Δc + ρh*Δh) / ρm

    # CMB source functions
    ST_SW ~ v * (δγ/4 + Ψ + Πγ/16)
    ST_ISW ~ exp(-κ) * Dτ(Ψ + Φ) |> expand_derivatives
    ST_Doppler ~ Dτ(v*θb) / k^2 |> expand_derivatives
    ST_polarization ~ 3/(16k^2) * Dτ(Dτ(v*Πγ)) |> expand_derivatives
    ST ~ ST_SW + ST_ISW + ST_Doppler + ST_polarization
    SE ~ 3/16 * v*Πγ / (k*χ)^2
    Sψ ~ -(Ψ + Φ)

    # total interactions on each species from the pairwise ones
    Qb ~ -Qbc + QbΛ
    Qc ~ +Qbc + QcΛ
    QΛ ~ -QcΛ - QbΛ
    δQb ~ -δQbc + δQbΛ
    δQc ~ +δQbc + δQcΛ
    δQΛ ~ -δQcΛ - δQbΛ
    fQb ~ -fQbc + fQbΛ
    fQc ~ +fQbc + fQcΛ
    fQΛ ~ -fQcΛ - fQbΛ
    θ ~ ((ρΛ+PΛ)*θΛ + (ρh+Ph)*θh + (ρν+Pν)*θν + (ρc+Pc)*θc + (ργ+Pγ)*θγ + (ρb+Pb)*θb) / ((ρΛ+PΛ) + (ρh+Ph) + (ρν+Pν) + (ρc+Pc) + (ργ+Pγ) + (ρb+Pb)) # total velocity
    Dτ(ρb) ~ -3ℋ*(1+wb)*ρb + a*Qb
    Dτ(ρc) ~ -3ℋ*(1+wc)*ρc + a*Qc
    Dτ(ρΛ) ~ -3ℋ*(1+wΛ)*ρΛ + a*QΛ
    Dτ(δb) ~ -θb - 3ℋ*csb2*δb + (a * Qb / ρb) * (Ψ - δb + 3ℋ*csb2*θb/k^2) + (a * δQb / ρb) + 3*Dτ(Φ)
    Dτ(θb) ~ -ℋ*θb + k^2*csb2*δb + k^2*Ψ + (a * Qb / ρb) * (θ - θb*(1+csb2)) + (a * k^2 / ρb) * fQb - 4/3*Dτ(κ)*ργ/ρb*(θγ-θb)
    Dτ(δc) ~ -θc + (a * Qc / ρc) * (Ψ - δc) + (a * δQc / ρc) + 3*Dτ(Φ)
    Dτ(θc) ~ -ℋ*θc + k^2*Ψ + (a * Qc / ρc) * (θ - θc) + (a * k^2 / ρc) * fQc
    Dτ(δΛ) ~ -(1+wΛ)*θΛ - 3ℋ*(cΛs2-wΛ)*δΛ - 9*(ℋ/k)^2*(1+wΛ)*(cΛs2-cΛa2)*θΛ + (a * QΛ / ρΛ) * (Ψ - δΛ + 3ℋ*(cΛs2-cΛa2)*θΛ/k^2) + (a * δQΛ / ρΛ) + 3*(1+wΛ)*Dτ(Φ)
    Dτ(θΛ) ~ -ℋ*(1-3*cΛs2)*θΛ + cΛs2/(1+wΛ)*k^2*δΛ + k^2*Ψ  + (a * QΛ / (ρΛ*(1+wΛ))) * (θ - θΛ*(1+cΛs2)) + (a * k^2 / (ρΛ*(1+wΛ))) * fQΛ
]

# 4) Equations for initial conditions (modify or add your own)
initialization_eqs = [
    # metric/gravity
    Ψ ~ 20C / (15 + 4fν)

    # baryons
    δb ~ -3/2 * Ψ
    θb ~ 1/2 * (k^2*τ) * Ψ

    # photons
    Fγ0 ~ -2Ψ
    Fγ[1] ~ 2/3 * k*τ*Ψ
    Fγ[2] ~ -8/15 * k/Dτ(κ) * Fγ[1]
    [Fγ[l] ~ -l/(2l+1) * k/Dτ(κ) * Fγ[l-1] for l in 3:lγmax]...
    Gγ0 ~ 5/16 * Fγ[2]
    Gγ[1] ~ -1/16 * k/Dτ(κ) * Fγ[2]
    Gγ[2] ~ 1/16 * Fγ[2]
    [Gγ[l] ~ -l/(2l+1) * k/Dτ(κ) * Gγ[l-1] for l in 3:lγmax]...

    # cold dark matter
    δc ~ -3/2 * Ψ
    θc ~ 1/2 * (k^2*τ) * Ψ

    # massless neutrinos
    δν ~ -2 * Ψ
    θν ~ 1/2 * (k^2*τ) * Ψ
    σν ~ 1/15 * (k*τ)^2 * Ψ
    [Fν[l] ~ l/(2l+1) * k*τ * Fν[l-1] for l in 3:lνmax]...

    # massive neutrinos
    [ψh0[i] ~ -1/4 * (-2Ψ) * dlnf₀_dlnx(x[i]) for i in 1:nx]...
    [ψh[i,1] ~ -1/3 * Eh[i]/x[i] * (1/2*k*τ*Ψ) * dlnf₀_dlnx(x[i]) for i in 1:nx]...
    [ψh[i,2] ~ -1/2 * (1/15*(k*τ)^2*Ψ) * dlnf₀_dlnx(x[i]) for i in 1:nx]...
    [ψh[i,l] ~ 0 for i in 1:nx, l in 3:lhmax]...

    # dark energy (w0wa)
    δΛ ~ -3/2 * (1+wΛ) * Ψ # for w0wa
    θΛ ~ 1/2 * (k^2*τ) * Ψ # for w0wa
]

# 5) Default numerical values for parameters and initial conditions (modify or add your own, remove to require explicit value when creating CosmologyProblem)
initial_conditions = [
    a => 1.0 # today
    τ => 1/ℋ # initialize conformal time from radiation-dominated solution
    ρb => 3/8π*Ωb0 # today
    ρc => 3/8π*Ωc0 # today
    ρΛ => 3/8π*ΩΛ0 # today
    H0SI => H100*h
    C => 1/2
    XHe⁺ => 1.0
    XH⁺ => 1.0
    κ => 0.0 # today
    ΔT => 0.0
    zre1 => 7.6711
    Δzre1 => 0.5
    nre1 => 3/2
    zre2 => 3.5
    Δzre2 => 0.5
    nre2 => 1
    Tν0 => (4/11)^(1/3) * Tγ0
    Ων0 => Nν * 7/8 * (4/11)^(4/3) * Ωγ0
    Th0 => (4/11)^(1/3) * Tγ0
    ΩΛ0 => 1 - Ωγ0 - Ωc0 - Ωb0 - Ων0 - Ωh0
    Ωγ0 => π^2/15 * (kB*Tγ0)^4 / (ħ^3*c^5) * 8π*GN / (3*H0SI^2)
    mh => mh_eV * eV/c^2
    yh0 => mh*c^2 / (kB*Th0)
    Iρh0 => ∫dx_x²_f₀(@. √(x^2 + yh0^2))
    Ωh0 => Nh * 8π/3 * 2/(2π^2) * (kB*Th0)^4 / (ħ*c)^3 * Iρh0 / ((H0SI*c)^2/GN)
    fHe => YHe / (mHe/mH*(1-YHe))
    cΛs2 => 1
]

# Equations are written with τ-derivatives; move the sign to the right side so the left sides are χ-derivatives
eqs = [Symbolics.is_derivative(Symbolics.unwrap(-eq.lhs)) ? -eq.lhs ~ -eq.rhs : eq for eq in eqs]

# 6) Pack everything down into a symbolic system (modify the name to fit your modified model)
return complete(System(eqs, χ, vars, pars; initialization_eqs, initial_conditions, name))
end

# 0) Create non-interacting base model to access variables
M = interacting_model()
@unpack a, ℋ, ρc, ρΛ, δc, δΛ, θ, θc, θΛ, Ψ, Φ, k, αcΛ = M

# 1) Energy transfer QcΛ = αcΛ*H*ρΛ with covariant perturbations (Interaction_via_H_rho_DE)
F = ρΛ
δF = ρΛ*δΛ
M1 = interacting_model(
    QcΛ = αcΛ * ℋ * F / a,
    δQcΛ = αcΛ / a * (ℋ * (δF - F*Ψ) + F * (θ/3 - Dτ(Φ))),
    fQcΛ = αcΛ * ℋ * F / (a*k^2) * (θc - θ),
)

# 2) Energy transfer QcΛ = αcΛ*H*(ρc+ρΛ) with covariant perturbations (Interaction_via_H_rho_sumas)
F = ρc + ρΛ
δF = ρc*δc + ρΛ*δΛ
M2 = interacting_model(
    QcΛ = αcΛ * ℋ * F / a,
    δQcΛ = αcΛ / a * (ℋ * (δF - F*Ψ) + F * (θ/3 - Dτ(Φ))),
    fQcΛ = αcΛ * ℋ * F / (a*k^2) * (θc - θ),
)

# 3) Energy transfer QcΛ = αcΛ*H*3ρcρΛ/(ρc+ρΛ) with covariant perturbations (Interaction_via_no_linear_prodcuto_sobre_suma_DE)
F = 3ρc*ρΛ / (ρc + ρΛ)
δF = 3ρc*ρΛ / (ρc + ρΛ)^2 * (ρΛ*δc + ρc*δΛ)
M3 = interacting_model(
    QcΛ = αcΛ * ℋ * F / a,
    δQcΛ = αcΛ / a * (ℋ * (δF - F*Ψ) + F * (θ/3 - Dτ(Φ))),
    fQcΛ = αcΛ * ℋ * F / (a*k^2) * (θc - θ),
)

# 4) Pure momentum transfer fQcΛ ∝ θΛ - θc (alphaCDM)
M4 = interacting_model(
    fQcΛ = αcΛ * 3/8π / k^2 * (θΛ - θc)
)

p1 = Dict(
    M.h => 0.7,
    M.Ωc0 => 0.25,
    M.Ωb0 => 0.05,
    M.YHe => 0.25,
    M.Tγ0 => 2.7255,
    M.Nν => 3.046 - 1,
    M.Nh => 1,
    M.mh_eV => 0.02,
    M.As => 2e-9,
    M.ns => 0.94,
    M.w0 => -1.32,
    M.wa => 0.0,
    M.αbc => 0.0,
    M.αcΛ => 0.0,
    M.αbΛ => 0.0,
)
p2 = p1
p3 = p1
p4 = merge(p1, Dict(M.w0 => -0.98, M.Ωc0 => 0.3)) # alphaCDM needs w0 > -1 (momentum drag flips sign for phantom w0)

tspan = (100.0, 0.0) # forwards in time is decreasing χ; backwards stages integrate from χ = 0 (today) until terminating at a = 1e-8
terminate = M.a ~ 1e-8
prob1 = CosmologyProblem(M1, p1; tspan, terminate)
prob2 = CosmologyProblem(M2, p2; tspan, terminate)
prob3 = CosmologyProblem(M3, p3; tspan, terminate)
prob4 = CosmologyProblem(M4, p4; tspan, terminate)

probf1 = remake_function(prob1, M1.αcΛ)
probf2 = remake_function(prob2, M2.αcΛ)
probf3 = remake_function(prob3, M3.αcΛ)
probf4 = remake_function(prob4, M4.αcΛ)

ks = [1.0, 10.0, 100.0, 1000.0]
@test all(issuccess(solve(probf1(α), ks)) for α in [-0.3, -0.1, -0.01, 0.0, 0.01, 0.1, 0.3, 1.0])
@test all(issuccess(solve(probf2(α), ks)) for α in [0.0, 0.001, 0.01, 0.03, 0.1, 0.3, 0.7])
@test all(issuccess(solve(probf3(α), ks)) for α in [-0.1, -0.01, 0.0, 0.01, 0.1, 0.3, 1.0])
@test all(issuccess(solve(probf4(α), ks)) for α in [0.0, 1.0, 10.0, 100.0])

# Matter power spectrum
sol = solve(probf1(0.1), ks)
@test spectrum_matter(sol, ks) == spectrum_matter(sol, ks, 0.0) # today is at χ = 0

# CMB power spectrum
jl = SphericalBesselCache(25:25:2000)
ls = 25:2000
Dls = spectrum_cmb([:TT, :TE, :EE], probf1(0.1), jl, ls; normalization = :Dl)
@test all(isfinite, Dls) && all(>(0), Dls[:, 1]) && all(>(0), Dls[:, 3]) # TT, EE > 0
@test ls[argmax(Dls[:, 1])] in 200:250 # first acoustic peak at ℓ ≈ 220
