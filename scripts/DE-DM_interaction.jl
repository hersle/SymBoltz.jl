using SymBoltz
using Plots
using CLASS
using Interpolations
using LinearAlgebra: BLAS
using Unitful
using UnitfulAstro
using Plots.PlotMeasures

function model(; analytical_noninteracting_continuity = true, name = :QΛCDM)
# Constants, some functions and atomic energy levels defined in internal files
@unpack kB, ħ, c, GN, H100, eV, me, mH, mHe, σT, aR, δkron, smoothifelse, λH2s1s, EH2s1s, EH∞2s, EHe2s1s, λHe2p1s, fHe2p1s, EHe2p2s, EHe∞2s, EHe⁺∞1s, EHet∞2s, λHet2p1s, fHet2p1s, EHet2s1s, EHet2p2s = SymBoltz
lγmax = 10
lνmax = 10
lhmax = 10
ϵ = 1e-9
ΛH = 8.2245809
ΛHe = 51.3
A2ps = 1.798287e9
A2pt = 177.58e0
αHfit(T; F=1.125, a=4.309, b=-0.6166, c=0.6703, d=0.5300, T₀=1e4) = F * 1e-19 * a * (T/T₀)^b / (1 + c * (T/T₀)^d)
αHefit(T; q=NaN, p=NaN, T1=10^5.114, T2=3.0) = q / (√(T/T2) * (1+√(T/T2))^(1-p) * (1+√(T/T1))^(1+p))
KHfitfactorfunc(a, A, z, w) = A*exp(-((log(a)+z)/w)^2)
γHe(; A=NaN, σ=NaN, f=NaN) = 3*A*fHe*(1-XHe⁺+ϵ)*c^2 / (8π*σ*√(2π/(β*mHe*c^2))*(1-XH⁺+ϵ)*f^3)

# Massive neutrino distribution function and quadrature momenta
nx = 4 # number of momenta
f₀(x) = 1 / (exp(x) + 1)
dlnf₀_dlnx(x) = -x / (1 + exp(-x))
x, W = SymBoltz.momentum_quadrature(f₀, nx)
x² = x .^ 2
∫dx_x²_f₀(f) = sum(collect(f .* W))

# 1) Independent variable for time evolution
@independent_variables τ # conformal time
D = Differential(τ) # derivative operator

# 2) Parameters (add your own)
pars = @parameters begin
    k, τ0, # wavenumber and conformal time today
    h, H0SI, # Hubble parameter in SI units (most equations have units where H0=1 and do not need these)
    Ωc0, # cold dark matter
    Ωb0, YHe, fHe, κ0, # baryons and recombination
    Tγ0, Ωγ0, # photons
    Ων0, Tν0, Neff, # massless neutrinos
    mh, mh_eV, Nh, Th0, Ωh0, yh0, Iρh0, # massive neutrinos
    ΩΛ0, w0, wa, cΛs2, # dark energy (cosmological constant or w0wa)
    zre1, Δzre1, nre1, # 1st reionization
    zre2, Δzre2, nre2, # 2nd reionization
    C, # integration constant in initial conditions
    As, ns, # primordial power spectrum
    αbc, αcΛ, αbΛ # interactions
end

# 3) Background (τ) and perturbation (τ,k) variables (add your own)
vars = @variables begin
    a(τ), z(τ), ℋ(τ), H(τ), Ψ(τ,k), Φ(τ,k), χ(τ), # metric
    ρ(τ), P(τ), δρ(τ,k), Π(τ,k), # gravity
    ρb(τ), [shoot=true], Pb(τ), wb(τ), Tb(τ), δb(τ,k), Δb(τ,k), θb(τ,k), # baryons
    κ(τ), _κ(τ), v(τ), csb2(τ), β(τ), ΔT(τ), DTb(τ), μc²(τ), Xe(τ), ne(τ), λe(τ), HSI(τ), # recombination
    XH⁺(τ), nH(τ), αH(τ), βH(τ), KH(τ), KHfitfactor(τ), CH(τ), # Hydrogen recombination
    nHe(τ), XHe⁺(τ), XHe⁺⁺(τ), αHe(τ), βHe(τ), RHe⁺(τ), τHe(τ), KHe(τ), invKHe0(τ), invKHe1(τ), invKHe2(τ), CHe(τ), DXHe⁺(τ), DXHet⁺(τ), γ2ps(τ), αHet(τ), βHet(τ), τHet(τ), pHet(τ), CHet(τ), CHetnum(τ), γ2pt(τ), # Helium recombination
    Xre1(τ), Xre2(τ), # reionization
    ργ(τ), Pγ(τ), wγ(τ), Tγ(τ), Fγ0(τ,k), Fγ(τ,k)[1:lγmax], Gγ0(τ,k), Gγ(τ,k)[1:lγmax], δγ(τ,k), θγ(τ,k), σγ(τ,k), Πγ(τ,k), # photons
    ρc(τ), [shoot=true], Pc(τ), wc(τ), δc(τ,k), Δc(τ,k), θc(τ,k), # cold dark matter
    ρν(τ), Pν(τ), wν(τ), Tν(τ), Fν0(τ,k), Fν(τ,k)[1:lνmax], δν(τ,k), θν(τ,k), σν(τ,k), # massless neutrinos
    ρh(τ), Ph(τ), wh(τ), Ωh(τ), Th(τ), yh(τ), csh2(τ,k), δh(τ,k), Δh(τ,k), σh(τ,k), uh(τ,k), θh(τ,k), Eh(τ)[1:nx], ψh0(τ,k)[1:nx], ψh(τ,k)[1:nx,1:lhmax], Iρh(τ), IPh(τ), Iδρh(τ,k), # massive neutrinos
    ρΛ(τ), [shoot=true], PΛ(τ), wΛ(τ), cΛa2(τ), δΛ(τ,k), θΛ(τ,k), ΔΛ(τ,k), # dark energy (cosmological constant or w0wa)
    Qb(τ), Qc(τ), QΛ(τ), Qbc(τ), QbΛ(τ), QcΛ(τ), δQb(τ,k), δQc(τ,k), δQΛ(τ,k), δQbc(τ,k), δQbΛ(τ,k), δQcΛ(τ,k), fQb(τ,k), fQc(τ,k), fQΛ(τ,k), fQbc(τ,k), fQbΛ(τ,k), fQcΛ(τ,k), θ(τ,k), # interactions
    fν(τ), # misc
    ρm(τ,k), Δm(τ,k), # matter source functions
    ST_SW(τ,k), ST_ISW(τ,k), ST_Doppler(τ,k), ST_polarization(τ,k), ST(τ,k), SE_kχ²(τ,k), Sψ(τ,k) # CMB source functions
end

# 4) Equations for time evolution (modify or add your own)
eqs = [
    # metric equations
    z ~ 1/a - 1
        ℋ ~ D(a) / a
    H ~ ℋ / a
    χ ~ τ0 - τ

    # gravity equations
    D(a) ~ √(8π/3 * ρ) * a^2 # 1st Friedmann equation
    D(Φ) ~ -4π/3*a^2/ℋ*δρ - k^2/(3ℋ)*Φ - ℋ*Ψ
    k^2 * (Φ - Ψ) ~ 12π * a^2 * Π
    ρ ~ ρc + ρb + ργ + ρν + ρh + ρΛ
    P ~ Pγ + Pν + Ph + PΛ
    δρ ~ δc*ρc + δb*ρb + δγ*ργ + δν*ρν + δh*ρh + δΛ*ρΛ
    Π ~ (1+wγ)*ργ*σγ + (1+wν)*ρν*σν + (1+wh)*ρh*σh

    # baryon recombination
    β ~ 1 / (kB*Tb)
    λe ~ 2π*ħ / √(2π*me/β)
    HSI ~ H0SI * H
    D(_κ) ~ -a/H0SI * ne * σT * c
    κ ~ _κ - κ0
    v ~ expand_derivatives(D(exp(-κ)))
    csb2 ~ kB/μc² * (Tb - D(Tb)/3ℋ)
    μc² ~ mH*c^2 / (1 + (mH/mHe-1)*YHe + Xe*(1-YHe))
    DTb ~ -2Tb*ℋ - a/h * 8/3*σT*aR/H100*Tγ^4 / (me*c) * Xe / (1+fHe+Xe) * ΔT
    D(ΔT) ~ DTb - D(Tγ)
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
    D(XH⁺) ~ -a/H0SI * CH * (αH*XH⁺*ne - βH*(1-XH⁺)*exp(-β*EH2s1s))

    # baryon He⁺ + e⁻ singlet recombination
    αHe ~ αHefit(Tb; q=10^(-16.744), p=0.711)
    βHe ~ 4 * αHe / λe^3 * exp(-β*EHe∞2s)
    KHe ~ 1 / (invKHe0 + invKHe1 + invKHe2)
    invKHe0 ~ 8π*HSI / λHe2p1s^3
    τHe ~ 3*A2ps*nHe*(1-XHe⁺+ϵ) / invKHe0
    invKHe1 ~ -exp(-τHe) * invKHe0
    γ2ps ~ γHe(A = A2ps, σ = 1.436289e-22, f = fHe2p1s)
    invKHe2 ~ A2ps/(1+0.36*γ2ps^0.86)*3*nHe*(1-XHe⁺)
    CHe ~ smoothifelse(XHe⁺ - 0.99, (exp(-β*EHe2p2s) + KHe*ΛHe*nHe*(1-XHe⁺)) / (exp(-β*EHe2p2s) + KHe*(ΛHe+βHe)*nHe*(1-XHe⁺)), 1; k = 1e3)
    DXHe⁺ ~ -a/H0SI * CHe * (αHe*XHe⁺*ne - βHe*(1-XHe⁺)*exp(-β*EHe2s1s))

    # baryon He⁺ + e⁻ triplet recombination
    αHet ~ αHefit(Tb; q=10^(-16.306), p=0.761)
    βHet ~ 4/3 * αHet / λe^3 * exp(-β*EHet∞2s)
    τHet ~ 3*A2pt*nHe*(1-XHe⁺+ϵ) * λHet2p1s^3/(8π*HSI)
    pHet ~ (1 - exp(-τHet)) / τHet
    γ2pt ~ γHe(A = A2pt, σ = 1.484872e-22, f = fHet2p1s)
    CHetnum ~ A2pt*(pHet+1/(1+0.66*γ2pt^0.9)/3)*exp(-β*EHet2p2s)
    CHet ~ (ϵ + CHetnum) / (ϵ + CHetnum + βHet)
    DXHet⁺ ~ -a/H0SI * CHet * (αHet*XHe⁺*ne - βHet*(1-XHe⁺)*3*exp(-β*EHet2s1s))

    # baryon He⁺ + e⁻ total recombination
    D(XHe⁺) ~ DXHe⁺ + DXHet⁺

    # baryon He⁺⁺ + e⁻ recombination
    RHe⁺ ~ exp(-β*EHe⁺∞1s) / (nH * λe^3)
    XHe⁺⁺ ~ 2RHe⁺*fHe / (1+fHe+RHe⁺) / (1 + √(1 + 4RHe⁺*fHe/(1+fHe+RHe⁺)^2))

    # reionization
    Xre1 ~ smoothifelse((1+zre1)^nre1 - (1+z)^nre1, 0, 1 + fHe; k = 1/(nre1*(1+zre1)^(nre1-1)*Δzre1))
    Xre2 ~ smoothifelse((1+zre2)^nre2 - (1+z)^nre2, 0, 0 + fHe; k = 1/(nre2*(1+zre2)^(nre2-1)*Δzre2))

    # baryons
    wb ~ 0
    Pb ~ wb*ρb
    Δb ~ δb + 3ℋ*θb/k^2

    # photons
    Tγ ~ Tγ0 / a
    ργ ~ 3/8π * Ωγ0 / a^4
    wγ ~ 1/3
    Pγ ~ wγ * ργ
    D(Fγ0) ~ -k*Fγ[1] + 4*D(Φ)
    D(Fγ[1]) ~ k/3*(Fγ0-2Fγ[2]+4Ψ) - 4/3 * D(κ)/k * (θb - θγ)
    [D(Fγ[l]) ~ k/(2l+1) * (l*Fγ[l-1] - (l+1)*Fγ[l+1]) + D(κ) * (Fγ[l] - δkron(l,2)/10*Πγ) for l in 2:lγmax-1]...
    D(Fγ[lγmax]) ~ k*Fγ[lγmax-1] - (lγmax+1) / τ * Fγ[lγmax] + D(κ) * Fγ[lγmax]
    δγ ~ Fγ0
    θγ ~ 3k*Fγ[1]/4
    σγ ~ Fγ[2]/2
    Πγ ~ Fγ[2] + Gγ0 + Gγ[2]
    D(Gγ0) ~ k * (-Gγ[1]) + D(κ) * (Gγ0 - Πγ/2)
    D(Gγ[1]) ~ k/(2*1+1) * (1*Gγ0 - 2*Gγ[2]) + D(κ) * Gγ[1]
    [D(Gγ[l]) ~ k/(2l+1) * (l*Gγ[l-1] - (l+1)*Gγ[l+1]) + D(κ) * (Gγ[l] - δkron(l,2)/10*Πγ) for l in 2:lγmax-1]...
    D(Gγ[lγmax]) ~ k*Gγ[lγmax-1] - (lγmax+1) / τ * Gγ[lγmax] + D(κ) * Gγ[lγmax]

    # cold dark matter
    wc ~ 0
    Pc ~ wc*ρc
    Δc ~ δc + 3ℋ*θc/k^2

    # massless neutrinos
    ρν ~ 3/8π * Ων0 / a^4
    wν ~ 1/3
    Pν ~ wν * ρν
    Tν ~ Tν0 / a
    D(Fν0) ~ -k*Fν[1] + 4*D(Φ)
    D(Fν[1]) ~ k/3*(Fν0-2Fν[2]+4Ψ)
    [D(Fν[l]) ~ k/(2l+1) * (l*Fν[l-1] - (l+1)*Fν[l+1]) for l in 2:lνmax-1]...
    D(Fν[lνmax]) ~ k*Fν[lνmax-1] - (lνmax+1) / τ * Fν[lνmax]
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
    [D(ψh0[i]) ~ -k * x[i]/Eh[i] * ψh[i,1] - D(Φ) * dlnf₀_dlnx(x[i]) for i in 1:nx]...
    [D(ψh[i,1]) ~ k/3 * x[i]/Eh[i] * (ψh0[i] - 2ψh[i,2]) - k/3 * Eh[i]/x[i] * Ψ * dlnf₀_dlnx(x[i]) for i in 1:nx]...
    [D(ψh[i,l]) ~ k/(2l+1) * x[i]/Eh[i] * (l*ψh[i,l-1] - (l+1) * ψh[i,l+1]) for i in 1:nx, l in 2:lhmax-1]...
    [D(ψh[i,lhmax]) ~ k/(2lhmax+1) * x[i]/Eh[i] * (lhmax*ψh[i,lhmax-1] - (lhmax+1) * ((2lhmax+1) * Eh[i]/x[i] * ψh[i,lhmax] / (k*τ) - ψh[i,lhmax-1])) for i in 1:nx]...

    # dark energy (cosmological constant or w0wa)
    wΛ ~ w0 + wa*(1-a)
    PΛ ~ wΛ*ρΛ
    cΛa2 ~ wΛ + ρΛ * D(wΛ) / D(ρΛ) # completely general
    ΔΛ ~ δΛ + 3ℋ*(1+wΛ)*θΛ/k^2

    # neutrino-to-radiation fraction
    fν ~ (ρν + ρh) / (ρν + ρh + ργ)

    # matter source functions
    ρm ~ ρb + ρc + ρh
    Δm ~ (ρb*Δb + ρc*Δc + ρh*Δh) / ρm

    # CMB source functions
    ST_SW ~ v * (δγ/4 + Ψ + Πγ/16)
    ST_ISW ~ exp(-κ) * D(Ψ + Φ) |> expand_derivatives
    ST_Doppler ~ D(v*θb) / k^2 |> expand_derivatives
    ST_polarization ~ 3/(16k^2) * D(D(v*Πγ)) |> expand_derivatives
    ST ~ ST_SW + ST_ISW + ST_Doppler + ST_polarization
    SE_kχ² ~ 3/16 * v*Πγ
    Sψ ~ 0 # ifelse(τ ≥ τrec, -(g.Ψ+g.Φ) * (τ-τrec)/(τ0-τrec)/(τ0-τ), 0) # TODO

    # Modified equations below:
    # Interactions: c-b-Λ triangle (the following is extremely model dependent and extremely frame dependent).
    # Pert equations for components now are general so only change here the coupling
    # Slight different symbol convection, inspired by DM being the only comoving observers that one could define in Cosmology (and in syncronous gauge):
    #                                -B_DM (bc) interaction:  ∇^μ T^μν_c = + Qbc^μ   &&   ∇^μ T^μν_b = - Qbc^μ
    #                                -b_DE (bΛ) interaction:  ∇^μ T^μν_b = + QbΛ^μ   &&   ∇^μ T^μν_Λ = - QbΛ^μ    
    #                                -c_DE (cΛ) interaction:  ∇^μ T^μν_c = + QcΛ^μ   &&   ∇^μ T^μν_Λ = - QcΛ^μ    
    #Background
    Qbc ~ 0 # b-c interaction
    QbΛ ~ 0 # b-Λ interaction
    QcΛ ~ αcΛ * ℋ * ρΛ / a  # c-Λ interaction
    Qb  ~ -Qbc + QbΛ # total Q interaction on b
    Qc  ~ +Qbc + QcΛ # total Q interaction on c
    QΛ  ~ -QcΛ - QbΛ # total Q interaction on Λ
    #Perturbation: δQ
    δQbc ~ 0
    δQbΛ ~ 0
    δQcΛ ~ αcΛ * ℋ * ρΛ / a * ( δΛ - Ψ ) + αcΛ * ρΛ / a * ( (θ/3) - D(Φ) )
    δQb  ~ -δQbc + δQbΛ # total δQ interaction on bp
    δQc  ~ +δQbc + δQcΛ # total δQ interaction on c
    δQΛ  ~ -δQcΛ - δQbΛ # total δQ interaction on Λ
    #Perturbation: fQ
    fQbc ~ 0 
    fQbΛ ~ 0 
    fQcΛ ~  ( (αcΛ * ℋ * ρΛ) / (a*k^2) ) * ( θΛ - θ ) 
    fQb  ~ -fQbc + fQbΛ # total δQ interaction on b
    fQc  ~ +fQbc + fQcΛ # total δQ interaction on c
    fQΛ  ~ -fQcΛ - fQbΛ # total δQ interaction on Λ
    #theta_frame
    θ ~ ((ρΛ+PΛ)*θΛ + (ρh+Ph)*θh + (ρν+Pν)*θν + (ρc+Pc)*θc + (ργ+Pγ)*θγ + (ρb+Pb)*θb) / ((ρΛ+PΛ) + (ρh+Ph) + (ρν+Pν) + (ρc+Pc) + (ργ+Pγ) + (ρb+Pb)) # general variable regardless of the interaction, it is a "averaged" velocity: θframe= SUM[(rho_i+P_i)*theta_i]/SUM[rho_i+P_i]
    D(ρb) ~ -3ℋ *(1+wb)*ρb + a*Qb
    D(ρc) ~ -3ℋ *(1+wc)*ρc #+ a*Qc
    D(ρΛ) ~ -3ℋ *(1+wΛ)*ρΛ #+ a*QΛ
    D(δb) ~ -θb - 3ℋ*csb2*δb + 3*D(Φ) + (a * Qb / ρb) * (Ψ - δb + 3ℋ*csb2*θb/k^2) + (a * δQb / ρb) 
    D(θb) ~ -ℋ*θb + k^2*csb2*δb + k^2*Ψ - 4/3*D(κ)*ργ/ρb*(θγ-θb) + (a * Qb / ρb) * (θ - θb*(1+csb2)) + (a * k^2 / ρb) * fQb
    D(δc) ~ -θc + 3*D(Φ) #+ (a * Qc / ρc) * (Ψ - δc) + (a * δQc / ρc)
    D(θc) ~ -ℋ*θc + k^2*Ψ #+ (a * Qc / ρc) * (θ - θc) + (a * k^2 / ρc) * fQc
    D(δΛ) ~ -(1+wΛ)*(θΛ-3*D(Φ)) - 3ℋ*(cΛs2-wΛ)*δΛ - 9*(ℋ/k)^2*(1+wΛ)*(cΛs2-cΛa2)*θΛ #+ (a * QΛ / (ρΛ*(1+wΛ))) * (Ψ - δΛ + 3ℋ*(cΛs2-cΛa2)*θΛ/k^2) + (a * δQΛ / (ρΛ*(1+wΛ)))
    D(θΛ) ~ -ℋ*(1-3*cΛs2)*θΛ + cΛs2/(1+wΛ)*k^2*δΛ + k^2*Ψ  #+ (a * QΛ / (ρΛ*(1+wΛ)  ) ) * (θ - θΛ*(1+cΛs2)) + (a * k^2 / (ρΛ*(1+wΛ) ) ) * fQΛ
]

# 5) Equations for initial conditions (modify or add your own)
initialization_eqs = [
    # metric/gravity
    Ψ ~ 20C / (15 + 4fν)
    D(a) ~ a / τ

    # baryons
    δb ~ -3/2 * Ψ
    θb ~ 1/2 * (k^2*τ) * Ψ

    # photons
    Fγ0 ~ -2Ψ
    Fγ[1] ~ 2/3 * k*τ*Ψ
    Fγ[2] ~ -8/15 * k/D(κ) * Fγ[1]
    [Fγ[l] ~ -l/(2l+1) * k/D(κ) * Fγ[l-1] for l in 3:lγmax]...
    Gγ0 ~ 5/16 * Fγ[2]
    Gγ[1] ~ -1/16 * k/D(κ) * Fγ[2]
    Gγ[2] ~ 1/16 * Fγ[2]
    [Gγ[l] ~ -l/(2l+1) * k/D(κ) * Gγ[l-1] for l in 3:lγmax]...

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

# 6) Initial guess for variables solved for in initial conditions and shooting method (modify or add your own)
guesses = [
    a => τ # a(τini) is solved for in a nonlinear system constrained to ℋ(aini) ~ 1/τini (see initialization_eqs)
    ρb => τ^(-3)
    ρc => τ^(-3)
    ρΛ => τ^(-3(1+w0+wa)) * exp(-3wa*(1-τ))
]

# 7) Shooting constraints (evaluated today)
constraints = [
    ρb ~ 3/8π*Ωb0
    ρc ~ 3/8π*Ωc0
    ρΛ ~ 3/8π*ΩΛ0
]

# 8) Default numerical values for parameters and initial conditions (modify or add your own, remove to require explicit value when creating CosmologyProblem)
initial_conditions = [
    H0SI => H100*h
    τ0 => NaN
    C => 1/2
    XHe⁺ => 1.0
    XH⁺ => 1.0
    _κ => 0.0
    κ0 => NaN
    ΔT => 0.0
    zre1 => 7.6711
    Δzre1 => 0.5
    nre1 => 3/2
    zre2 => 3.5
    Δzre2 => 0.5
    nre2 => 1
    Tν0 => (4/11)^(1/3) * Tγ0
    Ων0 => Neff * 7/8 * (4/11)^(4/3) * Ωγ0
    Nh => 3
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

# Optional: use analytical solutions for noninteracting continuity equations
if analytical_noninteracting_continuity
isequal(expandeq(eqs, Qb), 0) && push!(eqs, ρb ~ 3/8π * Ωb0 / a^3)
isequal(expandeq(eqs, Qc), 0) && push!(eqs, ρc ~ 3/8π * Ωc0 / a^3)
isequal(expandeq(eqs, QΛ), 0) && push!(eqs, ρΛ ~ 3/8π * ΩΛ0 * abs(a)^(-3*(1+w0+wa)) * exp(-3wa*(1-a)))
anal = intersect(Set([ρΛ, ρb, ρc]), Set(eq.lhs for eq in eqs)) # which energy densities do we have the analytical solution for?
Danal = Set(D.(anal))
filter!(eq -> !(eq.lhs in Danal), eqs) # remove ODEs where we have the analytical solution
filter!(eq -> !(eq.lhs in anal), constraints) # remove shooting constraint
filter!(guess -> !(guess[1] in anal), guesses) # remove shooting guess
end

# 9) Pack everything down into a symbolic system (modify the name to fit your modified model)
return complete(System(eqs, τ, vars, pars; initialization_eqs, initial_conditions, guesses, constraints, name))
end

# ============================================================
# MAIN COMPARISON WITH CLASS: alphaCDM scan
# ============================================================

# This part replaces the original quick test at the end of the file.
# It keeps the SymBoltz model above unchanged and adds:
#   1) background comparison against CLASS
#   2) matter power spectrum P(k) comparison against CLASS
#   3) CMB TT comparison against CLASS (disabled by default here)
#
# The model is prepared for an alphaCDM scan over αcΛ values.
# The functions Q, δQ and fQ are intentionally left untouched in the model above.
# Edit them there if you want to activate a specific interacting model.

# -----------------------------
# CLASS executable path
# -----------------------------
const CLASS_EXEC = "/mnt/wwn-0x50014ee2c101e66c-part1/Zona_de_Trabajo/Proyecto_SymBoltz/class/CLASS_DE-DM_interacting_marcel-2-5/class"

# If your CLASS version accepts an interaction parameter for this model, set it here.
# If you are using standard CLASS/wCDM, keep it as nothing.
const CLASS_ALPHA_PARAMETER = "delta_de_copling"

# -----------------------------
# Output paths
# -----------------------------
const OUT_BACKGROUND = joinpath(pwd(), "alphaCDM_background_CLASS_vs_SymBoltz.png")
const OUT_PK         = joinpath(pwd(), "alphaCDM_Pk_CLASS_vs_SymBoltz.png")
const OUT_TT         = joinpath(pwd(), "alphaCDM_TT_CLASS_vs_SymBoltz.png")

# -----------------------------
# User switches
# -----------------------------
run_background = true
run_pk         = true
run_tt         = false

# -----------------------------
# Alpha scan
# -----------------------------
alphas =  [0.0]#[0]

# -----------------------------
# P(k) grid
# -----------------------------
kmin = 1e-4
kmax = 1.0
nk   = 100

# -----------------------------
# CMB grid
# -----------------------------
lmin = 25
lmax = 1000
dl_cache = 25
ls_tt = lmin:lmax
jl_tt = SphericalBesselCache(lmin:dl_cache:lmax)

# -----------------------------
# Perturbation wavenumber used for the time-evolution plots
# -----------------------------
k_evolution = 1e3

# -----------------------------
# Helper: build labels for the SymBoltz parameter dictionaries
# -----------------------------
function interaction_label(M, p)
    return "α=$(p[M.αcΛ])"
end

# -----------------------------
# Helper: robust access to CLASS tables
# -----------------------------
function class_table(class, candidates)
    for name in candidates
        try
            return class[name]
        catch
        end
        try
            return class[Symbol(name)]
        catch
        end
    end
    error("Could not find any CLASS table among: " * join(string.(candidates), ", "))
end

# -----------------------------
# Helper: robust access to CLASS columns
# -----------------------------
function table_column(tbl, candidates)
    nms = names(tbl)
    snms = string.(nms)

    for candidate in candidates
        j = findfirst(==(candidate), snms)
        if j !== nothing
            return Float64.(tbl[!, nms[j]])
        end
    end

    # fallback: substring match, useful because CLASS column names can vary
    for candidate in candidates
        j = findfirst(x -> occursin(candidate, x), snms)
        if j !== nothing
            return Float64.(tbl[!, nms[j]])
        end
    end

    error("Could not find any column among: " * join(candidates, ", ") *
          ". Available columns are: " * join(snms, ", "))
end

# -----------------------------
# Helper: read background data from CLASS solution
# -----------------------------
function find_col(tab, candidates)
    names_tab = names(tab)

    for c in candidates
        if c in names_tab
            return c
        end
    end

    error("None of these columns were found: $(candidates). Available columns are: $(names_tab)")
end

function read_class_background(class, h)
    bg = class[:background]

    z_col  = find_col(bg, ["z"])
    H_col  = find_col(bg, ["H [1/Mpc]", "H"])
    rhoc_col = find_col(bg, ["(.)rho_cdm", "rho_cdm", "rho_cdm [Mpc^-2]"])
    rhob_col = find_col(bg, ["(.)rho_b", "rho_b", "rho_b [Mpc^-2]"])
    rhoX_col = find_col(bg, ["(.)rho_fld", "rho_fld", "rho_fld [Mpc^-2]", "(.)rho_lambda", "rho_lambda"])

    z_cls = Float64.(bg[!, z_col])
    a_cls = 1.0 ./ (1.0 .+ z_cls)

    H_cls_raw = Float64.(bg[!, H_col])
    H0_cls_raw = H_cls_raw[argmin(abs.(z_cls))]
    HoverH0_cls = H_cls_raw ./ H0_cls_raw

    # CLASS and SymBoltz use different normalizations for the background densities.
    #
    # CLASS background columns are in units proportional to H0^2.
    # SymBoltz uses the internal normalization:
    #
    #     ρ_i = 3/(8π) * Ω_i(a)
    #
    # Therefore, keep the CLASS time dependence exactly as it is read from CLASS,
    # but convert its normalization to the SymBoltz convention.
    ρc_cls_raw = Float64.(bg[!, rhoc_col])
    ρb_cls_raw = Float64.(bg[!, rhob_col])
    ρΛ_cls_raw = Float64.(bg[!, rhoX_col])

    ρc_cls = 3/(8π) .* ρc_cls_raw ./ H0_cls_raw^2
    ρb_cls = 3/(8π) .* ρb_cls_raw ./ H0_cls_raw^2
    ρΛ_cls = 3/(8π) .* ρΛ_cls_raw ./ H0_cls_raw^2

    return a_cls, ρc_cls, ρb_cls, ρΛ_cls, HoverH0_cls
end

# -----------------------------
# Helper: CLASS setup matching the SymBoltz parameter dictionary
# -----------------------------
function make_class_params(M, p; alpha_value = 0.0)
    params = Dict(
        "h" => p[M.h],
        "Omega_cdm" => p[M.Ωc0],
        "Omega_b" => p[M.Ωb0],
        "YHe" => p[M.YHe],
        "T_cmb" => p[M.Tγ0],

        # Neutrino setup matching the SymBoltz script:
        # one ncdm species with degeneracy 3 and mass mh_eV.
        "N_ur" => p[M.Neff],
        "N_ncdm" => 1,
        "deg_ncdm" => 3,
        "m_ncdm" => p[M.mh_eV],
        "T_ncdm" => (4/11)^(1/3),

        # w0-wa fluid dark energy
        "Omega_Lambda" => 0.0,
        "Omega_scf" => 0.0,
        "fluid_equation_of_state" => "CLP",
        "w0_fld" => p[M.w0],
        "wa_fld" => p[M.wa],
        "cs2_fld" => 1.0,
        "use_ppf" => "no",

        # Primordial spectrum
        "A_s" => p[M.As],
        "n_s" => p[M.ns],

        # Output settings
        "output" => vcat(run_pk ? ["mPk"] : String[], run_tt ? ["tCl"] : String[]),
        "P_k_max_h/Mpc" => kmax,
        "z_pk" => 0,
        "delta_dm_copling" => 0.0,
        
        

        # Ask CLASS to write background so CLASS.jl can expose it as a table.
        "write background" => "yes",

        # Use Newtonian gauge to match the SymBoltz equations.
        "gauge" => "newtonian",
    )

    if CLASS_ALPHA_PARAMETER !== nothing
        params[CLASS_ALPHA_PARAMETER] = alpha_value / 3 #En class esta sin el 1/3, cogemos y dividimos o lo absorvemos por el parametro de la interacion 
    end

    return params
end

# -----------------------------
# Helper: ratio for background lower panels
# -----------------------------
function background_ratio(a_sym, y_sym, x_cls, y_cls)
    x_sym_all = log10.(a_sym)
    mask_sym = isfinite.(x_sym_all) .& isfinite.(y_sym) .& (a_sym .> 0)

    x_sym_all = x_sym_all[mask_sym]
    y_sym_all = y_sym[mask_sym]

    x_lo = max(minimum(x_cls), minimum(x_sym_all))
    x_hi = min(maximum(x_cls), maximum(x_sym_all))

    mask_overlap = (x_sym_all .>= x_lo) .& (x_sym_all .<= x_hi)
    x_sym = x_sym_all[mask_overlap]
    y_sym = y_sym_all[mask_overlap]

    itp = LinearInterpolation(x_cls, y_cls, extrapolation_bc=Line())
    y_cls_on_sym = itp.(x_sym)

    mask_ok = isfinite.(y_cls_on_sym) .& (abs.(y_cls_on_sym) .> 0)
    x_sym = x_sym[mask_ok]
    y_sym = y_sym[mask_ok]
    y_cls_on_sym = y_cls_on_sym[mask_ok]

    relerr = (y_sym .- y_cls_on_sym) ./ y_cls_on_sym
    return x_sym, relerr
end

# -----------------------------
# Helper: analytic DE/DM background expressions
# -----------------------------
function rhoDE_analytic_user(a; a0, rhoDE0, w0, wa, α)
    x = a ./ a0
    return rhoDE0 .* x.^(-3*(1+w0+α)) 
end

function rhoDM_analytic_user(a; a0, rhoDM0, rhoDE0, w0, wa, α)
    x = a ./ a0
    return x.^(-3) .* (rhoDM0 .+ ( (α)/(w0+α) ) .* rhoDE0 .* ( 1 .- x.^(-α-w0) ))
end


# -----------------------------
# Helper: analytic-vs-SymBoltz lower panel error
# -----------------------------
function analytic_log_relative_error(y_analytic, y_sym)
    mask_ok =
        isfinite.(y_analytic) .&
        isfinite.(y_sym) .&
        (abs.(y_analytic) .> 0)

    relerr = (y_sym[mask_ok] .- y_analytic[mask_ok]) ./ y_analytic[mask_ok]

    err_floor = 1e-12
    logerr = log10.(max.(abs.(relerr), err_floor))

    return mask_ok, logerr
end
# -----------------------------
# MAIN
# -----------------------------
function main()
    BLAS.set_num_threads(1)

    # -----------------------------
    # Build SymBoltz model
    # -----------------------------
    M = model(name = :alphaCDM)

    # -----------------------------
    # Base SymBoltz parameter dictionary
    # -----------------------------
    base_p = Dict(
        M.h => 0.7,
        M.Ωc0 => 0.3,
        M.Ωb0 => 0.05,
        M.YHe => 0.25,
        M.Tγ0 => 2.7255,
        M.Neff => 3.046,
        M.mh_eV => 0.02,
        M.As => 2e-9,
        M.ns => 0.94,
        M.w0 => -0.98,
        M.wa => 0.0,
        M.αbc => 0.0,
        M.αbΛ => 0.0,
        M.αcΛ => 0.0,
    )

    sym_cases = []
    for αval in alphas
        p = copy(base_p)
        p[M.αbc] = 0.0
        p[M.αbΛ] = 0.0
        p[M.αcΛ] = Float64(αval)
        push!(sym_cases, ("SymBoltz α=$(αval)", p, αval))
    end

    # -----------------------------
    # Cosmology problems
    # -----------------------------
    probs = [(label, CosmologyProblem(M, p), αval) for (label, p, αval) in sym_cases]

    # -----------------------------
    # Run CLASS for each alpha value.
    # If CLASS_ALPHA_PARAMETER = nothing, these are standard wCDM runs and are identical;
    # if you set CLASS_ALPHA_PARAMETER = "alpha_model", each CLASS run receives the same alpha.
    # -----------------------------
    class_solutions = Dict{Float64, Any}()
    for (label, p, αval) in sym_cases
        class_params = make_class_params(M, p; alpha_value = αval)
        println("Running CLASS for α=$(αval) with outputs = ", class_params["output"])
        class_solutions[αval] = solve(CLASSProblem(class_params...); exec = CLASS_EXEC)
        println("Loaded CLASS solution for α=$(αval)")
    end

    pal = palette(:auto)
    colors = pal[1:length(alphas)]
    color_for_alpha = Dict(alphas[i] => colors[i] for i in eachindex(alphas))

    # ============================================================
    # 1) Background comparison
    # ============================================================
    if run_background
        println("Computing SymBoltz background/evolution for comparison...")

        sol_data = []
        for (label, prob, αval) in probs
            println(label)
            sol = solve(prob, k_evolution; verbose = true)
            push!(sol_data, (label, sol, αval))

            println("ρc(τ₀) = ", sol[M.ρc][end])
            println("ρb(τ₀) = ", sol[M.ρb][end])
            println("ρΛ(τ₀) = ", sol[M.ρΛ][end])
            println("H(τ₀) = ", sol[M.H][end])
            println("a(τ₀) = ", sol[M.a][end])
            QcΛ_vals = Float64.(sol[M.QcΛ])
            Qc_vals  = Float64.(sol[M.Qc])
            QΛ_vals  = Float64.(sol[M.QΛ])
            println("max |QcΛ| = ", maximum(abs.(QcΛ_vals)))
            println("max |Qc|  = ", maximum(abs.(Qc_vals)))
            println("max |QΛ|  = ", maximum(abs.(QΛ_vals)))
            target_ρc = 3/(8π) * base_p[M.Ωc0]
            target_ρb = 3/(8π) * base_p[M.Ωb0]
            target_ρΛ = 3/(8π) * (1 - sol[M.Ωγ0][end] - base_p[M.Ωc0] - base_p[M.Ωb0] - sol[M.Ων0][end] - sol[M.Ωh0][end])
            println("relerr ρc0 shooting = ", (sol[M.ρc][end] - target_ρc) / target_ρc)
            println("relerr ρb0 shooting = ", (sol[M.ρb][end] - target_ρb) / target_ρb)
            println("relerr ρΛ0 shooting = ", (sol[M.ρΛ][end] - target_ρΛ) / target_ρΛ)
            
            
            
        end

        # Read and store CLASS backgrounds per alpha.
        class_bg = Dict{Float64, Any}()
        xmin_background = -8.0
        x_min_all = Inf
        x_max_all = -Inf

        for αval in alphas
            class = class_solutions[αval]
            a_cls, ρc_cls, ρb_cls, ρΛ_cls, H_cls = read_class_background(class, base_p[M.h])

            mask_cls =
                (a_cls .> 0) .&
                isfinite.(a_cls) .&
                (log10.(a_cls) .>= xmin_background)

            x_cls_all = log10.(a_cls[mask_cls])
            order_cls = sortperm(x_cls_all)
            x_cls = x_cls_all[order_cls]

            ρc_cls_plot = ρc_cls[mask_cls][order_cls]
            ρb_cls_plot = ρb_cls[mask_cls][order_cls]
            ρΛ_cls_plot = ρΛ_cls[mask_cls][order_cls]
            H_cls_plot  = H_cls[mask_cls][order_cls]

            class_bg[αval] = (x_cls, ρc_cls_plot, ρb_cls_plot, ρΛ_cls_plot, H_cls_plot)

            x_min_all = min(x_min_all, minimum(x_cls))
            x_max_all = max(x_max_all, maximum(x_cls))
        end

        # Helper for lower panels:
        # interpolate CLASS onto the SymBoltz a-grid and compute
        # log10(abs((SymBoltz - CLASS) / CLASS)).
        function background_log_relative_error(a_sym, y_sym, x_cls, y_cls)
            x_sym_all = log10.(a_sym)

            mask_sym =
                isfinite.(x_sym_all) .&
                isfinite.(y_sym) .&
                (a_sym .> 0)

            x_sym_all = x_sym_all[mask_sym]
            y_sym_all = y_sym[mask_sym]

            x_lo = max(minimum(x_cls), minimum(x_sym_all))
            x_hi = min(maximum(x_cls), maximum(x_sym_all))

            mask_overlap =
                (x_sym_all .>= x_lo) .&
                (x_sym_all .<= x_hi)

            x_sym = x_sym_all[mask_overlap]
            y_sym = y_sym_all[mask_overlap]

            itp = LinearInterpolation(x_cls, y_cls, extrapolation_bc = Line())
            y_cls_on_sym = itp.(x_sym)

            mask_ok =
                isfinite.(y_cls_on_sym) .&
                isfinite.(y_sym) .&
                (abs.(y_cls_on_sym) .> 0)

            x_sym = x_sym[mask_ok]
            y_sym = y_sym[mask_ok]
            y_cls_on_sym = y_cls_on_sym[mask_ok]

            relerr = (y_sym .- y_cls_on_sym) ./ y_cls_on_sym

            err_floor = 1e-12
            logerr = log10.(max.(abs.(relerr), err_floor))

            return x_sym, logerr
        end

        # Top row: background quantities.
        # Bottom row:
        #   solid curves  -> log10(abs((SymBoltz - CLASS) / CLASS))
        #   dotted curves -> log10(abs((Analytic - SymBoltz) / SymBoltz))
        pltρc = plot(
            xlabel = "",
            xticks = false,
            bottom_margin = 0mm,
            ylabel = "log10(ρc)",
            title = "Cold dark matter",
            legend = :outertop,
        )

        pltρb = plot(
            xlabel = "",
            xticks = false,
            bottom_margin = 0mm,
            ylabel = "log10(ρb)",
            title = "Baryons",
            legend = false,
        )

        pltρΛ = plot(
            xlabel = "",
            xticks = false,
            bottom_margin = 0mm,
            ylabel = "log10(ρΛ)",
            title = "Dark energy",
            legend = false,
        )

        pltH = plot(
            xlabel = "",
            xticks = false,
            bottom_margin = 0mm,
            ylabel = "log10(H/H0)",
            title = "Hubble rate",
            legend = false,
        )

        pltρcErr = plot(
            top_margin = 0mm,
            xlabel = "log10(a)",
            ylabel = "log10|ερc|",
            legend = false,
            xlims = (x_min_all, x_max_all),
        )

        pltρbErr = plot(
            top_margin = 0mm,
            xlabel = "log10(a)",
            ylabel = "log10|ερb|",
            legend = false,
            xlims = (x_min_all, x_max_all),
        )

        pltρΛErr = plot(
            top_margin = 0mm,
            xlabel = "log10(a)",
            ylabel = "log10|ερΛ|",
            legend = false,
            xlims = (x_min_all, x_max_all),
        )

        pltHErr = plot(
            top_margin = 0mm,
            xlabel = "log10(a)",
            ylabel = "log10|εE|",
            legend = false,
            xlims = (x_min_all, x_max_all),
        )

        # CLASS curves first: dashed, same color as the corresponding alpha.
        for αval in alphas
            col = color_for_alpha[αval]
            x_cls, ρc_cls_plot, ρb_cls_plot, ρΛ_cls_plot, H_cls_plot = class_bg[αval]

            plot!(pltρc, x_cls, log10.(abs.(ρc_cls_plot)),
                  color = col, linestyle = :dash, label = "CLASS ε=$(αval)")

            plot!(pltρb, x_cls, log10.(abs.(ρb_cls_plot)),
                  color = col, linestyle = :dash, label = "CLASS ε=$(αval)")

            plot!(pltρΛ, x_cls, log10.(abs.(ρΛ_cls_plot)),
                  color = col, linestyle = :dash, label = "CLASS ε=$(αval)")

            plot!(pltH,  x_cls, log10.(abs.(H_cls_plot)),
                  color = col, linestyle = :dash, label = "CLASS ε=$(αval)")
        end

        for (label, sol, αval) in sol_data
            col = color_for_alpha[αval]
            x_cls, ρc_cls_plot, ρb_cls_plot, ρΛ_cls_plot, H_cls_plot = class_bg[αval]

            a_sym  = Float64.(sol[M.a])
            ρc_sym = Float64.(sol[M.ρc])
            ρb_sym = Float64.(sol[M.ρb])
            ρΛ_sym = Float64.(sol[M.ρΛ])
            H_sym  = Float64.(sol[M.H])

            x_sym = log10.(a_sym)

            # ------------------------------------------------------------
            # SymBoltz curves: solid
            # ------------------------------------------------------------
            plot!(pltρc, x_sym, log10.(abs.(ρc_sym)),
                  linestyle = :solid, color = col, label = label)

            plot!(pltρb, x_sym, log10.(abs.(ρb_sym)),
                  linestyle = :solid, color = col, label = label)

            plot!(pltρΛ, x_sym, log10.(abs.(ρΛ_sym)),
                  linestyle = :solid, color = col, label = label)

            plot!(pltH,  x_sym, log10.(abs.(H_sym)),
                  linestyle = :solid, color = col, label = label)

            # ------------------------------------------------------------
            # Analytic DE/DM curves: dotted
            # Only plotted in DM and DE panels.
            # Normalized using today's SymBoltz values.
            # ------------------------------------------------------------
            a0_sym  = a_sym[end]
            ρc0_sym = ρc_sym[end]
            ρΛ0_sym = ρΛ_sym[end]

            w0_val = p_for_alpha = nothing

            # Recover the parameter dictionary corresponding to this alpha.
            # This avoids relying on base_p if later you change w0/wa per case.
            p_case = nothing
            for (label2, p2, αval2) in sym_cases
                if αval2 == αval
                    p_case = p2
                    break
                end
            end

            if p_case === nothing
                error("Could not find parameter dictionary for α=$(αval)")
            end

            w0_val = p_case[M.w0]
            wa_val = p_case[M.wa]

            ρΛ_analytic = rhoDE_analytic_user(
                a_sym;
                a0 = a0_sym,
                rhoDE0 = ρΛ0_sym,
                w0 = w0_val,
                wa = wa_val,
                α = αval/3, #En class esta sin el 1/3, cogemos y dividimos o lo absorvemos por el parametro de la interacion 
            )

            ρc_analytic = rhoDM_analytic_user(
                a_sym;
                a0 = a0_sym,
                rhoDM0 = ρc0_sym,
                rhoDE0 = ρΛ0_sym,
                w0 = w0_val,
                wa = wa_val,
                α = αval/3, #En class esta sin el 1/3, cogemos y dividimos o lo absorvemos por el parametro de la interacion 
            )

            mask_ana_c =
                isfinite.(x_sym) .&
                isfinite.(ρc_analytic) .&
                (a_sym .> 0) .&
                (ρc_analytic .!= 0)

            mask_ana_Λ =
                isfinite.(x_sym) .&
                isfinite.(ρΛ_analytic) .&
                (a_sym .> 0) .&
                (ρΛ_analytic .!= 0)

            plot!(pltρc, x_sym[mask_ana_c], log10.(abs.(ρc_analytic[mask_ana_c])),
                  linestyle = :dot, color = col, label = "Analytic ε=$(αval)")

            plot!(pltρΛ, x_sym[mask_ana_Λ], log10.(abs.(ρΛ_analytic[mask_ana_Λ])),
                  linestyle = :dot, color = col, label = "Analytic ε=$(αval)")

            # ------------------------------------------------------------
            # Existing lower panels: SymBoltz vs CLASS
            # ------------------------------------------------------------
            xerr, logerr = background_log_relative_error(a_sym, ρc_sym, x_cls, ρc_cls_plot)
            plot!(pltρcErr, xerr, logerr,
                  linestyle = :solid, color = col, label = "")

            xerr, logerr = background_log_relative_error(a_sym, ρb_sym, x_cls, ρb_cls_plot)
            plot!(pltρbErr, xerr, logerr,
                  linestyle = :solid, color = col, label = "")

            xerr, logerr = background_log_relative_error(a_sym, ρΛ_sym, x_cls, ρΛ_cls_plot)
            plot!(pltρΛErr, xerr, logerr,
                  linestyle = :solid, color = col, label = "")

            xerr, logerr = background_log_relative_error(a_sym, H_sym, x_cls, H_cls_plot)
            plot!(pltHErr, xerr, logerr,
                  linestyle = :solid, color = col, label = "")

            # ------------------------------------------------------------
            # New lower panels: Analytic vs SymBoltz
            # Only plotted in DM and DE error panels.
            # ------------------------------------------------------------
            mask_err_c, logerr_c_analytic = analytic_log_relative_error(ρc_analytic, ρc_sym)
            plot!(pltρcErr, x_sym[mask_err_c], logerr_c_analytic,
                  linestyle = :dot, color = col, label = "")

            mask_err_Λ, logerr_Λ_analytic = analytic_log_relative_error(ρΛ_analytic, ρΛ_sym)
            plot!(pltρΛErr, x_sym[mask_err_Λ], logerr_Λ_analytic,
                  linestyle = :dot, color = col, label = "")
        end

        figBg = plot(
            pltρc, pltρb, pltρΛ, pltH,
            pltρcErr, pltρbErr, pltρΛErr, pltHErr,
            layout = (2, 4),
            size = (1800, 900),
            left_margin = 10mm,
            bottom_margin = 10mm,
        )

        savefig(figBg, OUT_BACKGROUND)
        println("Saved background comparison to: ", OUT_BACKGROUND)

        pltρc = nothing
        pltρb = nothing
        pltρΛ = nothing
        pltH = nothing
        pltρcErr = nothing
        pltρbErr = nothing
        pltρΛErr = nothing
        pltHErr = nothing
        figBg = nothing
        sol_data = nothing
        class_bg = nothing
        GC.gc()
    end

    # ============================================================
    # 2) Matter power spectrum comparison
    # ============================================================
    if run_pk
        println("Computing P(k) from SymBoltz...")

        kconv = 100 * SymBoltz.km / SymBoltz.c
        ks_sym = (10 .^ range(log10(kmin), log10(kmax), length=nk)) ./ kconv
        k_sym_hmpc = Float64.(ks_sym .* kconv)

        pltPk = plot(
            xlabel="",
            xticks=false,
            ylabel="log10(P / (Mpc/h)^3)",
            legend=:bottomleft,
            bottom_margin=0mm
        )

        pltPkErr = plot(
            xlabel="log10(k / (h/Mpc))",
            ylabel="(SymBoltz - CLASS) / CLASS",
            legend=false,
            xlims=(log10(kmin), log10(kmax)),
            top_margin=0mm
        )
        hline!(pltPkErr, [0.0], linestyle=:dash, color=:black, label="")

        for (label, prob, αval) in probs
            col = color_for_alpha[αval]
            class = class_solutions[αval]

            k_cls = Float64.(class[:pk][!, "k (h/Mpc)"])
            P_cls = Float64.(class[:pk][!, "P (Mpc/h)^3"])

            mask_cls = (k_cls .>= kmin) .& (k_cls .<= kmax)
            k_cls = k_cls[mask_cls]
            P_cls = P_cls[mask_cls]

            Ps = spectrum_matter(prob, ks_sym; shootopts = (alg = SymBoltz.shootalg(), abstol = 1e-8, reltol = 1e-8 ) )
            P_sym = Float64.(Ps ./ kconv^3)

            plot!(pltPk, log10.(k_cls), log10.(P_cls),
                  label="CLASS α=$(αval)", color=col, linestyle=:dash)
            plot!(pltPk, log10.(k_sym_hmpc), log10.(P_sym),
                  label=label, color=col, linestyle=:solid)

            k_lo = max(minimum(k_cls), minimum(k_sym_hmpc), kmin)
            k_hi = min(maximum(k_cls), maximum(k_sym_hmpc), kmax)

            mask_sym = (k_sym_hmpc .>= k_lo) .& (k_sym_hmpc .<= k_hi)
            k_sym2 = k_sym_hmpc[mask_sym]
            P_sym2 = P_sym[mask_sym]

            itpPk = LinearInterpolation(log.(k_cls), log.(P_cls), extrapolation_bc=Line())
            P_cls_on_sym = exp.(itpPk.(log.(k_sym2)))

            relerr = (P_sym2 .- P_cls_on_sym) ./ P_cls_on_sym

            plot!(pltPkErr, log10.(k_sym2), relerr,
                  label="", color=col, linestyle=:solid)
        end

        figPk = plot(pltPk, pltPkErr, layout=(2,1), size=(900,900),
                     left_margin=8mm, bottom_margin=8mm)
        savefig(figPk, OUT_PK)
        println("Saved P(k) comparison to: ", OUT_PK)

        pltPk = nothing
        pltPkErr = nothing
        figPk = nothing
        GC.gc()
    end

    # ============================================================
    # 3) CMB TT comparison
    # ============================================================
    if run_tt
        println("Computing TT from SymBoltz...")

        pltTT = plot(
            xlabel="",
            xticks=false,
            ylabel="10¹² DℓTT",
            legend=:bottomleft,
            bottom_margin=0mm
        )

        pltTTErr = plot(
            xlabel="ℓ",
            ylabel="(SymBoltz - CLASS) / CLASS",
            legend=false,
            xlims=(minimum(ls_tt), maximum(ls_tt)),
            top_margin=0mm
        )
        hline!(pltTTErr, [0.0], linestyle=:dash, color=:black, label="")

        for (label, prob, αval) in probs
            col = color_for_alpha[αval]
            class = class_solutions[αval]

            ell_cls = Float64.(class[:cl][!, "l"])
            TT_cls = Float64.(class[:cl][!, "TT"]) .* 1e12

            Dls = spectrum_cmb(
                [:TT], prob, jl_tt, ls_tt;
                ptopts = (alg = SymBoltz.Rodas5P(linsolve = SymBoltz.RFLUFactorization()), reltol = 1e-4, abstol = 1e-4),
                sourceopts = (refine = false,),
                normalization = :Dl,
                kτ0s = 0.05*jl_tt.l[begin]:2π/2:1.8*jl_tt.l[end],
                coarse_length = 300,
                verbose = true
            )

            ell_sym = Float64.(collect(ls_tt))
            TT_sym = Float64.(Dls[:, 1] .* 1e12)

            plot!(pltTT, ell_cls, TT_cls,
                  label="CLASS α=$(αval)", color=col, linestyle=:dash)
            plot!(pltTT, ell_sym, TT_sym,
                  label=label, color=col, linestyle=:solid)

            l_lo = max(minimum(ell_cls), minimum(ell_sym))
            l_hi = min(maximum(ell_cls), maximum(ell_sym))

            mask_sym = (ell_sym .>= l_lo) .& (ell_sym .<= l_hi)
            ell_sym2 = ell_sym[mask_sym]
            TT_sym2 = TT_sym[mask_sym]

            itpTT = LinearInterpolation(ell_cls, TT_cls, extrapolation_bc=Line())
            TT_cls_on_sym = itpTT.(ell_sym2)

            relerr = (TT_sym2 .- TT_cls_on_sym) ./ TT_cls_on_sym

            plot!(pltTTErr, ell_sym2, relerr,
                  label="", color=col, linestyle=:solid)
        end

        figTT = plot(pltTT, pltTTErr, layout=(2,1), size=(900,900),
                     left_margin=8mm, bottom_margin=8mm)
        savefig(figTT, OUT_TT)
        println("Saved TT comparison to: ", OUT_TT)

        pltTT = nothing
        pltTTErr = nothing
        figTT = nothing
        GC.gc()
    end

    class_solutions = nothing
    probs = nothing
    M = nothing
    closeall()
    GC.gc()
    GC.gc()
end

main()
GC.gc()
GC.gc()
