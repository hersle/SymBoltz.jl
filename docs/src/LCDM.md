# Minimal ΛCDM model

In addition to [component-based models](@ref "Cosmologies (full models)"),
SymBoltz provides a minimal ΛCDM model for those that prefer everything in one large system.
It includes all variables, parameters and equations in the background, thermodynamics and perturbations,
and makes it **very easy** to modify the model!

!!! tip "Modifying the model"
    Copy-paste the code below, for example into a Jupyter notebook.
    Modify the equations as you wish, but look out for collisions between variable names!
    Comment out equations and species to disable them.
    Equations to generalize the cosmological constant to $w_0 w_a$ dark energy are included below, but commented by default.
    You can also adapt other components already implemented in the [`src/models/`](https://github.com/hersle/SymBoltz.jl/tree/main/src/models) library in SymBoltz.

```@example LCDM
using SymBoltz

# Constants, some functions and atomic energy levels defined in internal files
@unpack kB, ħ, c, GN, H100, eV, me, mH, mHe, σT, aR, δkron, smoothifelse, λH2s1s, EH2s1s, EH∞2s, EHe2s1s, λHe2p1s, fHe2p1s, EHe2p2s, EHe∞2s, EHe⁺∞1s, EHet∞2s, λHet2p1s, fHet2p1s, EHet2s1s, EHet2p2s = SymBoltz
lγmax = 6
lνmax = 6
lhmax = 6
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
    As, ns # primordial power spectrum
end

# 3) Background (τ) and perturbation (τ,k) variables (add your own)
vars = @variables begin
    a(τ), z(τ), ℋ(τ), H(τ), Ψ(τ,k), Φ(τ,k), χ(τ), # metric
    ρ(τ), P(τ), δρ(τ,k), Π(τ,k), # gravity
    ρb(τ), Tb(τ), δb(τ,k), Δb(τ,k), θb(τ,k), # baryons
    κ(τ), _κ(τ), v(τ), csb2(τ), β(τ), ΔT(τ), DTb(τ), μc²(τ), Xe(τ), nH(τ), nHe(τ), ne(τ), Xe(τ), ne(τ), λe(τ), HSI(τ), # recombination
    XH⁺(τ), nH(τ), αH(τ), βH(τ), KH(τ), KHfitfactor(τ), CH(τ), # Hydrogen recombination
    nHe(τ), XHe⁺(τ), XHe⁺⁺(τ), αHe(τ), βHe(τ), RHe⁺(τ), τHe(τ), KHe(τ), invKHe0(τ), invKHe1(τ), invKHe2(τ), CHe(τ), DXHe⁺(τ), DXHet⁺(τ), γ2ps(τ), αHet(τ), βHet(τ), τHet(τ), pHet(τ), CHet(τ), CHetnum(τ), γ2pt(τ), # Helium recombination
    Xre1(τ), Xre2(τ), # reionization
    ργ(τ), Pγ(τ), wγ(τ), Tγ(τ), Fγ0(τ,k), Fγ(τ,k)[1:lγmax], Gγ0(τ,k), Gγ(τ,k)[1:lγmax], δγ(τ,k), θγ(τ,k), σγ(τ,k), Πγ(τ,k), # photons
    ρc(τ), δc(τ,k), Δc(τ,k), θc(τ,k), # cold dark matter
    ρν(τ), Pν(τ), wν(τ), Tν(τ), Fν0(τ,k), Fν(τ,k)[1:lνmax], δν(τ,k), θν(τ,k), σν(τ,k), # massless neutrinos
    ρh(τ), Ph(τ), wh(τ), Ωh(τ), Th(τ), yh(τ), csh2(τ,k), δh(τ,k), Δh(τ,k), σh(τ,k), uh(τ,k), θh(τ,k), Eh(τ)[1:nx], ψh0(τ,k)[1:nx], ψh(τ,k)[1:nx,1:lhmax], Iρh(τ), IPh(τ), Iδρh(τ,k), # massive neutrinos
    ρΛ(τ), PΛ(τ), wΛ(τ), cΛa2(τ), δΛ(τ,k), θΛ(τ,k), ΔΛ(τ,k), # dark energy (cosmological constant or w0wa)
    fν(τ), # misc
    ρm(τ,k), Δm(τ,k), # matter source functions
    ST_SW(τ,k), ST_ISW(τ,k), ST_Doppler(τ,k), ST_polarization(τ,k), ST(τ,k), SE_kχ²(τ,k), Sψ(τ,k) # CMB source functions
end

# 4) Equations for time evolution (modify or add your own)
# TODO: warn/error if equations include variable not in vars to catch common mistakes # hide
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
    ρb ~ 3/8π * Ωb0 / a^3
    D(δb) ~ -θb - 3ℋ*csb2*δb + 3*D(Φ)
    D(θb) ~ -ℋ*θb + k^2*csb2*δb + k^2*Ψ - 4/3*D(κ)*ργ/ρb*(θγ-θb)
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
    ρc ~ 3/8π * Ωc0 / a^3
    D(δc) ~ -θc + 3*D(Φ)
    D(θc) ~ -ℋ*θc + k^2*Ψ
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
    ρΛ ~ 3/8π*ΩΛ0 * abs(a)^(-3*(1+w0+wa)) * exp(-3wa*(1-a))
    PΛ ~ wΛ*ρΛ
    δΛ ~ 0 # for CC
    θΛ ~ 0 # for CC
    #cΛa2 ~ wΛ - D(wΛ)/(3ℋ*(1+wΛ)) # for w0wa
    #D(δΛ) ~ -(1+wΛ)*(θΛ-3*D(Φ)) - 3ℋ*(cΛs2-wΛ)*δΛ - 9*(ℋ/k)^2*(1+wΛ)*(cΛs2-cΛa2)*θΛ # for w0wa
    #D(θΛ) ~ -ℋ*(1-3*cΛs2)*θΛ + cΛs2/(1+wΛ)*k^2*δΛ + k^2*Ψ # for w0wa
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
    Sψ ~ 0 # ifelse(τ ≥ τrec, -(g.Ψ+g.Φ) * (τ-τrec)/(τ0-τrec)/(τ0-τ), 0) # TODO # hide
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
    #δΛ ~ -3/2 * (1+wΛ) * Ψ # for w0wa
    #θΛ ~ 1/2 * (k^2*τ) * Ψ # for w0wa
]

# 6) Initial guess for variables solved for in initial conditions (modify or add your own)
# TODO: also use for shooting method # hide
guesses = [
    a => τ # a(τini) is solved for in a nonlinear system constrained to ℋ(aini) ~ 1/τini (see initialization_eqs)
]

# 7) Default numerical values for parameters and initial conditions (modify or add your own, remove to require explicit value when creating CosmologyProblem)
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
    w0 => -1
    wa => 0
    cΛs2 => 1
]

# 8) Pack everything down into a symbolic system (modify the name to fit your modified model)
M = complete(System(eqs, τ, vars, pars; initialization_eqs, initial_conditions, guesses, name = :ΛCDM))
```

Now set remaining parameter values and compile the numerical problem:
```@example LCDM
p = Dict(
    M.h => 0.7,
    M.Ωc0 => 0.3,
    M.Ωb0 => 0.05,
    M.YHe => 0.25,
    M.Tγ0 => 2.7,
    M.Neff => 3.046,
    M.mh_eV => 0.02,
    M.As => 2e-9,
    M.ns => 1.0,
)
prob = CosmologyProblem(M, p; jac = true, sparse = true) # analytical+sparse Jacobian
```

Now solve it for some wavenumbers:
```@example LCDM
ks = [4e0, 4e1, 4e2, 4e3]
sol = solve(prob, ks)
```

Now plot the evolution of the variables you are interested in:
```@example LCDM
using Plots
p = plot(layout = (3, 1), size = (800, 1000))
plot!(p[1], sol, τ, a)
plot!(p[2], sol, log10(M.a), [M.Xe, M.XH⁺, M.XHe⁺, M.XHe⁺⁺]; legend_position = :left)
plot!(p[3], sol, log10(M.a), [M.Φ, M.Ψ], ks)
```

Now compute the matter power spectrum:
```@example LCDM
modes = [:bc, :m, :h]
ks = 10 .^ range(-1, 4, length=100)
Ps = spectrum_matter(modes, prob, ks)
plot(log10.(ks), log10.(transpose(Ps)), xlabel = "log10(k / (H₀/c))", ylabel = "log10(P / (H₀/c)⁻³)", ylims = (-10, -6), label = permutedims(string.(modes)))
```

Now compute the CMB power spectrum:
```@example LCDM
jl = SphericalBesselCache(25:25:3000)
ls = 25:3000
modes = [:TT, :EE, :TE]
Dls = spectrum_cmb(modes, prob, jl, ls; normalization = :Dl)
plot(ls, Dls[:,1]*1e12, ylabel = "10¹² D(ℓ)", label = "TT", subplot = 1, color = 1, layout = (3, 1), size = (600, 1000), left_margin=5*Plots.mm)
plot!(ls, Dls[:,2]*1e12, ylabel = "10¹² D(ℓ)", label = "EE", subplot = 2, color = 2)
plot!(ls, Dls[:,3]*1e12, ylabel = "10¹² D(ℓ)", label = "TE", subplot = 3, color = 3, xlabel = "ℓ")
```
