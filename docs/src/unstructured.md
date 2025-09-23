# Unstructured ΛCDM model

As an alternative to the modular structure in terms of physical components,
SymBoltz provides a full unstructured ΛCDM model for those that prefer to work with everything in one large system.
This includes *all* variables, parameters and equations in the background, thermodynamics and perturbations.
This makes it **extremely easy** to make direct changes anywhere in the model!

!!! tip
    Paste the code below into a notebook to freely edit the model and interactively inspect changes.

```@example unstructured
using SymBoltz

# constants and some functions
@unpack kB, h, ħ, c, GN, H100, eV, me, mH, mHe, σT, aR, δkron, smoothifelse = SymBoltz
lγmax = 6
ϵ = 1e-9
ΛH = 8.2245809
ΛHe = 51.3
λH∞1s   =  91.17534e-9; fH∞1s  = c/λH∞1s;  EH∞1s  = h*fH∞1s
λH2s1s  = 121.56700e-9; fH2s1s = c/λH2s1s; EH2s1s = h*fH2s1s
EH∞2s  = EH∞1s - EH2s1s # E(∞) - E(2s)
λHe∞1s  =  50.42590e-9; fHe∞1s  = c/λHe∞1s;  EHe∞1s  = h*fHe∞1s
λHe2s1s =  60.14045e-9; fHe2s1s = c/λHe2s1s; EHe2s1s = h*fHe2s1s
λHe2p1s =  58.43344e-9; fHe2p1s = c/λHe2p1s; EHe2p1s = h*fHe2p1s
EHe2p2s = EHe2p1s - EHe2s1s
EHe∞2s  = EHe∞1s - EHe2s1s
EHe⁺∞1s = 54.4178 * eV
λHet∞2s = 260.0463e-9; fHet∞2s = c/λHet∞2s; EHet∞2s = h*fHet∞2s
λHet2p1s = 59.1411e-9; fHet2p1s = c/λHet2p1s; EHet2p1s = h*fHet2p1s
λHet2s1s = 62.5563e-9; fHet2s1s = c/λHet2s1s; EHet2s1s = h*fHet2s1s
EHet2p2s = EHet2p1s - EHet2s1s
A2ps = 1.798287e9
A2pt = 177.58e0
αHfit(T; F=1.125, a=4.309, b=-0.6166, c=0.6703, d=0.5300, T₀=1e4) = F * 1e-19 * a * (T/T₀)^b / (1 + c * (T/T₀)^d)
αHefit(T; q=NaN, p=NaN, T1=10^5.114, T2=3.0) = q / (√(T/T2) * (1+√(T/T2))^(1-p) * (1+√(T/T1))^(1+p))
KHfitfactorfunc(a, A, z, w) = A*exp(-((log(a)+z)/w)^2)
γHe(; A=NaN, σ=NaN, f=NaN) = 3*A*fHe*(1-XHe⁺+ϵ)*c^2 / (8π*σ*√(2π/(β*mHe*c^2))*(1-XH⁺+ϵ)*f^3)
fν = 0

@independent_variables τ # conformal time
D = Differential(τ)

pars = @parameters begin
    k,
    h, # overwrites Planck constant above!
    Ωc0, 
    Ωb0, YHe, fHe,
    Tγ0, Ωγ0,
    ΩΛ0,
    zre1, Δzre1, nre1,
    zre2, Δzre2, nre2
    C,
    τ0
end

vars = @variables begin
    a(τ), z(τ), ℋ(τ), H(τ), Ψ(τ, k), Φ(τ, k), # metric
    ρ(τ), P(τ), δρ(τ, k), Π(τ, k), # gravity
    ρb(τ), Tb(τ), θb(τ, k), δb(τ, k), θb(τ, k), # baryons
    κ(τ), κ̇(τ), csb2(τ), β(τ), ΔT(τ), DTb(τ), DTγ(τ), μc²(τ), Xe(τ), nH(τ), nHe(τ), ne(τ), Xe(τ), ne(τ), λe(τ), Hrec(τ), # recombination
    XH⁺(τ), nH(τ), αH(τ), βH(τ), KH(τ), KHfitfactor(τ), CH(τ) # Hydrogen recombination
    nHe(τ), XHe⁺(τ), XHe⁺⁺(τ), αHe(τ), βHe(τ), RHe⁺(τ), τHe(τ), KHe(τ), invKHe0(τ), invKHe1(τ), invKHe2(τ), CHe(τ), DXHe⁺(τ), DXHet⁺(τ), γ2ps(τ), αHet(τ), βHet(τ), τHet(τ), pHet(τ), CHet(τ), CHetnum(τ), γ2pt(τ), # Helium recombination
    Xre1(τ), Xre2(τ), # reionization
    ργ(τ), Pγ(τ), wγ(τ), Tγ(τ), F0(τ, k), F(τ, k)[1:lγmax], G0(τ, k), G(τ, k)[1:lγmax], δγ(τ, k), θγ(τ, k), σγ(τ, k), Πγ(τ, k) # photons
    ρc(τ), δc(τ, k), θc(τ, k) # cold dark matter
    ρΛ(τ), PΛ(τ), wΛ(τ) # cosmological constant
end

eqs = [
    # metric equations
    z ~ 1/a - 1
    ℋ ~ D(a) / a
    H ~ ℋ / a

    # gravity equations
    D(a) ~ √(8*Num(π)/3 * ρ) * a^2 # 1st Friedmann equation
    D(Φ) ~ -4*Num(π)/3*a^2/ℋ*δρ - k^2/(3*ℋ)*Φ - ℋ*Ψ
    k^2 * (Φ - Ψ) ~ 12*Num(π) * a^2 * Π
    ρ ~ ρc + ρb + ργ + ρΛ
    P ~ Pγ + PΛ
    δρ ~ δc*ρc + δb*ρb + δγ*ργ
    Π ~ (1+wγ) * ργ * σγ

    # baryon recombination
    β ~ 1 / (kB*Tb)
    λe ~ 2π*ħ / √(2π*me/β)
    Hrec ~ H100 * h * H
    fHe ~ YHe / (mHe/mH*(1-YHe))
    D(κ) ~ -a/(H100*h) * ne * σT * c
    κ̇ ~ D(κ)
    csb2 ~ kB/μc² * (Tb - D(Tb)/3ℋ)
    μc² ~ mH*c^2 / (1 + (mH/mHe-1)*YHe + Xe*(1-YHe))
    DTb ~ -2*Tb*ℋ - a/h * 8/3*σT*aR/H100*Tγ^4 / (me*c) * Xe / (1+fHe+Xe) * ΔT
    DTγ ~ D(Tγ)
    D(ΔT) ~ DTb - DTγ
    Tb ~ ΔT + Tγ
    nH ~ (1-YHe) * ρb*(H100*h)^2/GN / mH
    nHe ~ fHe * nH
    ne ~ Xe * nH
    Xe ~ 1*XH⁺ + fHe*XHe⁺ + XHe⁺⁺ + Xre1 + Xre2

    # baryon H⁺ + e⁻ recombination
    αH ~ αHfit(Tb)
    βH ~ αH / λe^3 * exp(-β*EH∞2s)
    KHfitfactor ~ 1 + KHfitfactorfunc(a, -0.14, 7.28, 0.18) + KHfitfactorfunc(a, 0.079, 6.73, 0.33)
    KH ~ KHfitfactor/8π * λH2s1s^3 / Hrec
    CH ~ smoothifelse(XH⁺ - 0.99, (1 + KH*ΛH*nH*(1-XH⁺)) / (1 + KH*(ΛH+βH)*nH*(1-XH⁺)), 1; k = 1e3)
    D(XH⁺) ~ -a/(H100*h) * CH * (αH*XH⁺*ne - βH*(1-XH⁺)*exp(-β*EH2s1s))

    # baryon He⁺ + e⁻ singlet recombination
    αHe ~ αHefit(Tb; q=10^(-16.744), p=0.711)
    βHe ~ 4 * αHe / λe^3 * exp(-β*EHe∞2s)
    KHe ~ 1 / (invKHe0 + invKHe1 + invKHe2)
    invKHe0 ~ 8π*Hrec / λHe2p1s^3
    τHe ~ 3*A2ps*nHe*(1-XHe⁺+ϵ) / invKHe0
    invKHe1 ~ -exp(-τHe) * invKHe0
    γ2ps ~ γHe(A = A2ps, σ = 1.436289e-22, f = fHe2p1s)
    invKHe2 ~ A2ps/(1+0.36*γ2ps^0.86)*3*nHe*(1-XHe⁺)
    CHe ~ smoothifelse(XHe⁺ - 0.99, (exp(-β*EHe2p2s) + KHe*ΛHe*nHe*(1-XHe⁺)) / (exp(-β*EHe2p2s) + KHe*(ΛHe+βHe)*nHe*(1-XHe⁺)), 1; k = 1e3)
    DXHe⁺ ~ -a/(H100*h) * CHe * (αHe*XHe⁺*ne - βHe*(1-XHe⁺)*exp(-β*EHe2s1s))

    # baryon He⁺ + e⁻ triplet recombination
    αHet ~ αHefit(Tb; q=10^(-16.306), p=0.761)
    βHet ~ 4/3 * αHet / λe^3 * exp(-β*EHet∞2s)
    τHet ~ A2pt*nHe*(1-XHe⁺+ϵ)*3 * λHet2p1s^3/(8π*Hrec)
    pHet ~ (1 - exp(-τHet)) / τHet
    γ2pt ~ γHe(A = A2pt, σ = 1.484872e-22, f = fHet2p1s)
    CHetnum ~ A2pt*(pHet+1/(1+0.66*γ2pt^0.9)/3)*exp(-β*EHet2p2s)
    CHet ~ (ϵ + CHetnum) / (ϵ + CHetnum + βHet)
    DXHet⁺ ~ -a/(H100*h) * CHet * (αHet*XHe⁺*ne - βHet*(1-XHe⁺)*3*exp(-β*EHet2s1s))

    # baryon He⁺ + e⁻ total recombination
    D(XHe⁺) ~ DXHe⁺ + DXHet⁺

    # baryon He⁺⁺ + e⁻ recombination
    RHe⁺ ~ 1 * exp(-β*EHe⁺∞1s) / (nH * λe^3)
    XHe⁺⁺ ~ 2*RHe⁺*fHe / (1+fHe+RHe⁺) / (1 + √(1 + 4*RHe⁺*fHe/(1+fHe+RHe⁺)^2))

    # reionization
    Xre1 ~ smoothifelse((1+zre1)^nre1 - (1+z)^nre1, 0, 1 + fHe; k = 1/(nre1*(1+zre1)^(nre1-1)*Δzre1))
    Xre2 ~ smoothifelse((1+zre2)^nre2 - (1+z)^nre2, 0, 0 + fHe; k = 1/(nre2*(1+zre2)^(nre2-1)*Δzre2))

    # baryons
    ρb ~ 3/(8*Num(π)) * Ωb0 * a^(-3)
    D(δb) ~ -θb - 3*ℋ*csb2*δb + 3*D(Φ)
    D(θb) ~ -ℋ*θb + k^2*csb2*δb + k^2*Ψ - 4//3*κ̇*ργ/ρb*(θγ-θb)

    # photons
    Tγ ~ Tγ0 / a
    Ωγ0 ~ π^2/15 * (kB*Tγ0)^4 / (ħ^3*c^5) * 8π*GN / (3*(H100*h)^2)
    ργ ~ 3/(8*Num(π)) * Ωγ0 * a^(-4)
    wγ ~ 1//3
    Pγ ~ wγ * ργ
    D(F0) ~ -k*F[1] + 4*D(Φ)
    D(F[1]) ~ k/3*(F0-2*F[2]+4*Ψ) - 4//3 * κ̇/k * (θb - θγ)
    [D(F[l]) ~ k/(2l+1) * (l*F[l-1] - (l+1)*F[l+1]) + κ̇ * (F[l] - δkron(l,2)//10*Πγ) for l in 2:lγmax-1]...
    D(F[lγmax]) ~ k*F[lγmax-1] - (lγmax+1) / τ * F[lγmax] + κ̇ * F[lγmax]
    δγ ~ F0
    θγ ~ 3*k*F[1]/4
    σγ ~ F[2]/2
    Πγ ~ F[2] + G0 + G[2]
    D(G0) ~ k * (-G[1]) + κ̇ * (G0 - Πγ/2)
    D(G[1]) ~ k/(2*1+1) * (1*G0 - 2*G[2]) + κ̇ * G[1]
    [D(G[l]) ~ k/(2l+1) * (l*G[l-1] - (l+1)*G[l+1]) + κ̇ * (G[l] - δkron(l,2)//10*Πγ) for l in 2:lγmax-1]...
    D(G[lγmax]) ~ k*G[lγmax-1] - (lγmax+1) / τ * G[lγmax] + κ̇ * G[lγmax]

    # cold dark matter
    ρc ~ 3/(8*Num(π)) * Ωc0 * a^(-3)
    D(δc) ~ -(θc-3*D(Φ))
    D(θc) ~ -ℋ*θc + k^2*Ψ

    # cosmological constant
    ρΛ ~ 3/(8*Num(π)) * ΩΛ0
    wΛ ~ -1
    PΛ ~ wΛ * ρΛ
]

initialization_eqs = [
    # metric/gravity
    Ψ ~ 20C / (15 + 4fν)

    # baryons
    δb ~ -3//2 * Ψ
    θb ~ 1//2 * (k^2*τ) * Ψ

    # photons
    F0 ~ -2*Ψ
    F[1] ~ 2//3 * k*τ*Ψ
    F[2] ~ -8//15 * k/κ̇ * F[1]
    [F[l] ~ -l//(2*l+1) * k/κ̇ * F[l-1] for l in 3:lγmax]...
    G0 ~ 5//16 * F[2]
    G[1] ~ -1//16 * k/κ̇ * F[2]
    G[2] ~ 1//16 * F[2]
    [G[l] ~ -l//(2l+1) * k/κ̇ * G[l-1] for l in 3:lγmax]...

    # cold dark matter
    δc ~ -3//2 * Ψ
    θc ~ 1//2 * (k^2*τ) * Ψ
]

defaults = [
    XHe⁺ => 1.0
    XH⁺ => 1.0
    κ => 0.0
    ΔT => 0.0
    zre1 => 7.6711
    Δzre1 => 0.5
    nre1 => 3/2
    zre2 => 3.5
    Δzre2 => 0.5
    nre2 => 1
    C => 1//2
    τ0 => NaN
    D(a) => a / τ
    ΩΛ0 => 1 - Ωγ0 - Ωc0 - Ωb0
]

guesses = [
    a => τ
]

M = System(eqs, τ, vars, pars; initialization_eqs, defaults, guesses, name = :ΛCDM)
```

Now set parameter values and compile the numerical problem:
```@example unstructured
pars = Dict(h => 0.7, Ωc0 => 0.3, Ωb0 => 0.05, YHe => 0.25, Tγ0 => 2.7)
prob = CosmologyProblem(M, pars)
```

Now solve it for some wavenumbers:
```@example unstructured
ks = [4, 40, 400, 4000]
sol = solve(prob, ks)
```

Now plot the evolution of the variables you are interested in:
```@example unstructured
using Plots
p = plot(layout = (3, 1), size = (800, 1000))
plot!(p[1], sol, τ, a)
plot!(p[2], sol, log10(a), [Xe, XH⁺, XHe⁺, XHe⁺⁺]; legend_position = :left)
plot!(p[3], sol, log10(a), [Φ, Ψ], ks)
```

!!! warning
    This unstructured model is provided as a secondary way to work with SymBoltz in the hope that it is useful.
    It may not always be up-to-date with the equivalent structured model built from physical components, which is the primary and best supported interface to the code.
    Because it is unstructured, this model is only available for the standard ΛCDM model.
    It cannot and will not include all combinations of extended models in a way that scales well in model space.
    To enable and disable components, you have to manually edit them in and out yourself.
    For this it can be helpful to look at SymBoltz' library of components in the [src/models/](https://github.com/hersle/SymBoltz.jl/tree/main/src/models) directory.

