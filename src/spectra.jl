using Integrals
using Bessels: besselj!, sphericalbesselj
using DataInterpolations
using TwoFAST
using MatterPower
using ForwardDiff
using ForwardDiffChainRules
import ChainRulesCore

"""
    spectrum_primordial(k, h, As, ns=1.0; kp = k_dimensionless(0.05 / u"Mpc", h))

Compute the primordial power spectrum
```math
P₀(k) = 2π² Aₛ (k/kₚ)^{nₛ-1} / k³
```
with spectral amplitude `As`, spectral index `ns` and pivot scale wavenumber `kp` at the wavenumber(s) `k`.
"""
function spectrum_primordial(k, h, As, ns=1.0; kp = k_dimensionless(0.05 / u"Mpc", h))
    P = @. 2*π^2 / k^3 * As

    # compute k / kpivot with equal wavenumber units
    k = k_dimensionless.(k, h)
    @. P *= (k/kp)^(ns-1)

    return P
end
function spectrum_primordial(k, sol::CosmologySolution)
    M = sol.prob.M
    return spectrum_primordial(k, sol[M.g.h], sol[M.I.As], sol[M.I.ns])
end
function spectrum_primordial(k, M::System, pars::Dict)
    return spectrum_primordial(k, pars[M.g.h], pars[M.I.As], pars[M.I.ns])
end

"""
    spectrum_matter(sol::CosmologySolution, k, τ = sol[τ][end]; species = [:c, :b, :h])

Compute the power spectrum
```math
P(k,τ) = P₀(k) |Δ(k,τ)|²
```
of the total gauge-invariant overdensity
```math
Δ = δ + (3ℰ/k²) θ = (∑ₛδρ)/(∑ₛρₛ) + (3ℰ/k²) (∑ₛ(ρₛ+Pₛ)θₛ) / (∑ₛ(ρₛ+Pₛ))
```
for the given `species` at wavenumber(s) `k` and conformal time(s) `tau` (final, if omitted) from the solution `sol`.
By default, the species are cold dark matter, baryons and massive neutrinos, which are matter-like at late times in the ΛCDM model.
"""
function spectrum_matter(sol::CosmologySolution, k, τ = sol[τ][end]; species = [:c, :b, :h])
    M = sol.prob.M
    species = getproperty.(M, filter(s -> have(M, s), species))

    δρ = sum(s.ρ*s.δ for s in species)
    ρ = sum(s.ρ for s in species)
    δ = δρ / ρ # total gauge-dependent overdensity

    ρ_plus_P_θ = sum((1+s.w)*s.ρ*s.θ for s in species) # this is the additive perturbation of the energy-momentum tensor
    ρ_plus_P = sum((1+s.w)*s.ρ for s in species) # additive # TODO: meaningful varname for this property?
    θ = ρ_plus_P_θ / ρ_plus_P
    Δ = δ + 3*M.g.ℰ*θ/M.k^2 # total gauge-independent overdensity

    # convert P (through P0) to same units as 1/k^3
    #=
    P0 = sol(k, τ, M.I.P)
    kunit = only(unique(unit.(k)))
    if kunit != NoUnits
        H₀ = sol[M.g.h] * H100
        P0 *= (H₀/c / u"m")^(-3)
        Punit = unit(1 / kunit^3)
        P0 = uconvert.(Punit, P0)
    end
    =#

    P0 = spectrum_primordial(k, sol)
    P = P0 .* sol(k, τ, Δ^2) # Baumann (4.4.172)

    return P
end

"""
    spectrum_matter(prob::CosmologyProblem, k, τ = nothing; species = [:c, :b, :h], kwargs...)

Solve the problem `prob` with exact wavenumber(s) `k`, and then compute the power spectrum with the solution `sol`.
"""
function spectrum_matter(prob::CosmologyProblem, k, τ = nothing; species = [:c, :b, :h], kwargs...)
    sol = solve(prob, k; kwargs...) # TODO: just save endpoints
    τ = isnothing(τ) ? sol[prob.M.τ][end] : τ
    return spectrum_matter(sol, k, τ; species)
end

"""
    spectrum_matter_nonlinear(sol::CosmologySolution, k)

Compute the nonlinear matter power spectrum from the cosmology solution `sol` at wavenumber(s) `k` using halofit implemented in MatterPower.jl.
"""
function spectrum_matter_nonlinear(sol::CosmologySolution, k)
    P = spectrum_matter(sol, k)
    lgPspl = spline(log.(ustrip(P)), log.(ustrip(k)))
    Pf(k) = exp(lgPspl(log(k)))
    halofit_params = setup_halofit(Pf)
    M = sol.prob.M
    Ωm0 = sol[M.c.Ω₀ + M.b.Ω₀] # TODO: generalize to massive neutrinos, redshift etc.
    Pf_halofit(k) = MatterPower.halofit(Pf, halofit_params, Ωm0, ustrip(k))
    PNL = Pf_halofit.(k)

    # if the input add units, multiply it back into the output
    Punit = only(unique(unit.(P)))
    if !isnothing(Punit)
        PNL *= Punit
    end

    return PNL
end

# TODO: generalize to arbitrary field?
"""
    variance_matter(sol::CosmologySolution, R)

Compute the variance ``⟨δ²⟩`` of the *linear* matter density field with a top-hat filter with radius `R`.
Wraps the implementation in MatterPower.jl.
"""
function variance_matter(sol::CosmologySolution, R)
    M = sol.prob.M
    k = sol.ks
    P = spectrum_matter(sol, k)
    lgPspl = spline(log.(P), log.(k))
    Pf(k) = exp(lgPspl(log(k)))
    R = 1 / k_dimensionless(1 / R, sol[M.g.h]) # make dimensionless
    return MatterPower.sigma2(Pf, R)
end
"""
    stddev_matter(sol::CosmologySolution, R)

Compute the standard deviation ``√(⟨δ²⟩)`` of the *linear* matter density field with a top-hat filter with radius `R`.
"""
stddev_matter(sol::CosmologySolution, R) = √(variance_matter(sol, R))

# Out-of-place spherical Bessel function variants
jl(l, x) = sphericalbesselj(l, x) # for l ≥ 0, from Bessels.jl
jl′(l, x) = l/(2l+1)*jl(l-1,x) - (l+1)/(2l+1)*jl(l+1,x) # for l ≥ 1, analytical relation
jl_x2(l, x) = x == 0 ? l == 2 ? 1/15 : 0.0 : jl(l, x) / x^2 # for l ≥ 2, jₗ(x)/x² ≃ xˡ⁻²/(2l+1)!! as x → 0
jl_x2′(l, x) = x == 0 ? l == 3 ? 1/105 : 0.0 : jl′(l,x)/x^2 - 2*jl_x2(l,x)/x # for l ≥ 3; (jₗ(x)/x²)′ = jl′(x)/x² - 2jl(x)/x³ ≃ (l-2)/(2l+1)!! * x^(l-3) as x → 0

# In-place spherical Bessel function variants
# TODO: create SphericalBesselFunctionMachine-like struct, which can calculate value, derivative, ...
# TODO: contribute back to Bessels.jl
# TODO: update with chain rules, reuse jl(x) vector for computing jl′(x) vector
function jl!(out, l::AbstractRange, x::Number)
    besselj!(out, l .+ 0.5, x)
    if x == 0.0 && l[begin] == 0
        out[begin] = 1.0
    elseif x != 0.0
        out .*= √(π/(2*x))
    end
    return out
end
function jlsafe!(out, l::AbstractRange, x::Number)
    out .= jl.(l, x)
    return out
end
function jl′(l, ls::AbstractRange, Jls)
    i = 1 + l - ls[begin] # ls[i] == l (assuming step of ls is 1)
    return l/(2l+1)*Jls[i-1] - (l+1)/(2l+1)*Jls[i+1] # analytical result (see e.g. https://arxiv.org/pdf/astro-ph/9702170 eq. (13)-(15))
end
function jl_x2!(out, l::AbstractRange, x::Number)
    x = max(x, 1e-20) # avoid division by zero
    jl!(out, l, x)
    out ./= x .^ 2
    return out
end

# Overload chain rule for spherical Bessel function
ChainRulesCore.frule((_, _, Δx), ::typeof(jl), l, x) = jl(l, x), jl′(l, x) * Δx # (value, derivative)
ChainRulesCore.frule((_, _, Δx), ::typeof(jl_x2), l, x) = jl_x2(l, x), jl_x2′(l, x) * Δx
@ForwardDiff_frule jl(l::Integer, x::ForwardDiff.Dual) # define dispatch
@ForwardDiff_frule jl_x2(l::Integer, x::ForwardDiff.Dual)

function integrate(xs, ys; integrator = TrapezoidalRule())
    prob = SampledIntegralProblem(ys, xs)
    sol = solve(prob, integrator)
    return sol.u
end

# TODO: line-of-sight integrate Θl using ODE for evolution of Jl?
# TODO: spline sphericalbesselj for each l, from x=0 to x=kmax*(τ0-τini)
# TODO: integrate with ApproxFun? see e.g. https://discourse.julialang.org/t/evaluate-integral-on-many-points-cubature-jl/1723/2
# TODO: RombergEven() works with 513 or 1025 points (do Logging.disable_logging(Logging.Warn) first)
# TODO: gaussian quadrature with weight function? https://juliamath.github.io/QuadGK.jl/stable/weighted-gauss/
# line of sight integration
# TODO: use u = k*χ as integration variable, so oscillations of Bessel functions are the same for every k?
# TODO: define and document symbolic dispatch!
"""
    los_integrate(Ss::AbstractMatrix{T}, ls::AbstractVector, ks::AbstractVector, τs::AbstractVector, Rl::Function; integrator = TrapezoidalRule(), verbose = false) where {T <: Real}

For the given `ls` and `ks`, compute the line-of-sight-integrals
```math
Iₗ(k) = ∫dτ S(k,τ) Rₗ(k(τ₀-τ))
```
over the source function values `Ss` against the radial functions `Rl` (e.g. the spherical Bessel functions ``jₗ(x)``).
The element `Ss[i,j]` holds the source function value ``S(kᵢ, τⱼ)``.
"""
function los_integrate(Ss::AbstractMatrix{T}, ls::AbstractVector, ks::AbstractVector, τs::AbstractVector, Rl::Function = jl; integrator = TrapezoidalRule(), verbose = false) where {T <: Real}
    # Julia is column-major; make sure innermost loop indices appear first in slice expressions (https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-column-major)
    @assert size(Ss) == (length(τs), length(ks))
    χs = τs[end] .- τs
    Is = similar(Ss, length(ks), length(ls))

    # TODO: skip and set Rl to zero if l ≳ kτ0 or another cutoff?
    @tasks for il in eachindex(ls) # parallellize independent loop iterations
        @local begin # define task-local values (declared once for all loop iterations)
            ∂I_∂τ = similar(Ss, length(τs))
        end
        l = ls[il]
        for ik in eachindex(ks)
            k = ks[ik]
            for iτ in eachindex(τs)
                S = Ss[iτ,ik]
                χ = χs[iτ]
                kχ = k * χ
                ∂I_∂τ[iτ] = S * Rl(l, kχ) # TODO: rewrite LOS integral to avoid evaluating Rl with dual numbers? # TODO: rewrite LOS integral with y = kτ0 and x=τ/τ0 to cache jls independent of cosmology
            end
            Is[ik,il] = integrate(τs, ∂I_∂τ; integrator) # integrate over τ # TODO: add starting I(τini) to fix small l?
        end
    end

    return Is
end
function los_integrate(sol::CosmologySolution, ls::AbstractVector, ks::AbstractVector, τs::AbstractVector, S, Rl::Function = jl; ktransform = identity, kwargs...) # TODO: Ss
    Ss = [S]
    Ss = source_grid(sol, Ss, τs)
    Ss = source_grid(Ss, sol.ks, ks; ktransform)
    Ss = @view Ss[:, :, 1]
    return los_integrate(Ss, ls, ks, τs, Rl; kwargs...)
end

"""
    los_temperature(sol::CosmologySolution, ls::AbstractVector, ks::AbstractVector, τs::AbstractVector; ktransform = identity, kwargs...)

Calculate photon temperature multipoles today by line-of-sight integration.
"""
function los_temperature(sol::CosmologySolution, ls::AbstractVector, ks::AbstractVector, τs::AbstractVector; ktransform = identity, kwargs...)
    return los_integrate(sol, ls, ks, τs, sol.prob.M.ST0, jl)
end

"""
    los_polarization(sol::CosmologySolution, ls::AbstractVector, ks::AbstractVector, τs::AbstractVector; ktransform = identity, kwargs...)

Calculate photon E-mode polarization multipoles today by line-of-sight integration.
"""
function los_polarization(sol::CosmologySolution, ls::AbstractVector, ks::AbstractVector, τs::AbstractVector; ktransform = identity, kwargs...)
    return los_integrate(sol, ls, ks, τs, sol.prob.M.ST2_polarization, jl_x2; kwargs...) .* transpose(@. √((ls+2)*(ls+1)*(ls+0)*(ls-1)))
end

# TODO: integrate splines instead of trapz! https://discourse.julialang.org/t/how-to-speed-up-the-numerical-integration-with-interpolation/96223/5
# TODO: better name?
raw"""
    spectrum_cmb(ΘlAs::AbstractMatrix, ΘlBs::AbstractMatrix, P0s::AbstractVector, ls::AbstractVector, ks::AbstractVector; integrator = TrapezoidalRule(), normalization = :Cl)

Compute the angular power spectrum
```math
Cₗᴬᴮ = (2/π) ∫\mathrm{d}k \, k² P₀(k) Θₗᴬ(k,τ₀) Θₗᴮ(k,τ₀)
```
for the given `ls`.
If `normaliation == :Dl`, compute ``Dₗ = Cₗ l (l+1) / 2π`` instead.
"""
function spectrum_cmb(ΘlAs::AbstractMatrix, ΘlBs::AbstractMatrix, P0s::AbstractVector, ls::AbstractVector, ks::AbstractVector; integrator = TrapezoidalRule(), normalization = :Cl)
    size(ΘlAs) == size(ΘlBs) || error("ΘlAs and ΘlBs have different sizes")
    eltype(ΘlAs) == eltype(ΘlBs) || error("ΘlAs and ΘlBs have different types")

    Cls = similar(ΘlAs, length(ls))
    ks_with0 = [0.0; ks] # add dummy value with k=0 for integration

    @tasks for il in eachindex(ls)
        # TODO: skip kτ0 ≲ l?
        @local dCl_dks_with0 = zeros(eltype(ΘlAs), length(ks_with0)) # local task workspace (must zero first element)
        ΘlA = @view ΘlAs[:, il]
        ΘlB = @view ΘlBs[:, il]
        @. dCl_dks_with0[2:end] = 2/π * ks^2 * P0s * ΘlA * ΘlB
        Cls[il] = integrate(ks_with0, dCl_dks_with0; integrator) # integrate over k (_with0 adds one additional point at (0,0))
    end

    if normalization == :Cl
        return Cls
    elseif normalization == :Dl
        return @. Cls * ls * (ls+1) / 2π
    else
        error("Normalization $normalization is not :Cl or :Dl")
    end
end

"""
    spectrum_cmb(modes::AbstractVector, prob::CosmologyProblem, ls::AbstractVector; normalization = :Cl, unit = nothing, Δkτ0 = 2π/2, Δkτ0_S = 8.0, kτ0min = 0.1*ls[begin], kτ0max = 3*ls[end], u = (τ->tanh(τ)), u⁻¹ = (u->atanh(u)), Nlos = 768, integrator = TrapezoidalRule(), bgopts = (alg = Rodas4P(), reltol = 1e-8), ptopts = (alg = KenCarp4(), reltol = 1e-8), thread = true, kwargs...)

Compute the CMB power spectra `modes` (`:TT`, `:EE`, `:TE` or an array thereof) ``C_l^{AB}``'s at angular wavenumbers `ls` from the cosmological solution `sol`.
If `unit` is `nothing` the spectra are of dimensionless temperature fluctuations relative to the present photon temperature; while if `unit` is a temperature unit the spectra are of dimensionful temperature fluctuations.
"""
function spectrum_cmb(modes::AbstractVector, prob::CosmologyProblem, ls::AbstractVector; normalization = :Cl, unit = nothing, Δkτ0 = 2π/2, Δkτ0_S = 8.0, kτ0min = 0.1*ls[begin], kτ0max = 3*ls[end], u = (τ->tanh(τ)), u⁻¹ = (u->atanh(u)), Nlos = 768, integrator = TrapezoidalRule(), bgopts = (alg = Rodas4P(), reltol = 1e-8), ptopts = (alg = KenCarp4(), reltol = 1e-8), thread = true, kwargs...)
    kτ0s_coarse, kτ0s_fine = cmb_kτ0s(ls[begin], ls[end]; Δkτ0, Δkτ0_S, kτ0min, kτ0max)
    sol = solve(prob; bgopts)
    τ0 = getsym(sol, prob.M.τ0)(sol)
    ks_coarse = kτ0s_coarse ./ τ0
    τs = sol.bg.t # by default, use background (thermodynamics) time points for line of sight integration
    if Nlos != 0 # instead choose Nlos time points τ = τ(u) corresponding to uniformly spaced u
        τmin, τmax = extrema(τs)
        umin, umax = u(τmin), u(τmax)
        us = range(umin, umax, length = Nlos)
        τs = u⁻¹.(us)
        τs[end] = τ0
    end

    # Integrate perturbations to calculate source function on coarse k-grid
    iT = 'T' in join(modes) ? 1 : 0
    iE = 'E' in join(modes) ? iT + 1 : 0
    Ss = Num[]
    iT > 0 && push!(Ss, prob.M.ST0)
    iE > 0 && push!(Ss, prob.M.ST2_polarization)
    Ss_coarse = source_grid(prob, Ss, ks_coarse, τs) # TODO: pass kτ0 and x

    # Interpolate source function to finer k-grid
    ks_fine = collect(kτ0s_fine ./ τ0)
    ks_fine = clamp.(ks_fine, ks_coarse[begin], ks_coarse[end]) # TODO: ideally avoid
    Ss_fine = source_grid(Ss_coarse, ks_coarse, ks_fine)

    ΘlTs = iT > 0 ? los_integrate(@view(Ss_fine[:, :, iT]), ls, ks_fine, τs, jl; integrator, kwargs...) : nothing
    ΘlEs = iE > 0 ? los_integrate(@view(Ss_fine[:, :, iE]), ls, ks_fine, τs, jl_x2; integrator, kwargs...) .* transpose(@. √((ls+2)*(ls+1)*(ls+0)*(ls-1))) : nothing

    P0s = spectrum_primordial(ks_fine, sol) # more accurate

    if isnothing(unit)
        factor = 1.0 # keep dimensionless
    elseif dimension(unit) == dimension(u"K")
        factor = uconvert(unit, sol[sol.prob.M.γ.T₀] * u"K") # convert to a temperature unit
    else
        error("Requested unit $unit is not a temperature unit")
    end

    spectra = [] # Cls or Dls # TODO: make multidim array
    for mode in modes
        if mode == :TT
            spectrum = spectrum_cmb(ΘlTs, ΘlTs, P0s, ls, ks_fine; integrator, normalization)
        elseif mode == :EE
            spectrum = spectrum_cmb(ΘlEs, ΘlEs, P0s, ls, ks_fine; integrator, normalization)
        elseif mode == :TE
            spectrum = spectrum_cmb(ΘlTs, ΘlEs, P0s, ls, ks_fine; integrator, normalization)
        else
            error("Unknown CMB power spectrum mode $mode")
        end
        spectrum *= factor^2 # possibly make dimensionful
        push!(spectra, spectrum)
    end

    return spectra
end
spectrum_cmb(mode::Symbol, args...; kwargs...) = only(spectrum_cmb([mode], args...; kwargs...))

function cmb_kτ0s(lmin, lmax; Δkτ0 = 2π/2, Δkτ0_S = 8.0, kτ0min = 0.1*lmin, kτ0max = 3*lmax)
    kτ0s_fine = range(kτ0min, kτ0max, step=Δkτ0) # use integer multiple so endpoints are the same
    kτ0s_coarse = range(kτ0s_fine[begin], kτ0s_fine[end], length = Int(floor((kτ0max-kτ0min)/Δkτ0_S+1)))
    kτ0s_coarse[begin] == kτ0s_fine[begin] && kτ0s_coarse[end] == kτ0s_fine[end] || error("different wavenumber endpoints")
    return kτ0s_coarse, kτ0s_fine
end

"""
    correlation_function(sol::CosmologySolution; N = 2048, spline = true)

Compute the two-point correlation function in real space by Fourier transforming the matter power spectrum of `sol` with `N` points the FFTLog algorithm implemented in TwoFAST.
Returns `N` radii and correlation function values (e.g. `r`, `ξ`).
"""
function correlation_function(sol::CosmologySolution; N = 2048, spline = true)
    ks = sol.ks
    if spline
        P = SymBoltz.spline(spectrum_matter(sol, ks), ks) # create spline interpolation (fast)
    else
        P(k) = only(spectrum_matter(sol, k)) # use solution's built-in interpolation (elegant)
    end
    kmin, kmax = extrema(ks)
    rmin = 2π / kmax
    return xicalc(P, 0, 0; N, kmin, kmax, r0=rmin)
end

# TODO: not really a spectrum...
# TODO: add formula
"""
    distance_luminosity(sol::CosmologySolution, ivs = sol.bg.t, τ0 = sol[sol.prob.M.τ0])

Compute luminosity distances
```math
d_L = \\frac{r}{a} = \\chi \\, \\mathrm{sinc} (\\sqrt{K} (τ₀-τ)),
```
at the independent variable values `ivs` relative to the (present) time `τ0`.
"""
function distance_luminosity(sol::CosmologySolution, ivs = sol.bg.t, τ0 = sol[sol.prob.M.τ0])
    M = sol.prob.M
    χ = sol(ivs, M.χ)
    Ωk0 = have(M, :K) ? sol[M.K.Ω₀] : 0.0
    r = sinc.(√(-Ωk0+0im)*χ/π) .* χ |> real # Julia's sinc(x) = sin(π*x) / (π*x)
    H0 = H100 * sol[M.g.h]
    a = sol(ivs, M.g.a)
    return @. r / a * SymBoltz.c / H0 # to meters
end

# TODO: take in ks and N_skip_interp = 1, 2, 3, ... for more efficient? same for τ?
# TODO: source_grid(prob::CosmologyProblem seemed type-stable?
# TODO: return getter for (ks, τs)
function source_grid(sol::CosmologySolution, Ss::AbstractVector, τs)
    # Evaluate integrated perturbations on coarse grid
    ks = sol.ks
    getSs = map(S -> getsym(sol.prob.pt, S), Ss)
    Ss = similar(sol.bg, length(τs), length(ks), length(Ss))
    @tasks for ik in eachindex(ks)
        for iS in eachindex(getSs)
            Ss[:, ik, :] .= getSs[iS](sol.pts[ik]) # TODO: type infer getS?
        end
    end
    return Ss
end

function source_grid(Ss_coarse::AbstractArray, ks_coarse, ks_fine; ktransform = identity)
    size_coarse = size(Ss_coarse)
    size_fine = (size_coarse[1], length(ks_fine), size_coarse[3])
    Nτ, _, Ns = size(Ss_coarse)

    Ss_fine = similar(Ss_coarse, size_fine)
    xs_coarse = ktransform.(ks_coarse) # TODO: user should just pass different ks as input instead
    xs_fine = ktransform.(ks_fine)
    @tasks for iτ in 1:Nτ
        for iS in 1:Ns
            interp = LinearInterpolation(@view(Ss_coarse[iτ, :, iS]), xs_coarse)
            Ss_fine[iτ, :, iS] .= interp.(xs_fine)
        end
    end
    return Ss_fine
end

# TODO: take in kτ0s and xs
function source_grid(prob::CosmologyProblem, S::AbstractArray, ks, τs)
    bgsol = solvebg(prob.bg)
    getSs = map(s -> getsym(prob.pt, s), S)
    Ss = similar(bgsol, length(τs), length(ks), length(S))
    function output_func(sol, ik)
        for iS in eachindex(getSs)
            Ss[:, ik, iS] .= getSs[iS](sol)
        end
        return nothing, false
    end
    solvept(prob.pt, bgsol, ks, prob.var2spl; output_func, saveat = τs)
    return Ss
end

# TODO: test @inferred
function distance_luminosity_function(M::System, pars_fixed, pars_varying, zs; bgopts = (alg = SymBoltz.Tsit5(), reltol = 1e-5, maxiters = 1e3))
    isequal(independent_variable(M), M.g.a) || error("Independent variable must be $(M.g.a)")

    pars = merge(pars_fixed, Dict(pars_varying .=> NaN))
    as = @. 1 / (zs + 1)
    prob = CosmologyProblem(M, pars; pt = false, ivspan = (minimum(as), 1.0))
    probgen = SymBoltz.parameter_updater(prob, pars_varying; build_initializeprob = Val{false})

    geta = getsym(prob, M.g.a)
    getτ = getsym(prob, M.τ)
    geth = getsym(prob, M.g.h)
    getΩk0 = getsym(prob, M.K.Ω₀)

    return p -> begin
        prob = probgen(p)
        sol = solve(prob; bgopts, saveat = as, save_end = true)
        a = geta(sol)
        τ = getτ(sol)
        h = geth(sol)
        Ωk0 = getΩk0(sol)
        τ0 = τ[end] # time today
        χ = τ0 .- τ
        r = @. real(sinc(√(-Ωk0+0im)*χ/π) * χ) # Julia's sinc(x) = sin(π*x) / (π*x)
        H0 = SymBoltz.H100 * h
        return @. r / a * SymBoltz.c / H0 # luminosity distance in meters
    end
end
