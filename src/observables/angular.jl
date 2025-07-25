using Integrals
using Bessels: besselj!, sphericalbesselj
using DataInterpolations
using MatterPower
using ForwardDiff
using ForwardDiffChainRules
import ChainRulesCore

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

# TODO: line-of-sight integrate Θl using ODE for evolution of Jl?
# TODO: spline sphericalbesselj for each l, from x=0 to x=kmax*(τ0-τini)
# TODO: integrate with ApproxFun? see e.g. https://discourse.julialang.org/t/evaluate-integral-on-many-points-cubature-jl/1723/2
# TODO: RombergEven() works with 513 or 1025 points (do Logging.disable_logging(Logging.Warn) first)
# TODO: gaussian quadrature with weight function? https://juliamath.github.io/QuadGK.jl/stable/weighted-gauss/
# line of sight integration
# TODO: use u = k*χ as integration variable, so oscillations of Bessel functions are the same for every k?
# TODO: define and document symbolic dispatch!
"""
    los_integrate(Ss::AbstractMatrix{T}, ls::AbstractVector, τs::AbstractVector, ks::AbstractVector, Rl::Function; integrator = TrapezoidalRule(), verbose = false) where {T <: Real}

For the given `ls` and `ks`, compute the line-of-sight-integrals
```math
Iₗ(k) = ∫dτ S(k,τ) Rₗ(k(τ₀-τ))
```
over the source function values `Ss` against the radial functions `Rl` (e.g. the spherical Bessel functions ``jₗ(x)``).
The element `Ss[i,j]` holds the source function value ``S(kᵢ, τⱼ)``.
"""
function los_integrate(Ss::AbstractMatrix{T}, ls::AbstractVector, τs::AbstractVector, ks::AbstractVector, Rl::Function = jl; integrator = TrapezoidalRule(), verbose = false) where {T <: Real}
    # Julia is column-major; make sure innermost loop indices appear first in slice expressions (https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-column-major)
    @assert size(Ss) == (length(τs), length(ks)) "size(Ss) = $(size(Ss)) ≠ $((length(τs), length(ks)))"
    τs = collect(τs) # force array to avoid floating point errors with ranges in following χs due to (e.g. tiny negative χ)
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
function los_integrate(sol::CosmologySolution, ls::AbstractVector, τs::AbstractVector, ks::AbstractVector, S, Rl::Function = jl; ktransform = identity, kwargs...) # TODO: Ss
    Ss = [S]
    Ss = source_grid(sol, Ss, τs)
    Ss = source_grid(Ss, sol.ks, ks; ktransform)
    Ss = @view Ss[1, :, :]
    return los_integrate(Ss, ls, τs, ks, Rl; kwargs...)
end

"""
    los_temperature(sol::CosmologySolution, ls::AbstractVector, τs::AbstractVector, ks::AbstractVector; ktransform = identity, kwargs...)

Calculate photon temperature multipoles today by line-of-sight integration.
"""
function los_temperature(sol::CosmologySolution, ls::AbstractVector, τs::AbstractVector, ks::AbstractVector; ktransform = identity, kwargs...)
    return los_integrate(sol, ls, τs, ks, sol.prob.M.ST0, jl)
end

"""
    los_polarization(sol::CosmologySolution, ls::AbstractVector, τs::AbstractVector, ks::AbstractVector; ktransform = identity, kwargs...)

Calculate photon E-mode polarization multipoles today by line-of-sight integration.
"""
function los_polarization(sol::CosmologySolution, ls::AbstractVector, τs::AbstractVector, ks::AbstractVector; ktransform = identity, kwargs...)
    return los_integrate(sol, ls, τs, ks, sol.prob.M.ST2_polarization, jl_x2; kwargs...) .* transpose(@. √((ls+2)*(ls+1)*(ls+0)*(ls-1)))
end

# TODO: integrate splines instead of trapz! https://discourse.julialang.org/t/how-to-speed-up-the-numerical-integration-with-interpolation/96223/5
# TODO: better name?
@doc raw"""
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
    spectrum_cmb(modes::AbstractVector, prob::CosmologyProblem, ls::AbstractVector; normalization = :Cl, unit = nothing, Δkτ0 = 2π/2, Δkτ0_S = 8.0, kτ0min = 0.1*ls[begin], kτ0max = 3*ls[end], u = (τ->tanh(τ)), u⁻¹ = (u->atanh(u)), Nlos = 768, integrator = TrapezoidalRule(), bgopts = (alg = Rodas4P(), reltol = 1e-9, abstol = 1e-9), ptopts = (alg = KenCarp4(), reltol = 1e-8, abstol = 1e-8), thread = true, kwargs...)

Compute the CMB power spectra `modes` (`:TT`, `:EE`, `:TE` or an array thereof) ``C_l^{AB}``'s at angular wavenumbers `ls` from the cosmological solution `sol`.
If `unit` is `nothing` the spectra are of dimensionless temperature fluctuations relative to the present photon temperature; while if `unit` is a temperature unit the spectra are of dimensionful temperature fluctuations.
"""
function spectrum_cmb(modes::AbstractVector, prob::CosmologyProblem, ls::AbstractVector; normalization = :Cl, unit = nothing, Δkτ0 = 2π/2, Δkτ0_S = 8.0, kτ0min = 0.1*ls[begin], kτ0max = 3*ls[end], u = (τ->tanh(τ)), u⁻¹ = (u->atanh(u)), Nlos = 768, integrator = TrapezoidalRule(), bgopts = (alg = Rodas4P(), reltol = 1e-9, abstol = 1e-9), ptopts = (alg = KenCarp4(), reltol = 1e-8, abstol = 1e-8), thread = true, kwargs...)
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
        τs[begin] = 0.0
        τs[end] = τ0
    end

    # Integrate perturbations to calculate source function on coarse k-grid
    iT = 'T' in join(modes) ? 1 : 0
    iE = 'E' in join(modes) ? iT + 1 : 0
    Ss = Num[]
    iT > 0 && push!(Ss, prob.M.ST0)
    iE > 0 && push!(Ss, prob.M.ST2_polarization)
    Ss_coarse = source_grid(prob, Ss, τs, ks_coarse; bgopts, ptopts, thread) # TODO: pass kτ0 and x # TODO: pass bgsol

    # Interpolate source function to finer k-grid
    ks_fine = collect(kτ0s_fine ./ τ0)
    ks_fine = clamp.(ks_fine, ks_coarse[begin], ks_coarse[end]) # TODO: ideally avoid
    Ss_fine = source_grid(Ss_coarse, ks_coarse, ks_fine)

    ΘlTs = iT > 0 ? los_integrate(@view(Ss_fine[iT, :, :]), ls, τs, ks_fine, jl; integrator, kwargs...) : nothing
    ΘlEs = iE > 0 ? los_integrate(@view(Ss_fine[iE, :, :]), ls, τs, ks_fine, jl_x2; integrator, kwargs...) .* transpose(@. √((ls+2)*(ls+1)*(ls+0)*(ls-1))) : nothing

    P0s = spectrum_primordial(ks_fine, sol) # more accurate

    if isnothing(unit)
        factor = 1.0 # keep dimensionless
    elseif dimension(unit) == dimension(u"K")
        factor = uconvert(unit, sol[sol.prob.M.γ.T₀] * u"K") # convert to a temperature unit
    else
        error("Requested unit $unit is not a temperature unit")
    end

    spectra = zeros(eltype(Ss_fine[1,1,1] * P0s[1] * factor^2), length(ls), length(modes)) # Cls or Dls
    for (i, mode) in enumerate(modes)
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
        spectra[:, i] .= spectrum
    end

    return spectra
end

"""
    spectrum_cmb(modes::AbstractVector, prob::CosmologyProblem, ls::AbstractVector, ls_interpolate_from::AbstractVector; kwargs...)

Same, but compute the spectrum properly only for `ls_interpolate_from` and then interpolate the results to all `ls`.
"""
function spectrum_cmb(modes::AbstractVector, prob::CosmologyProblem, ls::AbstractVector, ls_interpolate_from::AbstractVector; kwargs...)
    extrema(ls) == extrema(ls_interpolate_from) || error("ls and ls_interpolate_from have different extrema $(extrema(ls)) != $(ls_interpolate_from)")
    spectra_coarse = spectrum_cmb(modes, prob, ls_interpolate_from; kwargs...)
    spectra_fine = similar(spectra_coarse, (length(ls), size(spectra_coarse)[2]))
    for imode in eachindex(modes)
        spectra_fine[:,imode] = QuadraticSpline(spectra_coarse[:,imode], ls_interpolate_from)(ls)
    end
    return spectra_fine
end

function spectrum_cmb(mode::Symbol, args...; kwargs...)
    return spectrum_cmb([mode], args...; kwargs...)[:, begin]
end

function cmb_kτ0s(lmin, lmax; Δkτ0 = 2π/2, Δkτ0_S = 8.0, kτ0min = 0.1*lmin, kτ0max = 3*lmax)
    kτ0s_fine = range(kτ0min, kτ0max, step=Δkτ0) # use integer multiple so endpoints are the same
    kτ0s_coarse = range(kτ0s_fine[begin], kτ0s_fine[end], length = Int(floor((kτ0max-kτ0min)/Δkτ0_S+1)))
    kτ0s_coarse[begin] == kτ0s_fine[begin] && kτ0s_coarse[end] == kτ0s_fine[end] || error("different wavenumber endpoints")
    return kτ0s_coarse, kτ0s_fine
end
