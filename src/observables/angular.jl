using Bessels: besselj!, sphericalbesselj
using DataInterpolations
using MatterPower
using ForwardDiff
using ForwardDiffChainRules
import ChainRulesCore

struct SphericalBesselCache{Tl, Tdy <: Union{Matrix{Float64}, Nothing}}
    l::Tl
    y::Matrix{Float64}
    dy::Tdy
    dx::Float64
    invdx::Float64
    x::Vector{Float64}
end

"""
    SphericalBesselCache(ls; xmax = 20*maximum(ls), dx = 2π/15, hermite = true, thread = true)

Create interpolation cache for the spherical Bessel function ``jₗ(x)`` for orders `ls` for `0 ≤ x ≤ xmax` with uniform spacing `dx`.
If `hermite`, cubic Hermite interpolation is used with the analytical derivative ``jₗ′(x)`` instead of linear interpolation.
The computation uses fast recurrence relations when `ls` contains `Integer` orders only,
and otherwise falls back to explicit ``jₗ(x)`` evaluation for every ``l`` and ``x``.
If `thread`, the tabulation is parallellized over independent ``x``.
"""
function SphericalBesselCache(ls; xmax = 20*maximum(ls), dx = 2π/15, hermite = true, thread = true)
    xmin = 0.0
    xs = range(xmin, xmax, length = trunc(Int, (xmax - xmin) / dx)) # fixed length (so endpoints are exact) that gives step as close to dx as possible
    invdx = 1.0 / step(xs) # using the resulting step, which need not be exactly dx
    xs = collect([xs; xs[end]]) # pad with 1 extra duplicate point to avoid bounds check during interpolation
    ys = Matrix{Float64}(undef, length(ls), length(xs)) # contiguous in l
    dys = hermite ? similar(ys) : nothing
    jl_table!(ys, dys, ls, xs; thread)
    return SphericalBesselCache{typeof(ls), typeof(dys)}(ls, ys, dys, dx, invdx, xs)
end

# First argument is the cache index il, not the multipole l
@inline Base.@propagate_inbounds @fastmath function (jl::SphericalBesselCache{Tl, Nothing})(il::Int, x) where {Tl}
    w = x * jl.invdx # 0-based float index (assume x0 = 0)
    i = trunc(Int, w) # 0-based integer index of left interval point; faster than searchsortedfirst(jl.x, x)
    w = w - i # remainder ∈ [0, 1]
    y₋ = jl.y[il, i+1] # +1 for 1-based indexing
    y₊ = jl.y[il, i+2]
    return muladd(w, y₊ - y₋, y₋) # i.e. y₋ + (y₊ - y₋) * (x - x₋) * jl.invdx
end

@inline Base.@propagate_inbounds @fastmath function (jl::SphericalBesselCache{Tl, Matrix{Float64}})(il::Int, x) where {Tl}
    w = x * jl.invdx
    i = trunc(Int, w)
    w = w - i
    wm1 = w - 1.0
    y₋  = jl.y[il, i+1]
    y₊  = jl.y[il, i+2]
    dy₋ = jl.dy[il, i+1]
    dy₊ = jl.dy[il, i+2]
    return (1+2w)*wm1*wm1 * y₋ + w*w*(3-2w) * y₊ + w*wm1 * (wm1 * dy₋ + w * dy₊) * jl.dx # https://en.wikipedia.org/wiki/Cubic_Hermite_spline
end

function Base.show(io::IO, jl::SphericalBesselCache{Tl, Tdy}) where {Tl, Tdy}
    method = Tdy == Nothing ? "linear" : "Hermite"
    print(io, "jₗ(x) $method interpolation cache ")
    print(io, "for $(minimum(jl.l)) ≤ l ≤ $(maximum(jl.l)) and ")
    print(io, "$(jl.x[begin]) ≤ x ≤ $(jl.x[end]) ")
    print(io, "($(Base.format_bytes(Base.summarysize(jl))))\n")
end

# Out-of-place spherical Bessel function variants
jl(l, x) = sphericalbesselj(l, x) # for l ≥ 0, from Bessels.jl
jl′(l, x) = iszero(l) ? -jl(one(l), x) : l/(2l+1)*jl(l-1,x) - (l+1)/(2l+1)*jl(l+1,x) # analytical relation with special case for j₀′(x) = -j₁(x), where the general relation would evaluate j₋₁, which diverges at x = 0, with zero weight)

# In-place spherical Bessel function variants
# TODO: contribute back to Bessels.jl
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

# Overload chain rule for spherical Bessel function
ChainRulesCore.frule((_, _, Δx), ::typeof(jl), l, x) = jl(l, x), jl′(l, x) * Δx # (value, derivative)
@ForwardDiff_frule jl(l::Integer, x::ForwardDiff.Dual) # define dispatch

# TODO: line-of-sight integrate Θl using ODE for evolution of Jl?
# TODO: spline sphericalbesselj for each l, from x=0 to x=kmax*(τ0-τini)
# TODO: integrate with ApproxFun? see e.g. https://discourse.julialang.org/t/evaluate-integral-on-many-points-cubature-jl/1723/2
# TODO: RombergEven() works with 513 or 1025 points (do Logging.disable_logging(Logging.Warn) first)
# TODO: gaussian quadrature with weight function? https://juliamath.github.io/QuadGK.jl/stable/weighted-gauss/
# line of sight integration
# TODO: use u = k*χ as integration variable, so oscillations of Bessel functions are the same for every k?
# TODO: define and document symbolic dispatch!
"""
    los_integrate(Ss::AbstractMatrix{T}, ls::AbstractVector, τs::AbstractVector, ks::AbstractVector, jl::SphericalBesselCache, K::Real = 0.0; l_limber = typemax(Int), thread = true, verbose = false) where {T}

For the given `ls` and `ks`, compute the line-of-sight integrals
```math
Iₗ(k) = ∫dτ S(k,τ) jₗ(k(τ₀-τ))
```
over the source function values `Ss` against the spherical Bessel functions ``jₗ(x)`` cached in `jl`.
The element `Ss[i,j]` holds the source function value ``S(τᵢ, kⱼ)``.
The Limber approximation
```math
Iₗ ≈ √(π/(2l+1)) S(τ₀-(l+1/2)/k, k)
```
is used for `l ≥ l_limber`.
"""
function los_integrate(Ss::AbstractMatrix{T}, ls::AbstractVector, τs::AbstractVector, ks::AbstractVector, jl::SphericalBesselCache, K::Real = 0.0; l_limber = typemax(Int), thread = true, verbose = false) where {T}
    @assert K == 0 "the spherical Bessel function cache only covers a flat universe (K = 0), but got K = $K; pass the multipoles `ls` instead of a `SphericalBesselCache` to use the hyperspherical recursion"
    @assert size(Ss, 1) == length(τs) "size(Ss, 1) = $(size(Ss, 1)) and length(τs) = $(length(τs)) differ"
    @assert size(Ss, 2) == length(ks) "size(Ss, 2) = $(size(Ss, 2)) and length(ks) = $(length(ks)) differ"
    @assert collect(ls) == collect(jl.l) "ls must match the l-values stored in the Bessel cache"
    @assert jl.x[begin] ≤ 0 "jl.x[begin] < 0"
    @assert jl.x[end] ≥ ks[end]*τs[end] "jl.x[end] < kmax*τmax"
    @assert issorted(τs) "τs must be sorted in ascending order"
    @assert issorted(ks) "ks must be sorted in ascending order"
    @assert issorted(ls) "ls must be sorted in ascending order" # necessary for Limber indexing logic
    error_if_nonfinite(Ss)

    τs = collect(τs) # force array to avoid floating point errors with ranges in following χs due to (e.g. tiny negative χ)
    τ0 = τs[end]
    χs = τ0 .- τs
    nτ = length(τs)

    ws = similar(τs) # precompute trapezoidal rule weights
    ws[1] = 0.5 * (τs[2] - τs[1])
    @inbounds for iτ in 2:nτ-1
        ws[iτ] = 0.5 * (τs[iτ+1] - τs[iτ-1])
    end
    ws[nτ] = 0.5 * (τs[nτ] - τs[nτ-1])

    nl = length(ls)
    Is = similar(Ss, length(ks), nl)
    il_limber = searchsortedfirst(ls, l_limber) # First il index with l ≥ l_limber (=nl+1 when l_limber = typemax, i.e. no Limber modes)

    verbose && l_limber < typemax(Int) && println("Using Limber approximation for l ≥ $l_limber")

    # Loop order k → τ → l to get SIMD on the innermost l-loop
    @fastmath @inbounds @tasks for ik in eachindex(ks)
        @set scheduler = thread ? :dynamic : :serial
        @local tmp = zeros(T, nl) # l-contiguous storage for integrals (to help SIMD over l)
        k = ks[ik]
        verbose && print("\rLOS integrating k-mode $ik / $(length(ks))")

        # Full line-of-sight integrals for l < l_limber
        fill!(tmp, zero(T))
        @inbounds for iτ in eachindex(τs)
            kχ = k * χs[iτ]
            Sw = ws[iτ] * Ss[iτ, ik]
            @inbounds @simd for il in 1:il_limber-1
                tmp[il] += Sw * jl(il, kχ)
            end
        end

        # Limber approximation for l ≥ l_limber
        @inbounds for il in il_limber:nl
            l = ls[il]
            χ = (l + 1/2) / k
            if χ ≤ χs[1] # otherwise source is zero before recombination
                i₋ = searchsortedfirst(τs, τ0 - χ)
                χ₋ = χs[i₋]
                S₋ = Ss[i₋, ik]
                if i₋ == 1
                    S = S₋
                else
                    i₊ = i₋ - 1 # χs is descending, so χ₋ < χ < χ₊
                    χ₊ = χs[i₊]
                    S₊ = Ss[i₊, ik]
                    Δχ = χ₊ - χ₋
                    S′₋ = i₋ ≤ nτ-1 ? (Ss[i₋+1, ik] - S₊) / (χs[i₋+1] - χ₊) : (S₊ - S₋) / Δχ
                    S′₊ = i₊ ≥ 2    ? (S₋ - Ss[i₋-2, ik]) / (χ₋ - χs[i₋-2]) : (S₊ - S₋) / Δχ
                    t = (χ - χ₋) / Δχ
                    t² = t*t
                    t³ = t²*t
                    S = (2t³-3t²+1)*S₋ + (t³-2t²+t)*Δχ*S′₋ + (-2t³+3t²)*S₊ + (t³-t²)*Δχ*S′₊
                    tmp[il] = √(π/(2l+1)) * S / k
                end
            end
        end

        Is[ik, :] .= tmp
    end
    verbose && println()

    return Is
end

"""
    los_integrate(Ss::AbstractMatrix{T}, ls::AbstractVector, τs::AbstractVector, ks::AbstractVector, ::Nothing, K::Real = 0.0; kwargs...) where {T}

Same as `los_integrate(Ss, ls, τs, ks, K; ...)`. This lets callers that hold a `jl` which is either a
`SphericalBesselCache` or `nothing` dispatch on it directly, instead of branching on `isnothing(jl)`.
"""
function los_integrate(Ss::AbstractMatrix{T}, ls::AbstractVector, τs::AbstractVector, ks::AbstractVector, ::Nothing, K::Real = 0.0; kwargs...) where {T}
    return los_integrate(Ss, ls, τs, ks, K; kwargs...)
end

"""
    los_integrate(Ss::AbstractMatrix{T}, ls::AbstractVector, τs::AbstractVector, ks::AbstractVector, K::Real = 0.0; l_limber = typemax(Int), thread = true, verbose = false) where {T}

Same as `los_integrate(Ss, ls, τs, ks, jl::SphericalBesselCache; ...)`, but compute the hyperspherical Bessel function
``Φₗ^k(χ)`` on the fly with recurrences for integer `ls` (without a precomputed cache) and with curvature ``K``.
"""
function los_integrate(Ss::AbstractMatrix{T}, ls::AbstractVector, τs::AbstractVector, ks::AbstractVector, K::Real = 0.0; l_limber = typemax(Int), thread = true, verbose = false) where {T}
    K ≤ 0 || throw(ArgumentError("line-of-sight integration is only implemented for flat (K = 0) and open (K < 0) universes, but got K = $K (closed)"))
    @assert size(Ss, 1) == length(τs) "size(Ss, 1) = $(size(Ss, 1)) and length(τs) = $(length(τs)) differ"
    @assert size(Ss, 2) == length(ks) "size(Ss, 2) = $(size(Ss, 2)) and length(ks) = $(length(ks)) differ"
    @assert issorted(τs) "τs must be sorted in ascending order"
    @assert issorted(ks) "ks must be sorted in ascending order"
    @assert issorted(ls) "ls must be sorted in ascending order" # necessary for Limber indexing logic
    all(l -> l == round(l), ls) || throw(ArgumentError(
        "ls must be exact integer multipoles for the recursion-based los_integrate (got e.g. $(first(Iterators.filter(l -> l != round(l), ls)))). " *
        "A coarse/fractional l-grid for interpolation is only supported by the SphericalBesselCache-based method."))
    ls = round.(Int, ls)
    error_if_nonfinite(Ss)

    τs = collect(τs) # force array to avoid floating point errors with ranges in following χs due to (e.g. tiny negative χ)
    τ0 = τs[end]
    χs = τ0 .- τs
    nτ = length(τs)

    ws = similar(τs) # precompute trapezoidal rule weights
    ws[1] = 0.5 * (τs[2] - τs[1])
    @inbounds for iτ in 2:nτ-1
        ws[iτ] = 0.5 * (τs[iτ+1] - τs[iτ-1])
    end
    ws[nτ] = 0.5 * (τs[nτ] - τs[nτ-1])

    nl = length(ls)
    Is = similar(Ss, length(ks), nl)
    il_limber = searchsortedfirst(ls, l_limber) # First il index with l ≥ l_limber (=nl+1 when l_limber = typemax, i.e. no Limber modes)
    ls_nonlimber = @view ls[1:il_limber-1] # only these need the (non-Limber) recursion; Φl is sized/indexed to match
    lmax = il_limber > 1 ? ls[il_limber-1] : 0 # largest l that needs the full (non-Limber) recursion

    verbose && l_limber < typemax(Int) && println("Using Limber approximation for l ≥ $l_limber")

    il_limber == 1 || lmax ≥ 2 || throw(ArgumentError("the recursion needs multipoles l ≥ 2, but got lmax = $lmax below l_limber = $l_limber"))

    @fastmath @inbounds @tasks for ik in eachindex(ks)
        @set scheduler = thread ? :dynamic : :serial
        @local begin
            tmp = zeros(T, nl)
            Φl = zeros(Float64, nτ, length(ls_nonlimber)) # Φₗ(χᵢ) for every τᵢ and l ≤ lmax at once (see Φl_recurrence!)
            sqK = zeros(Float64, lmax+2) # recursion coefficients √Kₗ = √(k²-Kl²) and 1/√Kₗ for l = 0, …, lmax+1,
            invsqK = zeros(Float64, lmax+2) # which depend on k but not on χ (so are tabulated once per k)
        end
        k = ks[ik]
        sqrtK_table!(sqK, invsqK, K, k)
        verbose && print("\rLOS integrating k-mode $ik / $(length(ks))")

        # Full line-of-sight integrals for l < l_limber (skipped entirely when every requested l uses Limber)
        fill!(tmp, zero(T))
        if il_limber > 1
            Φl_recurrence!(Φl, ls_nonlimber, χs, k, K, sqK, invsqK) # one recursion sweep gives Φₗ(χᵢ) for every τᵢ and l ≤ lmax at once
            @inbounds for iτ in eachindex(τs)
                Sw = ws[iτ] * Ss[iτ, ik]
                @inbounds @simd for il in 1:il_limber-1
                    tmp[il] += Sw * Φl[iτ, il]
                end
            end
        end

        # Limber approximation for l ≥ l_limber (does not use Φₗ at all). Φₗ acts like a delta function at its
        # turning point, which now sits at k⋅sinK(K,χ) = l+½ instead of at kχ = l+½. Matching turning points is
        # exactly CLASS' flat approximation Φₗ(χ) ≈ A⋅jₗ(γkχ) with γk = (l+½)/χₗ (transfer.c, `rescale_argument`
        # and `rescale_amplitude`), so the flat Limber result carries over with 1/k → A⋅χₗ/(l+½).
        @inbounds for il in il_limber:nl
            l = ls[il]
            χ = asinK(K, (l + 1/2) / k) # turning point χₗ, where k⋅sinK(K,χ) = l+½
            if χ ≤ χs[1] # otherwise source is zero before recombination
                i₋ = searchsortedfirst(τs, τ0 - χ)
                χ₋ = χs[i₋]
                S₋ = Ss[i₋, ik]
                if i₋ == 1
                    S = S₋
                else
                    i₊ = i₋ - 1 # χs is descending, so χ₋ < χ < χ₊
                    χ₊ = χs[i₊]
                    S₊ = Ss[i₊, ik]
                    Δχ = χ₊ - χ₋
                    S′₋ = i₋ ≤ nτ-1 ? (Ss[i₋+1, ik] - S₊) / (χs[i₋+1] - χ₊) : (S₊ - S₋) / Δχ
                    S′₊ = i₊ ≥ 2    ? (S₋ - Ss[i₋-2, ik]) / (χ₋ - χs[i₋-2]) : (S₊ - S₋) / Δχ
                    t = (χ - χ₋) / Δχ
                    t² = t*t
                    t³ = t²*t
                    S = (2t³-3t²+1)*S₋ + (t³-2t²+t)*Δχ*S′₋ + (-2t³+3t²)*S₊ + (t³-t²)*Δχ*S′₊
                    A = (1 - K*l*(l+1)/k^2)^(-1/12) # amplitude of the same rescaling (CLASS' `rescale_amplitude`); → 1 when K = 0
                    tmp[il] = √(π/(2l+1)) * S * A * χ / (l + 1/2) # → the flat √(π/(2l+1))⋅S/k, since χ = (l+½)/k there
                end
            end
        end

        Is[ik, :] .= tmp
    end
    verbose && println()

    return Is
end

# TODO: integrate splines instead of trapz! https://discourse.julialang.org/t/how-to-speed-up-the-numerical-integration-with-interpolation/96223/5
@doc raw"""
    spectrum_cmb(ΘlAs::AbstractMatrix, ΘlBs::AbstractMatrix, P0s::AbstractVector, ls::AbstractVector, ks::AbstractVector; normalization = :Cl, thread = true)

Compute the angular power spectrum
```math
Cₗᴬᴮ = (2/π) ∫\mathrm{d}k \, k² P₀(k) Θₗᴬ(k,τ₀) Θₗᴮ(k,τ₀)
```
for the given `ls`.
If `normalization == :Dl`, compute ``Dₗ = Cₗ l (l+1) / 2π`` instead.
"""
function spectrum_cmb(ΘlAs::AbstractMatrix, ΘlBs::AbstractMatrix, P0s::AbstractVector, ls::AbstractVector, ks::AbstractVector; normalization = :Cl, thread = true)
    size(ΘlAs) == size(ΘlBs) || error("ΘlAs and ΘlBs have different sizes")
    eltype(ΘlAs) == eltype(ΘlBs) || error("ΘlAs and ΘlBs have different types")

    Cls = similar(ΘlAs, length(ls))
    ks_with0 = [0.0; ks] # add dummy value with k=0 for integration

    @tasks for il in eachindex(ls)
        # TODO: skip kτ0 ≲ l?
        @set scheduler = thread ? :dynamic : :static
        @local dCl_dks_with0 = zeros(eltype(ΘlAs), length(ks_with0)) # local task workspace (must zero first element)
        ΘlA = @view ΘlAs[:, il]
        ΘlB = @view ΘlBs[:, il]
        @. dCl_dks_with0[2:end] = 2/π * ks^2 * P0s * ΘlA * ΘlB
        spline = CubicSpline(dCl_dks_with0, ks_with0)
        Cls[il] = DataInterpolations.integral(spline, ks_with0[begin], ks_with0[end]) # integrate over k (_with0 adds one additional point at (0,0))
    end

    return normalize_spectrum_cmb(normalization, ls, Cls)
end

fk_tanh(k, k0=2000.0) = tanh(k/k0)
fk⁻¹_tanh(k, k0=2000.0) = k0*atanh(k)

"""
    spectrum_cmb(modes::AbstractVector{<:Symbol}, prob::CosmologyProblem, jl::SphericalBesselCache; normalization = :Cl, unit = nothing, kinterp = nothing, Δkτ0 = 2π/4, xs = cosgrid(0.0, 1.0; length=300), τcut = 1e-2, l_limber = 10, bgopts = (alg = bgalg(prob), reltol = 1e-7, abstol = 1e-7), ptopts = (alg = ptalg(prob), reltol = 1e-5, abstol = 1e-5), thread = true, verbose = false, kwargs...)

Compute angular CMB power spectra ``Cₗᴬᴮ`` at angular wavenumbers `ls` from the cosmological problem `prob`.
The requested `modes` are specified as a vector of symbols in the form `:AB`, where `A` and `B` are `T` (temperature), `E` (E-mode polarization) or `ψ` (lensing).
If `unit` is `nothing` the spectra are of dimensionless temperature fluctuations relative to the present photon temperature; while if `unit` is a temperature unit the spectra are of dimensionful temperature fluctuations.
Returns a matrix of ``Cₗ`` if `normalization` is `:Cl`, or ``Dₗ = l(l+1)/2π`` if `normalization` is `:Dl`.

# Precision parameters

- `xs`: Grid of ``(τ-τᵢ)/(τ₀-τᵢ)`` specifying the ``τ``-points that will be sampled in line-of-sight integration.
- `τcut`: Remove all earlier times from the line-of-sight integral sampling time points.
- `kinterp`: Interpolator that decides which ``k``-modes the perturbation ODEs will be solved explicitly for, and then interpolated in-between to a finer grid set by `Δkτ0`.
- `Δkτ0`: Grid spacing to use when integrating over ``k`` to project to ``ℓ``-space.
- `l_limber`: Use Limber approximation for lensing line-of-sight integrals with equal or greater ``ℓ``.
- `bgopts`: Background ODE precision parameters passed to `solvebg`.
- `ptopts`: Perturbation ODE precision parameters passed to `solvept`.

# Examples

```julia
using SymBoltz, Unitful
M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)

ls = 10:10:1000
jl = SphericalBesselCache(ls)
modes = [:TT, :TE, :ψψ, :ψT]
Dls = spectrum_cmb(modes, prob, jl; normalization = :Dl, unit = u"μK")
```

Alternatively, pass a plain vector of integer multipoles instead of a `jl::SphericalBesselCache`:
```julia
Dls = spectrum_cmb(modes, prob, ls; normalization = :Dl, unit = u"μK")
```
to compute ``Φ_l(χ)`` on the fly with [`Φl_recurrence!`](@ref) instead of from a precomputed cache (see that method's docstring for the tradeoffs).
"""
function spectrum_cmb(modes::AbstractVector{<:Symbol}, prob::CosmologyProblem, jl::SphericalBesselCache; kwargs...)
    return spectrum_cmb(modes, prob, collect(jl.l); jl, kwargs...) # the cache already knows which l's it holds
end

"""
    spectrum_cmb(modes::AbstractVector{<:Symbol}, prob::CosmologyProblem, ls::AbstractVector; jl = nothing, kwargs...)

Same as `spectrum_cmb(modes, prob, jl::SphericalBesselCache; kwargs...)`, but without a precomputed spherical
Bessel function cache: ``Φ_l(χ)`` is instead computed on the fly with [`Φl_recurrence!`](@ref) for every
wavenumber and time, which also covers curved universes. `ls` must consist of actual (sorted, ascending) integer
multipoles, since a recursion in `l` can only ever produce values at integer `l` (unlike the cache-based method,
which also supports a coarse, fractional `ls` grid for later interpolation).

This method carries the implementation shared with the cache-based one, which calls it with `ls = jl.l` and the
cache passed as `jl`.
"""
function spectrum_cmb(modes::AbstractVector{<:Symbol}, prob::CosmologyProblem, ls::AbstractArray; jl::Union{SphericalBesselCache, Nothing} = nothing, normalization = :Cl, unit = nothing, kinterp = nothing, Δkτ0 = 2π/4, xs = cosgrid(0.0, 1.0; length=300), τcut = 1e-2, l_limber = 10, bgopts = (alg = bgalg(prob), reltol = 1e-7, abstol = 1e-7), ptopts = (alg = ptalg(prob), reltol = 1e-5, abstol = 1e-5), thread = true, verbose = false, kwargs...)
    ls = collect(ls)
    # Define 1-2-3 indices corresponding for present modes
    iT = 'T' in join(modes) ? 1 : 0
    iE = 'E' in join(modes) ? iT + 1 : 0
    iψ = 'ψ' in join(modes) ? max(iE, iT) + 1 : 0

    sol = solve(prob; bgopts, verbose)
    K = curvature(sol) # spatial curvature constant of the model, passed on to the line-of-sight integration
    τ0 = getsym(sol, prob.M.τ0)(sol)

    # An open universe has no super-curvature modes: its scalar spectrum is bounded below by k² ≥ -K.
    # Including k < √(-K) anyway inflates the lowest multipoles by an order of magnitude.
    kmin = max(1e-2, K < 0 ? √(-K) : 0.0)

    # Automatically determine grid if not provided manually
    if isnothing(kinterp)
        if iψ > 0
            kinterp = ChebyshevInterpolator(kmin, 1e4, 130; f = fk_tanh, f⁻¹ = fk⁻¹_tanh) # higher kmax for lensing; f that stretches acoustic oscillations for k ≲ 2000 with higher sampling density
        else
            kinterp = ChebyshevInterpolator(kmin, 2e3, 60) # lower kmax for T/E-only; sample uniform acoustic oscillations in linear k
        end
    elseif minimum(kinterp) < kmin
        throw(ArgumentError("k-grid starts at k = $(minimum(kinterp)), below the smallest wavenumber √(-K) = $kmin supported by this open universe"))
    end
    ks_fine = lingrid(minimum(kinterp), maximum(kinterp); step=Δkτ0/τ0) # for k-quadrature after LOS integration

    τs = sol.bg.t # by default, use background (thermodynamics) time points for line of sight integration
    τs = τs[τs .≥ τcut]
    if xs isa AbstractArray
        # explicit fractional grid x = (τ-τi)/(τ0-τi) ∈ [0,1]
        xs[begin] == 0 || error("xs begins with $(xs[begin]), but should begin with 0")
        xs[end] == 1 || error("xs ends with $(xs[end]), but should end with 1")
        τs = τs[begin] .+ (τs[end] .- τs[begin]) .* xs
    elseif xs isa Int
        # interpolate xs points from background time grid, preserving its density structure
        τs = LinearInterpolation(τs, 1.0:length(τs)).(range(1.0, length(τs), length = xs))
    end

    # Integrate perturbations to calculate source function on coarse k-grid
    Ss = [S for (S, i) in [(prob.M.k*prob.M.ST, iT), (prob.M.k^2*prob.M.SE, iE), (prob.M.Sψ, iψ)] if i > 0]
    Ss = SVector{length(Ss), eltype(Ss)}(Ss) # turn into SVector
    Ss = source_grid(prob, Ss, τs, ks_fine, kinterp, sol.bg; ptopts, verbose, thread)
    Ss[end, :] .= Ref(zero(eltype(Ss))) # remove any Inf/NaN at last time χ=0; weighted by jₗ(0)=0 anyway

    # Integrate all sources simultaneously without Limber approximation
    Θls = los_integrate(Ss, ls, τs, ks_fine, jl, K; verbose, thread, kwargs...)
    Θls = stack(Θls) # to 3D array
    if iT > 0
        Θls[iT, :, :] ./= ks_fine
    end
    if iE > 0
        Θls[iE, :, :] .*= transpose(@. √((ls+2)*(ls+1)*(ls+0)*(ls-1))) ./ (ks_fine .^ 2)
    end
    if iψ > 0 && l_limber ≤ ls[end]
        Θls[iψ, :, :] .= los_integrate(getindex.(Ss, iψ), ls, τs, ks_fine, jl, K; l_limber, verbose, thread, kwargs...) # overwrite with Limber result
    end

    P0s = spectrum_primordial(ks_fine, sol) # more accurate

    if isnothing(unit)
        factor = 1.0 # keep dimensionless
    elseif dimension(unit) == dimension(u"K")
        factor = uconvert(unit, sol[sol.prob.M.γ.T₀] * u"K") # convert to a temperature unit
    else
        error("Requested unit $unit is not a temperature unit")
    end

    function geti(mode)
        mode == :T && return iT
        mode == :E && return iE
        mode == :ψ && return iψ
        error("Unknown CMB power spectrum mode $mode")
    end

    spectra = zeros(eltype(first(first(Ss)) * P0s[1] * factor^2), length(ls), length(modes)) # Cls or Dls
    for (i, mode) in enumerate(modes)
        mode = String(mode)
        iA = geti(Symbol(mode[firstindex(mode)]))
        iB = geti(Symbol(mode[lastindex(mode)]))
        ΘlAs = @view(Θls[iA, :, :])
        ΘlBs = @view(Θls[iB, :, :])
        spectrum = spectrum_cmb(ΘlAs, ΘlBs, P0s, ls, ks_fine; normalization, thread)
        spectrum *= factor^2 # possibly make dimensionful
        spectra[:, i] .= spectrum
    end

    return spectra
end

"""
    spectrum_cmb(modes::AbstractVector, prob::CosmologyProblem, jl::Union{SphericalBesselCache, AbstractVector}, ls::AbstractVector; kwargs...)

Same, but compute the spectrum properly only for `jl` (a `SphericalBesselCache`, `AbstractInterpolator` over ``l`` or an `AbstractArray` of multipoles)
and then interpolate the results to all `ls`.
"""
function spectrum_cmb(modes::AbstractVector, prob::CosmologyProblem, ls_coarse::AbstractArray, ls::AbstractVector; normalization = :Cl, linterp_normalization = l -> l^5, kwargs...)
    minimum(ls) ≥ minimum(ls_coarse) && maximum(ls) ≤ maximum(ls_coarse) || throw(ArgumentError("l-range $(extrema(ls)) is outside the l-range $(extrema(ls_coarse)) of the coarse l-grid"))
    spectra_coarse = spectrum_cmb(modes, prob, ls_coarse; kwargs...)
    spectra_fine = similar(spectra_coarse, (length(ls), size(spectra_coarse)[2]))
    for imode in eachindex(modes)
        spectra_fine[:, imode] = interpolate(ls_coarse, spectra_coarse[:, imode] .* linterp_normalization.(ls_coarse), ls) ./ linterp_normalization.(ls)
        spectra_fine[:, imode] = normalize_spectrum_cmb(normalization, ls, spectra_fine[:, imode])
    end
    return spectra_fine
end
function spectrum_cmb(modes::AbstractVector, prob::CosmologyProblem, jl::SphericalBesselCache, ls::AbstractVector; normalization = :Cl, linterp_normalization = l -> l^5, kwargs...)
    ls_coarse = jl.l
    minimum(ls) ≥ minimum(ls_coarse) && maximum(ls) ≤ maximum(ls_coarse) || throw(ArgumentError("l-range $(extrema(ls)) is outside the l-range $(extrema(ls_coarse)) of the coarse l-grid"))
    spectra_coarse = spectrum_cmb(modes, prob, jl; kwargs...)
    spectra_fine = similar(spectra_coarse, (length(ls), size(spectra_coarse)[2]))
    for imode in eachindex(modes)
        spectra_fine[:, imode] = interpolate(ls_coarse, spectra_coarse[:, imode] .* linterp_normalization.(ls_coarse), ls) ./ linterp_normalization.(ls)
        spectra_fine[:, imode] = normalize_spectrum_cmb(normalization, ls, spectra_fine[:, imode])
    end
    return spectra_fine
end

function spectrum_cmb(mode::Symbol, args...; kwargs...)
    return spectrum_cmb([mode], args...; kwargs...)[:, begin]
end

normalize_spectrum_cmb(normalization::Nothing, l, Cl) = Cl
normalize_spectrum_cmb(normalization::Function, l, Cl) = normalization.(l) .* Cl
normalize_spectrum_cmb(normalization::Symbol, l, Cl) = normalization == :Dl ? normalize_spectrum_cmb(l -> l*(l+1)/2π, l, Cl) : normalization == :Cl ? Cl : throw(ArgumentError("Normalization symbol is not :Cl or :Dl"))
