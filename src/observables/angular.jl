using Integrals
using Bessels: besselj!, sphericalbesselj
using DataInterpolations
using MatterPower
using ForwardDiff
using ForwardDiffChainRules
import ChainRulesCore
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqVerner
using SpecialFunctions: loggamma, logabsgamma

struct SphericalBesselCache{Tl, Tdy <: Union{Matrix{Float64}, Nothing}}
    l::Tl
    y::Matrix{Float64}
    dy::Tdy
    dx::Float64
    invdx::Float64
    x::Vector{Float64}
end

"""
    SphericalBesselCache(ls; xmax = 20*maximum(ls), dx = 2π/15, hermite = true, method = :interp, odekwargs = (;))

Build a cache of the spherical Bessel functions ``jₗ(x)`` for `l ∈ ls` and `x ∈ [0, xmax]` on a uniform grid with spacing close to `dx`, for fast repeated interpolated evaluation (linear if `hermite = false`, cubic Hermite otherwise).

If `method = :interp` (the default), the grid is filled by direct evaluation of `jₗ(x)` (using Bessels.jl).
If `method = :ode`, each row is instead computed by solving the defining ODE of ``jₗ(x)``,
```math
jₗ'' + (2/x) jₗ' + (1 - l(l+1)/x²) jₗ = 0 ,
```
forward in `x` with `saveat` set to the cache's `x`-grid, using the solver and tolerances in `odekwargs` (default `OrdinaryDiffEqTsit5.Tsit5()` with `reltol = 1e-8, abstol = 1e-8`).
See [`spherical_bessel_ode_grid`](@ref) for details on how each `l`'s ODE is seeded and started away from the singular point `x=0`.
"""
function SphericalBesselCache(ls; xmax = 20*maximum(ls), dx = 2π/15, hermite = true, method = :interp, odekwargs = (;))
    xmin = 0.0
    xs = range(xmin, xmax, length = trunc(Int, (xmax - xmin) / dx)) # fixed length (so endpoints are exact) that gives step as close to dx as possible
    invdx = 1.0 / step(xs) # using the resulting step, which need not be exactly dx
    xs = collect([xs; xs[end]]) # pad with 1 extra duplicate point to avoid bounds check during interpolation
    if method == :interp
        ys  = jl.(ls, xs') # contiguous in l
        dys = hermite ? jl′.(ls, xs') : nothing
    elseif method == :ode
        ys, dys = spherical_bessel_ode_grid(ls, xs; hermite, odekwargs...)
    else
        throw(ArgumentError("Unknown SphericalBesselCache method $method (expected :interp or :ode)"))
    end
    return SphericalBesselCache{typeof(ls), typeof(dys)}(ls, ys, dys, dx, invdx, xs)
end

# ---------------------------------------------------------------------------
# ODE-based construction (method = :ode above)
# ---------------------------------------------------------------------------
function spherical_bessel_rhs!(du, u, l, x)
    du[1] = u[2]
    du[2] = -(1 - l*(l+1)/x^2)*u[1] - 2/x*u[2]
end

# WKB tunneling exponent ∫ₓ₀^ν √(ν²/x²-1) dx from x0 to the turning point ν=l+1/2 (closed form).
# Governs how much any admixture of the "other" solution gets amplified/damped between x0 and ν.
function _tunneling_exponent(l::Real, x0::Real)
    ν = l + 0.5
    x0 >= ν && return 0.0
    s = sqrt(max(ν^2 - x0^2, 0.0))
    return ν*log((ν+s)/x0) - s
end

# Starting point for the regular solution J (∝ jₗ): moderately deep in the classically forbidden
# region (bounded WKB exponent from the turning point), not deep enough to lose all precision, but
# deep enough that an arbitrary seed there "washes out" to the correct (unnormalized) shape once
# integrated forward past the turning point (see spherical_bessel_ode_grid).
_x0(l::Real) = max(0.01, l - 5*sqrt(l))

# Starting point for the independent second solution N: chosen close to x0 (a small, fixed WKB
# exponent gap dE away) so that forward-integrating N from xa to x0 also stays accurate — N is huge
# and shrinking there, so integrating it too far risks the same precision loss in reverse.
function _xa(l::Real, x0::Real; dE::Real = 3.0)
    target = _tunneling_exponent(l, x0) + dE
    lo, hi = 1e-8, x0*0.999999
    for _ in 1:200
        mid = sqrt(lo*hi)
        if _tunneling_exponent(l, mid) > target
            lo = mid
        else
            hi = mid
        end
    end
    return sqrt(lo*hi)
end

# Small-x series of the independent second solution N(x) = x^(-(l+1)) * [1 + O(x²)], generalized to
# real l via the Gamma function. Unlike jₗ's own small-x series, N is not exponentially suppressed
# (it is the dominant solution there), so this converges in plain Float64 with no cancellation issues
# — no arbitrary-precision arithmetic needed even though N(xa) itself can be astronomically large.
function _N_series(l::Real, x::Real; maxterms::Integer = 2000)
    loggam, siggam = logabsgamma(-l + 0.5)
    log_a0 = log(sqrt(π)/2) + (l+1)*log(2/x) - loggam
    a0 = -siggam * exp(log_a0)
    term = a0
    y, dy = term, term*(-(l+1))/x
    x2 = x^2
    for k in 1:maxterms
        ratio = -x2/4 / (k * (-l + 0.5 + (k-1)))
        term *= ratio
        dy += term * (-(l+1) + 2k) / x
        y += term
        abs(term) < 1e-15*abs(y) && k > 3 && break
    end
    return y, dy
end

# Exact, analytic Wronskian constant x²(J·N′ - J′·N) for the (J, N) pair above, computed from their
# leading small-x coefficients (valid since the Wronskian of two ODE solutions in Sturm-Liouville
# form x²(J N′ - J′ N) is x-independent, so it can be evaluated in the x→0 limit using only the
# leading terms of either series).
function _wronskian(l::Real)
    log_jlead = log(sqrt(π)/2) - l*log(2) - loggamma(l + 1.5)
    loggam, siggam = logabsgamma(-l + 0.5)
    log_nlead = log(sqrt(π)/2) + (l+1)*log(2) - loggam
    logW = log(2l+1) + log_jlead + log_nlead
    return siggam * exp(logW)
end

"""
    spherical_bessel_ode_grid(ls, xs::AbstractVector; hermite = true, odealg = Tsit5(), odereltol = 1e-8, odeabstol = 1e-8)

For each `l` in `ls`, compute `jₗ(x)` (and, if `hermite`, `jₗ′(x)`) on the grid `xs` by solving the
defining ODE of the spherical Bessel function forward in `x`, without relying on Bessels.jl anywhere
(needed since `l` need not be an integer, e.g. for Chebyshev interpolation over `l`):

Forward-integrating from very close to the singular point `x=0` is numerically unstable for large `l`:
any infinitesimal error acquired deep inside the classically forbidden region `x≪l` is amplified
exponentially by the time it reaches the oscillatory region past the turning point `x≈l`. Instead, the
regular solution `J` (proportional to `jₗ`) is seeded with an *arbitrary* value at `x0` from
[`_x0`](@ref) (moderately deep, bounded WKB exponent) and integrated forward: this "washes out" to the
correct shape up to an unknown overall scale `B`. That scale is then fixed via the Wronskian: an
independent second solution `N` is seeded *accurately* from its own small-`x` series ([`_N_series`](@ref),
which — unlike `jₗ`'s — has no cancellation problem and needs no arbitrary precision) at a nearby point
`xa` from [`_xa`](@ref), integrated forward to `x0`, and combined with `J` there through the Wronskian
invariant (analytic, [`_wronskian`](@ref)) to solve for `B`. Grid points below `x0` are negligible and
filled with `0` (or `1` for `l=0` at `x=0`).
"""
function spherical_bessel_ode_grid(ls, xs::AbstractVector; hermite::Bool = true, odealg = Tsit5(), odereltol::Real = 1e-8, odeabstol::Real = 1e-8)
    nl, nx = length(ls), length(xs)
    ys = Matrix{Float64}(undef, nl, nx)
    dys = hermite ? Matrix{Float64}(undef, nl, nx) : nothing

    for (il, l) in enumerate(ls)
        x0 = _x0(l)
        x0 < xs[end] || throw(ArgumentError("Starting point x0=$x0 for l=$l lies beyond the grid's maximum x=$(xs[end]); increase xmax"))
        xa = _xa(l, x0)

        probJ = ODEProblem(spherical_bessel_rhs!, [1.0, 0.0], (x0, xs[end]), l)

        y0N, dy0N = _N_series(l, xa)
        probN = ODEProblem(spherical_bessel_rhs!, [y0N, dy0N], (xa, x0), l)
        solN = solve(probN, odealg; reltol = odereltol, abstol = odeabstol)

        i0 = searchsortedfirst(xs, x0)
        saveat = @view xs[i0:end]
        solJ = solve(probJ, odealg; saveat, save_start = false, reltol = odereltol, abstol = odeabstol)
        length(solJ.u) == length(saveat) || error("ODE solve for l=$l did not save at every requested x (retcode $(solJ.retcode))")

        Jx0, dJx0 = 1.0, 0.0 # J's own seed, at t=x0
        Nx0, dNx0 = solN.u[end]
        W_raw = x0^2 * (Jx0*dNx0 - dJx0*Nx0)
        B = W_raw / _wronskian(l)

        for (j, i) in enumerate(i0:nx)
            ys[il, i] = solJ.u[j][1] / B
        end
        j0 = l == 0 ? 1.0 : 0.0 # jₗ(0); the region below x0 is negligible
        @views ys[il, 1:i0-1] .= j0

        if hermite
            for (j, i) in enumerate(i0:nx)
                dys[il, i] = solJ.u[j][2] / B
            end
            @views dys[il, 1:i0-1] .= 0.0
        end
    end

    return ys, dys
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
jl′(l, x) = l/(2l+1)*jl(l-1,x) - (l+1)/(2l+1)*jl(l+1,x) # for l ≥ 1, analytical relation

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
    los_integrate(Ss::AbstractMatrix{T}, ls::AbstractVector, τs::AbstractVector, ks::AbstractVector, jl::SphericalBesselCache; l_limber = typemax(Int), integrator = TrapezoidalRule(), thread = true, verbose = false) where {T}

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
function los_integrate(Ss::AbstractMatrix{T}, ls::AbstractVector, τs::AbstractVector, ks::AbstractVector, jl::SphericalBesselCache; l_limber = typemax(Int), integrator = TrapezoidalRule(), thread = true, verbose = false) where {T}
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

# TODO: integrate splines instead of trapz! https://discourse.julialang.org/t/how-to-speed-up-the-numerical-integration-with-interpolation/96223/5
@doc raw"""
    spectrum_cmb(ΘlAs::AbstractMatrix, ΘlBs::AbstractMatrix, P0s::AbstractVector, ls::AbstractVector, ks::AbstractVector; integrator = TrapezoidalRule(), normalization = :Cl, thread = true)

Compute the angular power spectrum
```math
Cₗᴬᴮ = (2/π) ∫\mathrm{d}k \, k² P₀(k) Θₗᴬ(k,τ₀) Θₗᴮ(k,τ₀)
```
for the given `ls`.
If `normalization == :Dl`, compute ``Dₗ = Cₗ l (l+1) / 2π`` instead.
"""
function spectrum_cmb(ΘlAs::AbstractMatrix, ΘlBs::AbstractMatrix, P0s::AbstractVector, ls::AbstractVector, ks::AbstractVector; integrator = TrapezoidalRule(), normalization = :Cl, thread = true)
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
    spectrum_cmb(modes::AbstractVector{<:Symbol}, prob::CosmologyProblem, jl::SphericalBesselCache; normalization = :Cl, unit = nothing, kinterp = nothing, Δkτ0 = 2π/4, xs = cosgrid(0.0, 1.0; length=300), τcut = 1e-2, l_limber = 10, integrator = TrapezoidalRule(), bgopts = (alg = bgalg(prob), reltol = 1e-7, abstol = 1e-7), ptopts = (alg = ptalg(prob), reltol = 1e-5, abstol = 1e-5), thread = true, verbose = false, kwargs...)

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
"""
function spectrum_cmb(modes::AbstractVector{<:Symbol}, prob::CosmologyProblem, jl::SphericalBesselCache; normalization = :Cl, unit = nothing, kinterp = nothing, Δkτ0 = 2π/4, xs = cosgrid(0.0, 1.0; length=300), τcut = 1e-2, l_limber = 10, integrator = TrapezoidalRule(), bgopts = (alg = bgalg(prob), reltol = 1e-7, abstol = 1e-7), ptopts = (alg = ptalg(prob), reltol = 1e-5, abstol = 1e-5), thread = true, verbose = false, kwargs...)
    # Define 1-2-3 indices corresponding for present modes
    iT = 'T' in join(modes) ? 1 : 0
    iE = 'E' in join(modes) ? iT + 1 : 0
    iψ = 'ψ' in join(modes) ? max(iE, iT) + 1 : 0

    # Automatically determine grid if not provided manually
    if isnothing(kinterp)
        if iψ > 0
            kinterp = ChebyshevInterpolator(1e-2, 1e4, 130; f = fk_tanh, f⁻¹ = fk⁻¹_tanh) # higher kmax for lensing; f that stretches acoustic oscillations for k ≲ 2000 with higher sampling density
        else
            kinterp = ChebyshevInterpolator(1e-2, 2e3, 60) # lower kmax for T/E-only; sample uniform acoustic oscillations in linear k
        end
    end

    ls = collect(jl.l)
    sol = solve(prob; bgopts, verbose)
    τ0 = getsym(sol, prob.M.τ0)(sol)
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
    Θls = los_integrate(Ss, ls, τs, ks_fine, jl; integrator, verbose, thread, kwargs...)
    Θls = stack(Θls) # to 3D array
    if iT > 0
        Θls[iT, :, :] ./= ks_fine
    end
    if iE > 0
        Θls[iE, :, :] .*= transpose(@. √((ls+2)*(ls+1)*(ls+0)*(ls-1))) ./ (ks_fine .^ 2)
    end
    if iψ > 0 && l_limber ≤ ls[end]
        Θls[iψ, :, :] .= los_integrate(getindex.(Ss, iψ), ls, τs, ks_fine, jl; l_limber, integrator, verbose, thread, kwargs...) # overwrite with Limber result
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
        spectrum = spectrum_cmb(ΘlAs, ΘlBs, P0s, ls, ks_fine; integrator, normalization, thread)
        spectrum *= factor^2 # possibly make dimensionful
        spectra[:, i] .= spectrum
    end

    return spectra
end

"""
    spectrum_cmb(modes::AbstractVector, prob::CosmologyProblem, jl::SphericalBesselCache, ls::AbstractVector; kwargs...)

Same, but compute the spectrum properly only for `jl.l` and then interpolate the results to all `ls`.
"""
function spectrum_cmb(modes::AbstractVector, prob::CosmologyProblem, jl::SphericalBesselCache, ls::AbstractVector; normalization = :Cl, linterp_normalization = l -> l^5, kwargs...)
    minimum(ls) ≥ minimum(jl.l) && maximum(ls) ≤ maximum(jl.l) || throw(ArgumentError("l-range $(extrema(ls)) is outside the l-range $(extrema(jl.l)) of the spherical Bessel function"))
    spectra_coarse = spectrum_cmb(modes, prob, jl; kwargs...)
    spectra_fine = similar(spectra_coarse, (length(ls), size(spectra_coarse)[2]))
    for imode in eachindex(modes)
        spectra_fine[:, imode] = interpolate(jl.l, spectra_coarse[:, imode] .* linterp_normalization.(jl.l), ls) ./ linterp_normalization.(ls) # interpolate l⁵*Cₗ (by default) for smoothness
        spectra_fine[:, imode] = normalize_spectrum_cmb(normalization, ls, spectra_fine[:, imode]) # normalize AFTER interpolation
    end
    return spectra_fine
end

function spectrum_cmb(mode::Symbol, args...; kwargs...)
    return spectrum_cmb([mode], args...; kwargs...)[:, begin]
end

normalize_spectrum_cmb(normalization::Nothing, l, Cl) = Cl
normalize_spectrum_cmb(normalization::Function, l, Cl) = normalization.(l) .* Cl
normalize_spectrum_cmb(normalization::Symbol, l, Cl) = normalization == :Dl ? normalize_spectrum_cmb(l -> l*(l+1)/2π, l, Cl) : normalization == :Cl ? Cl : throw(ArgumentError("Normalization symbol is not :Cl or :Dl"))
