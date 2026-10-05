using Bessels: besselj!, sphericalbesselj
using DataInterpolations
using MatterPower
using ForwardDiff

struct SphericalBesselCache{Tl, Tdy <: Union{Matrix{Float64}, Nothing}}
    l::Tl
    y::Matrix{Float64}
    dy::Tdy
    dx::Float64
    invdx::Float64
    x::Vector{Float64}
end

function SphericalBesselCache(ls; xmax = 20*maximum(ls), dx = 2π/15, hermite = true)
    xmin = 0.0
    xs = range(xmin, xmax, length = trunc(Int, (xmax - xmin) / dx)) # fixed length (so endpoints are exact) that gives step as close to dx as possible
    dx = step(xs) # the resulting step, which need not be exactly the requested dx
    invdx = 1.0 / dx
    xs = collect([xs; xs[end]]) # pad with 1 extra duplicate point to avoid bounds check during interpolation
    ys  = jl.(ls, xs') # contiguous in l
    dys = hermite ? jl′.(ls, xs') : nothing
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

# jₗ″ from the spherical Bessel equation x² jₗ″ + 2x jₗ′ + (x² - l(l+1)) jₗ = 0
@inline @fastmath function jl″(l, x, y, dy)
    x == 0 && return l == 0 ? -1/3 : l == 1 ? 0.0 : l == 2 ? 2/15 : 0.0
    invx = 1/x
    return -2dy*invx - (1 - l*(l+1)*invx^2) * y
end

# Hermite interpolation of the cached jₗ′ using jₗ″ at the nodes; more accurate than differentiating the jₗ interpolant
@inline Base.@propagate_inbounds @fastmath function jl′(jl::SphericalBesselCache{Tl, Matrix{Float64}}, il::Int, x) where {Tl}
    w = x * jl.invdx
    i = trunc(Int, w)
    w = w - i
    wm1 = w - 1.0
    l = jl.l[il]
    y₋  = jl.y[il, i+1]
    y₊  = jl.y[il, i+2]
    dy₋ = jl.dy[il, i+1]
    dy₊ = jl.dy[il, i+2]
    ddy₋ = jl″(l, jl.x[i+1], y₋, dy₋)
    ddy₊ = jl″(l, jl.x[i+2], y₊, dy₊)
    return (1+2w)*wm1*wm1 * dy₋ + w*w*(3-2w) * dy₊ + w*wm1 * (wm1 * ddy₋ + w * ddy₊) * jl.dx
end

# Propagate the interpolated derivative through ForwardDiff Duals
@inline Base.@propagate_inbounds function (jl::SphericalBesselCache{Tl, Matrix{Float64}})(il::Int, x::ForwardDiff.Dual{T}) where {Tl, T}
    x₀ = ForwardDiff.value(x)
    return ForwardDiff.Dual{T}(jl(il, x₀), jl′(jl, il, x₀) * ForwardDiff.partials(x))
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

# TODO: line-of-sight integrate Θl using ODE for evolution of Jl?
# TODO: spline sphericalbesselj for each l, from x=0 to x=kmax*(τ0-τini)
# TODO: integrate with ApproxFun? see e.g. https://discourse.julialang.org/t/evaluate-integral-on-many-points-cubature-jl/1723/2
# TODO: RombergEven() works with 513 or 1025 points (do Logging.disable_logging(Logging.Warn) first)
# TODO: gaussian quadrature with weight function? https://juliamath.github.io/QuadGK.jl/stable/weighted-gauss/
# line of sight integration
# TODO: use u = k*χ as integration variable, so oscillations of Bessel functions are the same for every k?
# TODO: define and document symbolic dispatch!
"""
    los_integrate(Ss::AbstractMatrix{T}, ls::AbstractVector, χs::AbstractVector, ks::AbstractVector, jl::SphericalBesselCache, ws::AbstractVector; l_limber = typemax(Int), thread = true, verbose = false) where {T}

For the given `ls` and `ks`, compute the line-of-sight integrals
```math
Iₗ(k) = ∫dχ S(χ,k) jₗ(kχ)
```
over the source function values `Ss` against the spherical Bessel functions ``jₗ(x)`` cached in `jl`, using the quadrature weights `ws` for the descending conformal distances ``χ = τ₀ - τ``.
The element `Ss[i,j]` holds the source function value ``S(χᵢ, kⱼ)``.
The Limber approximation
```math
Iₗ ≈ √(π/(2l+1)) S((l+1/2)/k, k)
```
is used for `l ≥ l_limber`.
Contributions where ``|jₗ(x)|`` is below `jltol` (at small ``x ≪ l``) are skipped.
"""
function los_integrate(Ss::AbstractMatrix{T}, ls::AbstractVector, χs::AbstractVector, ks::AbstractVector, jl::SphericalBesselCache, ws::AbstractVector; l_limber = typemax(Int), jltol = 1e-20, thread = true, verbose = false) where {T}
    @assert size(Ss, 1) == length(χs) == length(ws) "size(Ss, 1) = $(size(Ss, 1)), length(χs) = $(length(χs)) and length(ws) = $(length(ws)) differ"
    @assert size(Ss, 2) == length(ks) "size(Ss, 2) = $(size(Ss, 2)) and length(ks) = $(length(ks)) differ"
    @assert collect(ls) == collect(jl.l) "ls must match the l-values stored in the Bessel cache"
    @assert jl.x[begin] ≤ 0 "jl.x[begin] < 0"
    @assert jl.x[end] ≥ ks[end]*χs[begin] "jl.x[end] < kmax*χmax"
    @assert issorted(χs; rev = true) "χs must be sorted in descending order"
    @assert issorted(ks) "ks must be sorted in ascending order"
    @assert issorted(ls) "ls must be sorted in ascending order" # necessary for Limber indexing logic
    error_if_nonfinite(Ss)

    nχ = length(χs)

    nl = length(ls)
    Is = similar(Ss, length(ks), nl)
    il_limber = searchsortedfirst(ls, l_limber) # First il index with l ≥ l_limber (=nl+1 when l_limber = typemax, i.e. no Limber modes)

    # Skip negligible jₗ(x) at small x (avoids wasted work and very slow subnormal arithmetic)
    ilend = [something(findlast(il -> abs(jl.y[il, ix]) ≥ jltol, 1:il_limber-1), 0) for ix in eachindex(jl.x)] # last il with non-negligible jₗ for each x-node
    # TODO: `accumulate` with `max` to guarantee that ilend is monotonically increasing/decreasing?

    verbose && l_limber < typemax(Int) && println("Using Limber approximation for l ≥ $l_limber")

    # Loop order k → χ → l to get SIMD on the innermost l-loop
    @fastmath @inbounds @tasks for ik in eachindex(ks)
        @set scheduler = thread ? :dynamic : :serial
        @local tmp = zeros(T, nl) # l-contiguous storage for integrals (to help SIMD over l)
        k = ks[ik]
        verbose && print("\rLOS integrating k-mode $ik / $(length(ks))")

        # Full line-of-sight integrals for l < l_limber
        fill!(tmp, zero(T))
        @inbounds for iχ in eachindex(χs)
            kχ = k * χs[iχ]
            Sw = ws[iχ] * Ss[iχ, ik]
            ix = trunc(Int, ForwardDiff.value(kχ) * jl.invdx) + 2 # rightmost interpolation node corresponding to kχ
            @inbounds @simd for il in 1:ilend[ix]
                tmp[il] += Sw * jl(il, kχ)
            end
        end

        # Limber approximation for l ≥ l_limber
        @inbounds for il in il_limber:nl
            l = ls[il]
            χ = (l + 1/2) / k
            if χ ≤ χs[1] # otherwise source is zero before recombination
                i₋ = searchsortedfirst(χs, χ; rev = true)
                χ₋ = χs[i₋]
                S₋ = Ss[i₋, ik]
                if i₋ == 1
                    S = S₋
                else
                    i₊ = i₋ - 1 # χs is descending, so χ₋ < χ < χ₊
                    χ₊ = χs[i₊]
                    S₊ = Ss[i₊, ik]
                    Δχ = χ₊ - χ₋
                    S′₋ = i₋ ≤ nχ-1 ? (Ss[i₋+1, ik] - S₊) / (χs[i₋+1] - χ₊) : (S₊ - S₋) / Δχ
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

@doc raw"""
    spectrum_cmb(ΘlAs::AbstractMatrix, ΘlBs::AbstractMatrix, P0s::AbstractVector, ls::AbstractVector, ks::AbstractVector, ws::AbstractVector; normalization = :Cl, thread = true)

Compute the angular power spectrum
```math
Cₗᴬᴮ = (2/π) ∫\mathrm{d}k \, k² P₀(k) Θₗᴬ(τ₀,k) Θₗᴮ(τ₀,k)
```
for the given `ls`, using the quadrature weights `ws` for the wavenumbers `ks`.
If `normalization == :Dl`, compute ``Dₗ = Cₗ l (l+1) / 2π`` instead.
"""
function spectrum_cmb(ΘlAs::AbstractMatrix, ΘlBs::AbstractMatrix, P0s::AbstractVector, ls::AbstractVector, ks::AbstractVector, ws::AbstractVector; normalization = :Cl, thread = true)
    size(ΘlAs) == size(ΘlBs) || error("ΘlAs and ΘlBs have different sizes")
    eltype(ΘlAs) == eltype(ΘlBs) || error("ΘlAs and ΘlBs have different types")

    Cls = similar(ΘlAs, length(ls))

    @tasks for il in eachindex(ls)
        # TODO: skip kτ0 ≲ l?
        @set scheduler = thread ? :dynamic : :static
        @local dCl_dks = zeros(eltype(ΘlAs), length(ks)) # local task workspace
        ΘlA = @view ΘlAs[:, il]
        ΘlB = @view ΘlBs[:, il]
        dCl_dks .= 2/π .* ks .^ 2 .* P0s .* ΘlA .* ΘlB
        Cls[il] = sum(ws[ik] * dCl_dks[ik] for ik in eachindex(ws)) # integrate over k
    end

    return normalize_spectrum_cmb(normalization, ls, Cls)
end

fk_tanh(k, k0=2000.0) = tanh(k/k0)
fk⁻¹_tanh(k, k0=2000.0) = k0*atanh(k)

"""
    default_τquad(τi, τrec, τ0; N = 600, A = 20, w = τrec/3)

Create a trapezoidal quadrature rule with `N` nodes for (line-of-sight) integration over ``τ ∈ [τᵢ, τ₀]``, mapped to ``[-1, 1]``.
The nodes have density ``d(τ) = 1 + A \\mathrm{sech}^2((τ-τ_\\mathrm{rec})/w)``: a uniform base plus a bump of height `A` and width `w` around recombination at `τrec`.
It applies the trapezoidal rule in the cumulative density ``u(τ) = ∫d(τ) dτ`` sampled uniformly,
after the change of variables ``∫f(τ) dτ = ∫f(τ(u))/d(τ) du``.
"""
function default_τquad(τi, τrec, τ0; N = 600, A = 20, w = τrec/3)
    d(τ) = 1 + A * sech((τ-τrec)/w)^2
    u(τ) = τ + A * w * tanh((τ-τrec)/w)
    τf = range(τi, τ0, length = 10_000)
    us = range(u(τi), u(τ0), length = N)
    τs = CubicHermiteSpline(1 ./ d.(τf), collect(τf), u.(τf))(us) # invert u(τ) with exact derivative dτ/du = 1/d
    τs[begin], τs[end] = τi, τ0
    ws = step(us) ./ d.(τs)
    ws[begin] /= 2
    ws[end] /= 2
    x = clamp.(2 .* (τs .- τi) ./ (τ0 - τi) .- 1, -1, 1)
    return Quadrature(x, 2 .* ws ./ sum(ws); name = :Recombination)
end

"""
    spectrum_cmb(modes::AbstractVector{<:Symbol}, prob::CosmologyProblem, jl::SphericalBesselCache; normalization = :Cl, kinterp = nothing, τquad = nothing, kquad = nothing, l_limber = 11, bgalg = default_bgalg(prob), bgreltol = 1e-7, bgabstol = 1e-7, bgopts = (), ptalg = default_ptalg(prob), ptreltol = 1e-5, ptabstol = 1e-5, ptopts = (), thread = true, verbose = false, kwargs...)

Compute angular CMB power spectra ``Cₗᴬᴮ`` at angular wavenumbers `ls` from the cosmological problem `prob`.
The requested `modes` are specified as a vector of symbols in the form `:AB`, where `A` and `B` are `T` (temperature), `E` (E-mode polarization) or `ψ` (lensing).
The spectra are of dimensionless temperature fluctuations relative to the present photon temperature ``T_{γ0}``; multiply by ``T_{γ0}^2`` to get dimensionful spectra.
Returns a matrix of ``Cₗ`` if `normalization` is `:Cl`, or ``Dₗ = l(l+1)/2π`` if `normalization` is `:Dl`.

# Precision parameters

- `τquad`: Quadrature rule for the line-of-sight integral over ``τ``; its nodes are mapped linearly to ``[τᵢ, τ₀]``. Defaults to [`default_τquad`](@ref).
- `kquad`: Quadrature rule for line-of-sight integration and the integral over ``k``; its nodes are mapped linearly to the ``k``-range of `kinterp`.
- `kinterp`: Interpolator that decides which ``k``-modes the perturbation ODEs will be solved explicitly for, and then interpolated in-between to the nodes of `kquad`.
- `l_limber`: Use Limber approximation for lensing line-of-sight integrals with equal or greater ``ℓ``.
- `bgalg`/`ptalg`, `bgreltol`/`ptreltol`, `bgabstol`/`ptabstol`: ODE algorithms and tolerances for the background/perturbation stages.
- `bgopts`/`ptopts`: extra options for the background/perturbation ODE solves.

# Examples

```julia
using SymBoltz
M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)

ls = 10:10:1000
jl = SphericalBesselCache(ls)
modes = [:TT, :TE, :ψψ, :ψT]
Dls = spectrum_cmb(modes, prob, jl; normalization = :Dl)
```
"""
function spectrum_cmb(modes::AbstractVector{<:Symbol}, prob::CosmologyProblem, jl::SphericalBesselCache; normalization = :Cl, kinterp = nothing, τquad = nothing, kquad = nothing, l_limber = 11, bgalg = default_bgalg(prob), bgreltol = 1e-7, bgabstol = 1e-7, bgopts = (), ptalg = default_ptalg(prob), ptreltol = 1e-5, ptabstol = 1e-5, ptopts = (), thread = true, verbose = false, kwargs...)
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
    sol = solve(prob; bgalg, bgreltol, bgabstol, bgopts, verbose)
    τbg = sol[prob.M.τ] # conformal times at background time points (also if τ is not the independent variable)
    τi, τ0 = τbg[begin], τbg[end]
    v = hasproperty(prob.M, :v) ? prob.M.v : prob.M.b.v # visibility function in unstructured or structured models
    τrec = τbg[argmax(sol[v])] # recombination at peak of the visibility function
    if isnothing(τquad)
        τquad = default_τquad(ForwardDiff.value.((τi, τrec, τ0))...) # shape of quadrature rule is parameter-independent
    end

    kmin, kmax = extrema(kinterp)
    if isnothing(kquad)
        s = 1e3 # uniform k-spacing for k ≲ s where T/E oscillate uniformly, but logarithmic after damping in the lensing tail k ≳ s
        Δk = 0.5 * π / ForwardDiff.value(τ0) # ≈ 2 points per period π/χ of the integrand ∝ jₗ(kχ)²; drop derivatives because k-limits are parameter-independent
        kgrid = asinhgrid(kmin, kmax, s; step = Δk/s) # uniform spacing Δk for k ≲ s; logarithmic spacing for k ≳ s
        kquad = TrapezoidalQuadrature(kgrid)
    end
    ks_fine = nodes(kquad, kmin, kmax) # for k-quadrature after LOS integration
    kws = weights(kquad, kmin, kmax)

    tbg = timeseries(sol) # independent variable at background time points
    tmin, tmax = extrema(tbg)
    ts = LinearInterpolation(tbg, τbg; extrapolation = ExtrapolationType.Extension)(nodes(τquad, τi, τ0)) # map τ-nodes to the independent variable (e.g. τ or ln(a)); extrapolate to avoid out-of-bounds errors from rounding at the boundaries
    ts = clamp.(ts, tmin, tmax) # avoid rounding errors at boundaries if rescaling pushes times outside the background timespan
    τs = sol(prob.M.τ, ts) # conformal times at the sampled points for line-of-sight integration
    χs = τ0 .- τs # conformal distances (descending)
    τws = weights(τquad, τi, τ0)

    # Integrate perturbations to calculate source function on coarse k-grid
    Ss = [S for (S, i) in [(prob.M.k*prob.M.ST, iT), (prob.M.k^2*prob.M.SE, iE), (prob.M.Sψ, iψ)] if i > 0]
    Ss = SVector{length(Ss), eltype(Ss)}(Ss) # turn into SVector
    Ss = source_grid(prob, Ss, ts, ks_fine, kinterp, sol.bg; ptalg, ptreltol, ptabstol, ptopts, verbose, thread)
    if iψ > 0
        # apply lensing kernel for a thin last scattering surface at the peak of the visibility function # TODO: use more accurate Hermite interpolation?
        Ws = [τ ≥ τrec ? (τ-τrec)/(τ0-τrec)/(τ0-τ) : zero(τ) for τ in τs]
        for iτ in eachindex(τs), ik in eachindex(ks_fine)
            Ss[iτ, ik] = Base.setindex(Ss[iτ, ik], Ss[iτ, ik][iψ] * Ws[iτ], iψ)
        end
    end
    if χs[end] == 0 # remove any Inf/NaN at last time χ=0; weighted by jₗ(0)=0 anyway
        Ss[end, :] .= Ref(zero(eltype(Ss)))
    end

    # Integrate all sources simultaneously without Limber approximation
    Θls = los_integrate(Ss, ls, χs, ks_fine, jl, τws; verbose, thread, kwargs...)
    Θls = stack(Θls) # to 3D array
    if iT > 0
        Θls[iT, :, :] ./= ks_fine
    end
    if iE > 0
        Θls[iE, :, :] .*= transpose(@. √((ls+2)*(ls+1)*(ls+0)*(ls-1))) ./ (ks_fine .^ 2)
    end
    if iψ > 0 && l_limber ≤ ls[end]
        Θls[iψ, :, :] .= los_integrate(getindex.(Ss, iψ), ls, χs, ks_fine, jl, τws; l_limber, verbose, thread, kwargs...) # overwrite with Limber result
    end

    P0s = spectrum_primordial(ks_fine, sol) # more accurate

    function geti(mode)
        mode == :T && return iT
        mode == :E && return iE
        mode == :ψ && return iψ
        error("Unknown CMB power spectrum mode $mode")
    end

    spectra = zeros(eltype(first(first(Ss)) * P0s[1]), length(ls), length(modes)) # Cls or Dls
    for (i, mode) in enumerate(modes)
        mode = String(mode)
        iA = geti(Symbol(mode[firstindex(mode)]))
        iB = geti(Symbol(mode[lastindex(mode)]))
        ΘlAs = @view(Θls[iA, :, :])
        ΘlBs = @view(Θls[iB, :, :])
        spectra[:, i] .= spectrum_cmb(ΘlAs, ΘlBs, P0s, ls, ks_fine, kws; normalization, thread)
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
