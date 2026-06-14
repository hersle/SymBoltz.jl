using Integrals
using Bessels: besselj!, sphericalbesselj
using DataInterpolations
using MatterPower
using ForwardDiff
using ForwardDiffChainRules
import ChainRulesCore

struct SphericalBesselCache{Tdy <: Union{Matrix{Float64}, Nothing}}
    l::Vector{Int}
    i::Vector{Int}
    y::Matrix{Float64}
    dy::Tdy
    dx::Float64
    invdx::Float64
    x::Vector{Float64}
end

function SphericalBesselCache(ls::AbstractVector; xmax = 20*ls[end], dx = 2π/15, hermite = true)
    xmin = 0.0
    xs = range(xmin, xmax, length = trunc(Int, (xmax - xmin) / dx)) # fixed length (so endpoints are exact) that gives step as close to dx as possible
    invdx = 1.0 / step(xs) # using the resulting step, which need not be exactly dx
    xs = collect([xs; xs[end]]) # pad with 1 extra duplicate point to avoid bounds check during interpolation

    is = zeros(Int, maximum(ls))
    for (i, l) in enumerate(ls)
        is[l] = i
    end

    ys = jl.(ls', xs)
    dys = hermite ? jl′.(ls', xs) : nothing

    return SphericalBesselCache{typeof(dys)}(ls, is, ys, dys, dx, invdx, xs)
end

# TODO: define chain rule like in https://github.com/JuliaDiff/ForwardDiff.jl/blob/master/src/dual.jl?
Base.@propagate_inbounds @fastmath function (jl::SphericalBesselCache{Nothing})(l, x)
    il = jl.i[l]
    w = x * jl.invdx # 0-based float index (assume x0 = 0)
    i = trunc(Int, w) # 0-based integer index of left interval point; faster than searchsortedfirst(jl.x, x)
    w = w - i # remainder ∈ [0, 1]
    y₋ = jl.y[i+1, il] # +1 for 1-based indexing
    y₊ = jl.y[i+2, il]
    return muladd(w, y₊ - y₋, y₋) # i.e. y₋ + (y₊ - y₋) * (x - x₋) * jl.invdx
end

Base.@propagate_inbounds @fastmath function (jl::SphericalBesselCache{Matrix{Float64}})(l, x)
    il = jl.i[l]
    w = x * jl.invdx
    i = trunc(Int, w)
    w = w - i
    wm1 = w - 1.0
    y₋ = jl.y[i+1, il]
    y₊ = jl.y[i+2, il]
    dy₋ = jl.dy[i+1, il]
    dy₊ = jl.dy[i+2, il]
    return (1+2w)*wm1*wm1 * y₋ + w*w*(3-2w) * y₊ + w*wm1 * (wm1 * dy₋ + w * dy₊) * jl.dx # https://en.wikipedia.org/wiki/Cubic_Hermite_spline
end

function Base.show(io::IO, jl::SphericalBesselCache{T}) where {T}
    method = T == Nothing ? "linear" : "Hermite"
    print(io, "jₗ(x) $method interpolation cache ")
    print(io, "for $(jl.l[begin]) ≤ l ≤ $(jl.l[end]) and ")
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
    # Julia is column-major; make sure innermost loop indices appear first in slice expressions (https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-column-major)
    @assert size(Ss, 1) == length(τs) "size(Ss, 1) = $(size(Ss, 1)) and length(τs) = $(length(τs)) differ"
    @assert size(Ss, 2) == length(ks) "size(Ss, 2) = $(size(Ss, 2)) and length(ks) = $(length(ks)) differ"
    @assert jl.x[begin] ≤ 0 "jl.x[begin] < 0"
    @assert jl.x[end] ≥ ks[end]*τs[end] "jl.x[end] < kmax*τmax"
    @assert all(all(isfinite.(S)) for S in Ss) "Ss contain NaN or Inf"
    @assert τs[2] > τs[1] "τs must be sorted in ascending order"
    @assert ks[2] > ks[1] "ks must be sorted in ascending order"
    τs = collect(τs) # force array to avoid floating point errors with ranges in following χs due to (e.g. tiny negative χ)
    τ0 = τs[end]
    χs = τ0 .- τs
    halfdτs = 0.5 .* (τs[begin+1:end] .- τs[begin:end-1]) # precompute before loops
    Is = similar(Ss, length(ks), length(ls))

    verbose && l_limber < typemax(Int) && println("Using Limber approximation for l ≥ $l_limber")

    # TODO: skip and set jl to zero if l ≳ kτ0 or another cutoff?
    @fastmath @inbounds @tasks for il in eachindex(ls) # parallellize independent loop iterations
        @set scheduler = thread ? :dynamic : :serial
        l = ls[il]
        verbose && print("\rLOS integrating with l = $l")
        for ik in eachindex(ks)
            k = ks[ik]
            I = zero(T)
            if l ≥ l_limber
                χ = (l+1/2) / k
                if χ ≤ χs[1] # otherwise χ > χini > χrec and source function is definitely zero
                    # cubic Hermite interpolation between two closest points
                    i₋ = searchsortedfirst(τs, τ0 - χ) # highest index
                    χ₋ = χs[i₋]
                    S₋ = Ss[i₋, ik]
                    if i₋ == 1
                        S = S₋
                    else
                        i₊ = i₋ - 1 # lowest index; χs is sorted in descending order, so χ₋ < χ < χ₊
                        χ₊ = χs[i₊]
                        S₊ = Ss[i₊, ik]
                        Δχ = χ₊ - χ₋
                        S′₋ = i₋ ≤ length(τs)-1 ? (Ss[i₋+1, ik] - S₊) / (χs[i₋+1] - χ₊) : (S₊ - S₋) / Δχ
                        S′₊ = i₊ ≥ 2            ? (S₋ - Ss[i₋-2, ik]) / (χ₋ - χs[i₋-2]) : (S₊ - S₋) / Δχ
                        t  = (χ - χ₋) / Δχ
                        t² = t * t
                        t³ = t² * t
                        S  = (2t³-3t²+1)*S₋ + (t³-2t²+t)*Δχ*S′₋ + (-2t³+3t²)*S₊ + (t³-t²)*Δχ*S′₊
                    end
                    I = √(π/(2l+1)) * S / k
                end
            else
                prev = Ss[1, ik] * jl(l, k*χs[1]) # set up first point
                for iτ in 2:length(τs)
                    kχ = k * χs[iτ]
                    halfdτ = halfdτs[iτ-1]
                    _jl = jl(l, kχ)
                    curr = Ss[iτ, ik] * _jl
                    I += halfdτ * (curr + prev)
                    #kχ < l && abs(_jl) < 1e-20 && break # time cut approximation (disabled; unreliable)
                    prev = curr
                end
            end
            Is[ik, il] = I
            #k*τ0 < l && maximum(abs.(I)) < 1e-20 && break # multipole cut approximation (disabled; unreliable)
        end
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

    if normalization == :Cl
        return Cls
    elseif normalization == :Dl
        return @. Cls * ls * (ls+1) / 2π
    else
        error("Normalization $normalization is not :Cl or :Dl")
    end
end

"""
    spectrum_cmb(modes::AbstractVector{<:Symbol}, prob::CosmologyProblem, jl::SphericalBesselCache; normalization = :Cl, unit = nothing, kgrid = nothing, xs = cosgrid(0.0, 1.0; length=300), τcut = 1e-2, l_limber = 10, integrator = TrapezoidalRule(), bgopts = (alg = bgalg(prob), reltol = 1e-7, abstol = 1e-7), ptopts = (alg = ptalg(prob), reltol = 1e-5, abstol = 1e-5), kinterpolate = :cubic, korder = 100, thread = true, verbose = false, kwargs...)

Compute angular CMB power spectra ``Cₗᴬᴮ`` at angular wavenumbers `ls` from the cosmological problem `prob`.
The requested `modes` are specified as a vector of symbols in the form `:AB`, where `A` and `B` are `T` (temperature), `E` (E-mode polarization) or `ψ` (lensing).
If `unit` is `nothing` the spectra are of dimensionless temperature fluctuations relative to the present photon temperature; while if `unit` is a temperature unit the spectra are of dimensionful temperature fluctuations.
Returns a matrix of ``Cₗ`` if `normalization` is `:Cl`, or ``Dₗ = l(l+1)/2π`` if `normalization` is `:Dl`.

# Precision parameters

- `xs`: Grid of ``(τ-τᵢ)/(τ₀-τᵢ)`` specifying the ``τ``-points that will be sampled in line-of-sight integration.
- `τcut`: Remove all earlier times from the line-of-sight integral sampling time points.
- `kgrid`: Grid of ``k``-modes on which the perturbation ODEs will be solved (and then interpolated to a finer grid depending on `Δkτ0`).
- `Δkτ0`: Grid spacing to use when integrating over ``k`` to project to ``ℓ``-space.
- `l_limber`: Use Limber approximation for lensing line-of-sight integrals with equal or greater ``ℓ``.
- `kinterpolate`: Source k-interpolation method: `:cubic` or `:chebyshev` for cubic splines or Chebyshev polynomials.
- `korder`: Number of Chebyshev k-nodes to use when `kinterpolate = :chebyshev`.
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
function spectrum_cmb(modes::AbstractVector{<:Symbol}, prob::CosmologyProblem, jl::SphericalBesselCache; normalization = :Cl, unit = nothing, kgrid = nothing, Δkτ0 = 2π/2, xs = cosgrid(0.0, 1.0; length=300), τcut = 1e-2, l_limber = 10, integrator = TrapezoidalRule(), bgopts = (alg = bgalg(prob), reltol = 1e-7, abstol = 1e-7), ptopts = (alg = ptalg(prob), reltol = 1e-5, abstol = 1e-5), kinterpolate = :cubic, korder = 100, thread = true, verbose = false, kwargs...)
    # Define 1-2-3 indices corresponding for present modes
    iT = 'T' in join(modes) ? 1 : 0
    iE = 'E' in join(modes) ? iT + 1 : 0
    iψ = 'ψ' in join(modes) ? max(iE, iT) + 1 : 0

    # Automatically determine grid if not provided manually
    if isnothing(kgrid)
        kmax = iψ > 0 ? 10000.0 : 1000.0 # higher for lensing
        kgrid = CubicSplineWavenumberGrid(kgrid_default(kmax))
    end

    ls = jl.l
    sol = solve(prob; bgopts, verbose)
    τ0 = getsym(sol, prob.M.τ0)(sol)
    ks_fine = lingrid(minimum(kgrid), maximum(kgrid); step=Δkτ0/τ0) # for k-quadrature after LOS integration

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
    Ss = source_grid(prob, Ss, τs, ks_fine, kgrid, sol.bg; ptopts, verbose, thread)
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
function spectrum_cmb(modes::AbstractVector, prob::CosmologyProblem, jl::SphericalBesselCache, ls::AbstractVector; kwargs...)
    extrema(jl.l) == extrema(ls) || error("jl.l and ls have different extrema $(extrema(jl.l)) != $(ls)")
    spectra_coarse = spectrum_cmb(modes, prob, jl; kwargs...)
    spectra_fine = similar(spectra_coarse, (length(ls), size(spectra_coarse)[2]))
    T = eltype(spectra_coarse)
    ls_coarse = T.(jl.l)
    ls = T.(ls)
    for imode in eachindex(modes)
        spectra_fine[:,imode] = CubicSpline(spectra_coarse[:,imode], ls_coarse)(ls)
    end
    return spectra_fine
end

function spectrum_cmb(mode::Symbol, args...; kwargs...)
    return spectrum_cmb([mode], args...; kwargs...)[:, begin]
end

function cmb_kτ0s(lmin, lmax; Δkτ0 = 2π/2, Δkτ0_S = 8.0, kτ0min = 0.1*lmin, kτ0max = 2*lmax)
    kτ0s_fine = range(kτ0min, kτ0max, step=Δkτ0) # use integer multiple so endpoints are the same
    kτ0s_coarse = range(kτ0s_fine[begin], kτ0s_fine[end], length = Int(floor((kτ0max-kτ0min)/Δkτ0_S+1)))
    kτ0s_coarse[begin] == kτ0s_fine[begin] && kτ0s_coarse[end] == kτ0s_fine[end] || error("different wavenumber endpoints")
    return kτ0s_coarse, kτ0s_fine
end
