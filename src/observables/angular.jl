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

function SphericalBesselCache(ls::AbstractVector; xmax = 10*ls[end], dx = 2π/15, hermite = true)
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
    los_integrate(Ss::AbstractMatrix{T}, ls::AbstractVector, τs::AbstractVector, ks::AbstractVector, jl::SphericalBesselCache; l_limber = typemax(Int), integrator = TrapezoidalRule(), thread = true, verbose = false) where {T <: Real}

For the given `ls` and `ks`, compute the line-of-sight-integrals
```math
Iₗ(k) = ∫dτ S(k,τ) jₗ(k(τ₀-τ))
```
over the source function values `Ss` against the spherical Bessel functions ``jₗ(x)`` cached in `jl`.
The element `Ss[i,j]` holds the source function value ``S(τᵢ, kⱼ)``.
The limber approximation
```math
Iₗ ≈ √(π/(2l+1)) S(τ₀-(l+1/2)/k, k)
```
is used for `l ≥ l_limber`.
"""
function los_integrate(Ss::AbstractMatrix{T}, ls::AbstractVector, τs::AbstractVector, ks::AbstractVector, jl::SphericalBesselCache; l_limber = typemax(Int), integrator = TrapezoidalRule(), thread = true, verbose = false) where {T <: Real}
    # Julia is column-major; make sure innermost loop indices appear first in slice expressions (https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-column-major)
    @assert size(Ss, 1) == length(τs) "size(Ss, 1) = $(size(Ss, 1)) and length(τs) = $(length(τs)) differ"
    @assert size(Ss, 2) == length(ks) "size(Ss, 2) = $(size(Ss, 2)) and length(ks) = $(length(ks)) differ"
    @assert jl.x[begin] ≤ 0 "jl.x[begin] < 0"
    @assert jl.x[end] ≥ ks[end]*τs[end] "jl.x[end] < kmax*τmax"
    @assert all(isfinite, Ss) "Ss contain NaN or Inf"
    τs = collect(τs) # force array to avoid floating point errors with ranges in following χs due to (e.g. tiny negative χ)
    τ0 = τs[end]
    χs = τ0 .- τs
    halfdτs = 0.5 .* (τs[begin+1:end] .- τs[begin:end-1]) # precompute before loops
    Is = zeros(T, length(ks), length(ls))

    verbose && l_limber < typemax(Int) && println("Using Limber approximation for l ≥ $l_limber")

    # TODO: skip and set jl to zero if l ≳ kτ0 or another cutoff?
    @fastmath @inbounds @tasks for il in eachindex(ls) # parallellize independent loop iterations
        @set scheduler = thread ? :dynamic : :serial
        l = ls[il]
        verbose && print("\rLOS integrating with l = $l")
        for ik in reverse(eachindex(ks))
            k = ks[ik]
            I = 0.0
            if l ≥ l_limber
                χ = (l+1/2) / k
                if χ ≤ χs[1] # otherwise χ > χini > χrec and source function is definitely zero
                    # interpolate between two closest points in saved array
                    i₋ = searchsortedfirst(τs, τ0 - χ)
                    i₊ = i₋ - 1 # χ is sorted in descending order
                    χ₋, χ₊ = χs[i₋], χs[i₊] # now χ₋ < χ < χ₊
                    S₋, S₊ = Ss[i₋, ik], Ss[i₊, ik]
                    S = S₋ + (S₊-S₋) * (χ-χ₋) / (χ₊-χ₋)
                    I = √(π/(2l+1)) * S / k
                end
            else
                prev = Ss[1, ik] * jl(l, k*χs[1]) # set up first point
                for iτ in 2:length(τs)
                    kχ = k * χs[iτ]
                    halfdτ = halfdτs[iτ-1]
                    _jl = jl(l, kχ)
                    curr = Ss[iτ, ik] * _jl
                    dI = halfdτ * (curr + prev)
                    I += dI
                    kχ < l && abs(_jl) < 1e-20 && break # time cut approximation
                    prev = curr
                end
            end
            Is[ik, il] = I
            k*τ0 < l && abs(I) < 1e-20 && break # multipole cut approximation
        end
    end
    verbose && println()

    return Is
end
function los_integrate(sol::CosmologySolution, ls::AbstractVector, τs::AbstractVector, ks::AbstractVector, S, jl::SphericalBesselCache; ktransform = identity, kwargs...) # TODO: Ss
    Ss = [S]
    Ss = source_grid(sol, Ss, τs)
    Ss = source_grid(Ss, sol.ks, ks; ktransform)
    Ss[:, end, :] .= 0.0 # may be NaNs today, but jl(0) = 0, so today is always 0 in the line-of-sight integral
    Ss = @view Ss[1, :, :]
    return los_integrate(Ss, ls, τs, ks, jl; kwargs...)
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
    spectrum_cmb(modes::AbstractVector{<:Symbol}, prob::CosmologyProblem, jl::SphericalBesselCache; normalization = :Cl, unit = nothing, kτ0s = 0.1*jl.l[begin]:2π/2:10*jl.l[end], xs = 0.0:0.0008:1.0, l_limber = 10, integrator = TrapezoidalRule(), bgopts = (alg = bgalg(prob), reltol = 1e-7, abstol = 1e-7), ptopts = (alg = ptalg(prob),, reltol = 1e-5, abstol = 1e-5), sourceopts = (rtol = 1e-3, atol = 0.9), downsampleopts = (Ttol = 4e-3, Etol = 2e-4, ψtol = 1e-3), coarse_length = 9, thread = true, verbose = false, kwargs...)

Compute angular CMB power spectra ``Cₗᴬᴮ`` at angular wavenumbers `ls` from the cosmological problem `prob`.
The requested `modes` are specified as a vector of symbols in the form `:AB`, where `A` and `B` are `T` (temperature), `E` (E-mode polarization) or `ψ` (lensing).
If `unit` is `nothing` the spectra are of dimensionless temperature fluctuations relative to the present photon temperature; while if `unit` is a temperature unit the spectra are of dimensionful temperature fluctuations.
Returns a matrix of ``Cₗ`` if `normalization` is `:Cl`, or ``Dₗ = l(l+1)/2π`` if `normalization` is `:Dl`.

The lensing line-of-sight integral uses the Limber approximation for `l ≥ l_limber`.

Source functions are computed on a ``k``-grid that is adaptively refined from an initial grid with size `coarse_length`.
The refinement criterion is controlled with `sourceopts`.

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
function spectrum_cmb(modes::AbstractVector{<:Symbol}, prob::CosmologyProblem, jl::SphericalBesselCache; normalization = :Cl, unit = nothing, kτ0s = 0.1*jl.l[begin]:2π/2:10*jl.l[end], xs = 0.0:0.0008:1.0, l_limber = 10, integrator = TrapezoidalRule(), bgopts = (alg = bgalg(prob), reltol = 1e-7, abstol = 1e-7), ptopts = (alg = ptalg(prob), reltol = 1e-5, abstol = 1e-5), sourceopts = (rtol = 1e-3, atol = 0.9), downsampleopts = (Ttol = 4e-3, Etol = 2e-4, ψtol = 1e-3), coarse_length = 9, thread = true, verbose = false, kwargs...)
    ls = jl.l
    sol = solve(prob; bgopts, verbose)
    τ0 = getsym(sol, prob.M.τ0)(sol)
    ks_fine = collect(kτ0s ./ τ0)

    τs = sol.bg.t # by default, use background (thermodynamics) time points for line of sight integration
    if !isnothing(xs)
        # use user's array of x = (τ-τi)/(τ0-τi)
        xs[begin] == 0 || error("xs begins with $(xs[begin]), but should begin with 0")
        xs[end] == 1 || error("xs ends with $(xs[end]), but should end with 1")
        τs = τs[begin] .+ (τs[end] .- τs[begin]) .* xs
    end

    # Integrate perturbations to calculate source function on coarse k-grid
    iT = 'T' in join(modes) ? 1 : 0
    iE = 'E' in join(modes) ? iT + 1 : 0
    iψ = 'ψ' in join(modes) ? max(iE, iT) + 1 : 0
    Ss = [prob.M.ST, prob.M.SE_kχ², prob.M.Sψ]
    ks_coarse = range(ks_fine[begin], ks_fine[end]; length = coarse_length)
    ks_coarse, Ss = source_grid_adaptive(prob, Ss, τs, ks_coarse, sol.bg; ptopts, verbose, thread, sourceopts...) # TODO: pass kτ0 and x

    Θls = zeros(eltype(Ss), max(iT, iE, iψ), length(ks_fine), length(ls))
    if iT > 0
        STs, τTs = source_grid_downsample(Ss[1, :, :], τs; tol = downsampleopts.Ttol) # downsample in τ
        verbose && println("Downsampled T source function from ", length(τs), " to ", length(τTs), " time points")
        STs = source_grid(STs, ks_coarse, ks_fine; thread) # upsample in k
        Θls[iT, :, :] .= los_integrate(STs, ls, τTs, ks_fine, jl; integrator, verbose, thread, kwargs...)
    end
    if iE > 0
        SEs, τEs = source_grid_downsample(Ss[2, :, :], τs; tol = downsampleopts.Etol) # downsample in τ
        verbose && println("Downsampled E source function from ", length(τs), " to ", length(τEs), " time points")
        @. SEs ./= (ks_coarse' * (τ0-τEs))^2
        SEs[end, :] .= 0.0 # can be Inf, but is always weighted by zero-valued spherical Bessel function in LOS integration
        SEs = source_grid(SEs, ks_coarse, ks_fine; thread) # upsample in k
        Θls[iE, :, :] .= transpose(@. √((ls+2)*(ls+1)*(ls+0)*(ls-1))) .* los_integrate(SEs, ls, τEs, ks_fine, jl; integrator, verbose, thread, kwargs...)
    end
    if iψ > 0
        Sψs, τψs = source_grid_downsample(Ss[3, :, :], τs; tol = downsampleopts.ψtol) # downsample in τ
        verbose && println("Downsampled ψ source function from ", length(τs), " to ", length(τψs), " time points")
        Sψs = source_grid(Sψs, ks_coarse, ks_fine; thread) # upsample in k
        Sψs[end, :] .= 0.0 # contains Inf/NaN, but will be weighted by 0 from jl in LOS integral
        Θls[iψ, :, :] .= los_integrate(Sψs, ls, τψs, ks_fine, jl; l_limber, integrator, verbose, thread, kwargs...)
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

    spectra = zeros(eltype(Ss[1,1,1] * P0s[1] * factor^2), length(ls), length(modes)) # Cls or Dls
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
    for imode in eachindex(modes)
        spectra_fine[:,imode] = QuadraticSpline(spectra_coarse[:,imode], jl.l)(ls)
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
