using DataInterpolations
using MatterPower
using TwoFAST

"""
    spectrum_primordial(k, h, As, ns=1.0; kp = 0.05/u"Mpc")

Compute the primordial power spectrum
```math
P₀(k) = 2π² Aₛ (k/kₚ)^{nₛ-1} / k³
```
with spectral amplitude `As`, spectral index `ns` and pivot scale wavenumber `kp` at the wavenumber(s) `k`.
"""
function spectrum_primordial(k, h, As, ns=1.0; kp = 0.05/u"Mpc")
    P = 2*π^2 * As ./ k.^3

    # ensure both k and kp to dimensionless wavenumbers k/(H₀/c) before taking ratio
    k = k_dimensionless.(k, h)
    kp = k_dimensionless(kp, h)
    P .*= (k./kp).^(ns-1)

    return P
end
function spectrum_primordial(k, sol::CosmologySolution)
    M = sol.prob.M
    h = have(M, :g) ? sol[M.g.h] : sol[M.h]
    As = have(M, :I) ? sol[M.I.As] : sol[M.As]
    ns = have(M, :I) ? sol[M.I.ns] : sol[M.ns]
    return spectrum_primordial(k, h, As, ns)
end
function spectrum_primordial(k, M::System, pars::Dict)
    return spectrum_primordial(k, pars[M.g.h], pars[M.I.As], pars[M.I.ns])
end
function spectrum_primordial(k, prob::CosmologyProblem)
    M = prob.M
    return spectrum_primordial(k, prob.bg.ps[M.g.h], prob.bg.ps[M.I.As], prob.bg.ps[M.I.ns])
end

function total_symbolic_gauge_invariant_overdensities(M::System, mode::Symbol)
    ρtot = 0
    Δρtot = 0
    mode = String(mode)
    for s in mode
        s = Symbol(s)
        Δ = have(M, s) ? getproperty(getproperty(M, s), :Δ) : getproperty(M, Symbol(:Δ, s))
        length(mode) == 1 && return Δ # short circuit, don't need to do weighting
        ρ = have(M, s) ? getproperty(getproperty(M, s), :ρ) : getproperty(M, Symbol(:ρ, s))
        ρtot += ρ
        Δρtot += ρ*Δ
    end
    return Δρtot / ρtot # e.g. (ρb*Δb+ρc*Δc)/(ρb+ρc)
end

"""
    spectrum_matter([modes,] prob::CosmologyProblem, k[, τ]; kwargs...)

Compute the matter power spectrum
```math
P(k,τ) = P₀(k) |Δ(k,τ)|²
```
of the total gauge-invariant overdensity
```math
Δ = (∑ₛρₛΔₛ) / (∑ₛρₛ)
```
for one or more `modes` at wavenumbers `k` and conformal time(s) `τ` from the problem `prob`.
The problem is solved for the given ``k``, and the matter power spectrum is saved at the given ``τ``.

- `modes` must be `:c` (CDM), `:b` (baryons), `:h` (massive neutrinos), `:m` (matter; equivalent to ``c+b+h``), a vector thereof, or unspecified to use `:m`.
- `k` must be a vector of wavenumbers.
- `τ` must be a single or a vector of conformal times, or unspecified to use ``τ = τ₀`` today.
- `kτini` and `τinimax` specify initial values of ``kτ`` for each perturbation mode, no later than `τinimax` and no earlier than the initial background time.
- `kwargs...` are keyword arguments that are forwarded to `solve(prob, k; kwargs...)`.
"""
function spectrum_matter(modes::AbstractVector, prob::CosmologyProblem, k, τ::AbstractVector; kτini = 1e-2, τinimax = 1e-4, kwargs...)
    ptextraopts = (saveat = τ,)
    sol = solve(prob, k; ptivini = k -> min(kτini / k, τinimax), ptextraopts, kwargs...)
    return spectrum_matter(modes, sol, k, τ)
end
function spectrum_matter(modes::AbstractVector, prob::CosmologyProblem, k; kτini = 1e-2, τinimax = 1e-4, kwargs...)
    ptextraopts = (save_everystep = false, save_start = false, save_end = true)
    sol = solve(prob, k; ptivini = k -> min(kτini / k, τinimax), ptextraopts, kwargs...)
    return spectrum_matter(modes, sol, k)
end

"""
    spectrum_matter([modes,] sol::CosmologySolution, k[, τ]; kwargs...)

Compute the matter power spectrum in the same way, but interpolate between wavenumbers and times already stored in the solution `sol`.
"""
function spectrum_matter(modes::AbstractVector, sol::CosmologySolution, k::AbstractVector, τ::AbstractVector)
    M = sol.prob.M
    S = map(mode -> total_symbolic_gauge_invariant_overdensities(M, mode), modes)
    P0 = spectrum_primordial(k, sol)
    P0 = reshape(P0, 1, 1, :)
    P = P0 .* sol(S, τ, k) .^ 2
    return P
end
spectrum_matter(modes::AbstractVector, sol::CosmologySolution, k; kwargs...) = spectrum_matter(modes, sol, k, sol.bg.t[end]; kwargs...) # fallback without time
spectrum_matter(modes::AbstractVector, probsol, k, τ::Number; kwargs...) = spectrum_matter(modes, probsol, k, [τ])[:, 1, :] # fallback with single time
spectrum_matter(mode::Symbol, probsol, args...; kwargs...) = selectdim(spectrum_matter([mode], probsol, args...; kwargs...), 1, 1) # fallback with single mode specified
spectrum_matter(probsol::Union{CosmologyProblem, CosmologySolution}, args...; kwargs...) = spectrum_matter(:m, probsol, args...; kwargs...) # fallback with modes unspecified

"""
    spectrum_matter([modes,] prob::CosmologyProblem, k::NTuple{2, Number}, τs = nothing; sourceopts = (atol = 4.0, rtol = 4e-3), coarse_length = 9, kwargs...)

Compute the matter power spectrum on the interval ``k`` with adaptively chosen wavenumbers.
Returns wavenumbers and power spectrum values.

The interval is first divided into a grid with `coarse_length` logarithmically spaced wavenumbers.
It is then adaptively refined with the absolute and relative tolerances in `sourceopts`.
"""
function spectrum_matter(modes::AbstractVector, prob::CosmologyProblem, k::NTuple{2, Number}, τs = nothing; sourceopts = (atol = 4.0, rtol = 4e-3), coarse_length = 9, kwargs...)
    # Initial coarse k-grid
    kmin, kmax = k
    ks = exp.(range(log(kmin), log(kmax), length = coarse_length))
    ks[begin] = kmin # exp(log(k)) ≠ k with floats; ensure ends are exactly what the user passed
    ks[end] = kmax # exp(log(k)) ≠ k with floats; ensure ends are exactly what the user passed

    ktransform = (log, exp)
    Ss = map(mode -> total_symbolic_gauge_invariant_overdensities(prob.M, mode), modes)
    ks, Δs = source_grid_adaptive(prob, Ss, τs, ks; ktransform, sourceopts..., kwargs...)
    Ps = stack(Δs[1, :]) .^ 2
    P0s = spectrum_primordial(ks, prob)
    for iS in eachindex(Ss)
        Ps[iS, :] .*= P0s
    end
    return ks, Ps
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
    Ωm0 = sol[M.m.Ω₀]
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
    R = 1 / k_dimensionless(1 / R, sol.bg) # make dimensionless
    return MatterPower.sigma2(Pf, R)
end
"""
    stddev_matter(sol::CosmologySolution, R)

Compute the standard deviation ``√(⟨δ²⟩)`` of the *linear* matter density field with a top-hat filter with radius `R`.
"""
stddev_matter(sol::CosmologySolution, R) = √(variance_matter(sol, R))

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

abstract type AbstractInterpolator{T} end

Base.extrema(interp::AbstractInterpolator) = (minimum(interp), maximum(interp))

struct CubicSplineInterpolator{T, F <: Function} <: AbstractInterpolator{T}
    xs::Vector{T} # points in input domain: x = f⁻¹(y) (e.g. wavenumbers k)
    ys::Vector{T} # points in interpolation domain: y = f(x)
    f::F
end

Base.minimum(interp::CubicSplineInterpolator) = interp.xs[begin]
Base.maximum(interp::CubicSplineInterpolator) = interp.xs[end]

function CubicSplineInterpolator(xs; f = identity)
    issorted(xs) || throw(ArgumentError("Input points must be sorted in ascending order"))
    xs = collect(xs) # to array
    ys = f.(xs)
    return CubicSplineInterpolator(xs, ys, f)
end
function CubicSplineInterpolator(xmin, xmax, n)
    xs = range(xmin, xmax, length = n+1)
    return CubicSplineInterpolator(xs)
end


"""
    source_grid(Ss_coarse::AbstractMatrix, ks_coarse, ks_fine; ktransform = identity, thread = true)

Interpolate values `Ss_coarse` of source functions ``S(τ,k)`` from a coarse wavenumber grid `ks_coarse` to a fine grid `ks_fine`.
The interpolation is cubic spline in `ktransform(k)` (e.g. `identity` for interpolation in ``k`` or `log` for interpolation in ``\\ln k``).
Conformal times are unchanged.
"""
function source_kinterp!(out::AbstractVector, Ss_coarse::AbstractVector, kinterp::CubicSplineInterpolator, ys_fine)
    interp = CubicSpline(Ss_coarse, kinterp.ys)
    out .= interp.(ys_fine)
    return out
end
function source_kinterp!(out::AbstractMatrix, Ss_coarse::AbstractMatrix, kinterp::AbstractInterpolator, ks_fine; thread = true)
    ks_coarse = kinterp.xs
    size(Ss_coarse, 1) == size(out, 1) || error("out has first dimension with length $(size(out, 1)), but Ss_coarse has $(size(Ss_coarse, 1))")
    size(Ss_coarse, 2) == length(ks_coarse) || error("Length of coarse k-grid does not match source array")
    ys_fine = kinterp.f.(ks_fine)
    @inbounds @tasks for i in 1:size(Ss_coarse, 1)
        @set scheduler = thread ? :dynamic : :static
        source_kinterp!(@view(out[i, :]), @view(Ss_coarse[i, :]), kinterp, ys_fine)
    end
    return out
end
function source_kinterp(Ss_coarse::AbstractMatrix, kinterp::AbstractInterpolator, ks_fine; kwargs...)
    Ss_fine = similar(Ss_coarse, size(Ss_coarse, 1), length(ks_fine))
    source_kinterp!(Ss_fine, Ss_coarse, kinterp, ks_fine; kwargs...)
    return Ss_fine
end

function source_eltype(Ss, T)
    if Ss isa StaticVector
        SVector{length(Ss), T}
    elseif Ss isa AbstractVector
        Vector{T}
    else
        T
    end
end

"""
    source_grid(prob::CosmologyProblem, Ss, τs, ks[, bgsol]; bgopts = (), ptopts = (), thread = true, verbose = false)

Compute and evaluate source functions ``S(τ,k)`` with symbolic expressions `Ss` on a grid with conformal times `τs` and wavenumbers `ks` from the problem `prob`.
Returns a matrix of size `(Nτ, Nk)`, where each element is a vector of length `NS = length(Ss)` holding all source values at that `(τ, k)` point.

The options `bgopts` and `ptopts` are passed to the background and perturbation solves.
"""
function source_grid(prob::CosmologyProblem, Ss, τs, ks, bgsol::ODESolution; ptopts = (), thread = true, verbose = false)
    getSs = getsym(prob.pt, Ss)
    T = source_eltype(Ss, eltype(bgsol))
    Ss = Matrix{T}(undef, length(τs), length(ks))
    minimum(τs) ≥ bgsol.t[begin] && maximum(τs) ≤ bgsol.t[end] || error("input τs and computed background solution have different timespans")
    function output_func(sol, ik)
        vals = getSs(sol)
        @inbounds for iτ in eachindex(vals)
            Ss[iτ, ik] = T(vals[iτ])
        end
        return nothing
    end
    solvept(prob.pt, bgsol, ks; output_func, saveat = τs, ptopts..., thread, verbose)
    return Ss
end
function source_grid(prob::CosmologyProblem, Ss, τs, ks; bgopts = (), verbose = false, kwargs...)
    bgsol = solvebg(prob.bg; bgopts..., verbose)
    return source_grid(prob, Ss, τs, ks, bgsol)
end

function source_grid(prob::CosmologyProblem, Ss, τs, ks, kinterp::AbstractInterpolator, args...; thread = true, kwargs...)
    Ss = source_grid(prob, Ss, τs, kinterp.xs, args...; thread, kwargs...) # solve perturbation on coarse k-grid
    Ss = source_kinterp(Ss, kinterp, ks; thread) # interpolate to requested fine k-grid
    return Ss
end


# TODO: Hermite interpolation
# TODO: create SourceFunction type that does k and τ interpolation?
"""
    source_grid_adaptive(prob::CosmologyProblem, Ss, τs, ks[, bgsol]; bgopts = (), ptopts = (), refine = true, ktransform = (identity, identity), sort = true, thread = true, verbose = false, kwargs...)

Adaptively compute and evaluate source functions ``S(τ,k)`` with symbolic expressions `Ss` on a grid with fixed conformal times `τs`, but adaptively refined grid of wavenumbers from the problem `prob`.
The source functions are first evaluated on the (coarse) initial grid `ks`.
Each subinterval ``(k₁, k₂)`` of `ks` is then adaptively refined until the linear interpolation ``Sᵢ = (S(k₁)+S(k₂))/2`` to the midpoint ``k=(k₁+k₂)/2`` approximates the actual value ``S(k)`` there within some tolerance.
The comparison ``Sᵢ ≈ S`` is done with `isapprox(Sᵢ, S; kwargs...)`, where `S` and `Sᵢ` are vectors of all source values at a given ``(τ, k)``.
It receives the keyword arguments `kwargs` passed to this function, so `atol`, `rtol` and/or `norm` can be specified to tune the tolerance.

If `τs` is nothing, the source function is evaluated at the final time only (today).

Returns the refined wavenumbers `ks` sorted in ascending order and a matrix of size `(length(τs), length(ks))` with the corresponding source function values, where each element is a vector of length `length(Ss)`.
If not `sort`, the wavenumbers and source function values are instead left in the order in which they were inserted in the refinement process.

The options `bgopts` and `ptopts` are passed to the background and perturbation solves.
"""
function source_grid_adaptive(prob::CosmologyProblem, Ss, τs, ks, bgsol::ODESolution; ptopts = (), refine = true, ktransform = (identity, identity), sort = true, thread = true, verbose = false, kwargs...)
    length(ks) ≥ 2 || error("Initial k-grid must have at least 2 values")

    if isnothing(τs)
        Nτs = 1 # final time only
        ptsaveopts = (save_everystep = false, save_start = false, save_end = true)
    else
        Nτs = length(τs)
        ptsaveopts = (saveat = τs,)
    end

    ptprobgen = setuppt(prob.pt, bgsol)

    getSs = getsym(prob.pt, Ss)
    T = source_eltype(Ss, eltype(bgsol))
    function sourcek!(k, ik, Ss)
        ptprob = ptprobgen(k)
        ptsol = solvept(ptprob; ptsaveopts..., ptopts...)
        vals = getSs(ptsol)
        @inbounds for iτ in eachindex(vals)
            Ss[iτ, ik] = T(vals[iτ])
        end
        return nothing
    end

    kmin, kmax = extrema(ks)
    f, f⁻¹ = ktransform
    kmin ≈ f⁻¹(f(kmin)) || error("ktransform is not a tuple of inverse functions")
    kmax ≈ f⁻¹(f(kmax)) || error("ktransform is not a tuple of inverse functions")

    counter = Atomic{Int}(length(ks)-1)
    idx = Atomic{Int}(length(ks))

    ninitks = length(ks)
    nmaxks = 1024
    Ss = Matrix{T}(undef, Nτs, nmaxks)
    ks = resize!(collect(ks), nmaxks)
    iτs = 1 : max(Nτs - 1, 1) # exclude χ=0 from refinement comparison (where some CMB sources with 1/χ diverge); except if it is the only point (e.g. matter power spectrum is well-defined)

    queue = Channel{Tuple{Int, Int}}(nmaxks)
    tasks = Task[]
    thread && sizehint!(tasks, nmaxks)
    for i in 1:ninitks
        task = @spawnif begin
        sourcek!(ks[i], i, Ss)
        i ≥ 2 && put!(queue, (i-1, i))
        verbose && println("Solved k = $(ks[i]) on thread $(threadid()) to $i total points")
        end thread
        thread && push!(tasks, task)
    end

    while refine # equivalent to for (i1, i3) in queue, but rely on sentinel value instead of closing queue
        i1, i3 = take!(queue) # indices corresponding to left and right of interval
        i1 == -1 && break # kill on sentinel value
        task = @spawnif begin
        i2 = atomic_add!(idx, +1) + 1 # index corresponding to midpoint of interval; atomic_add! returns old value of idx
        if i2 > nmaxks
            put!(queue, (-1, -1)) # kill with sentinel value if no more space for ks
            return
        end
        f1 = f(ks[i1]) # left point
        f3 = f(ks[i3]) # right points
        f2 = (f1 + f3) / 2 # middle point
        k2 = f⁻¹(f2) # f transforms e.g. k → log(k)
        ks[i2] = k2
        sourcek!(ks[i2], i2, Ss)

        verbose && println("Refined k-grid between [$(ks[i1]), $(ks[i3])] on thread $(threadid()) to $i2 total points")

        # compare exactly solved S2 with linear interpolation from S1 and S3
        # refine if interpolation is not good enough
        if refine && any(iτ -> !isapprox(Ss[iτ,i2], (Ss[iτ,i1]+Ss[iτ,i3])/2; kwargs...), iτs)
            atomic_add!(counter, +2)
            put!(queue, (i2, i3)) # refine right subinterval
            put!(queue, (i1, i2)) # refine left subinterval
        end

        atomic_add!(counter, -1) == 1 && put!(queue, (-1, -1)) # kill with sentinel value if counter reaches 0 (atomic_add! returns old value)
        end thread
        thread && push!(tasks, task)
    end

    thread && foreach(wait, tasks) # all tasks must finish
    close(queue) # clean up

    idx[] > nmaxks && error("Source function refinement needs more than $nmaxks k-refinements. Reduce refinement criteria.")

    # sort according to k
    ks = ks[1:idx[]]
    if sort
        is = sortperm(ks)
        ks = ks[is]
        Ss = Ss[:, is]
    else
        Ss = Ss[:, 1:length(ks)]
    end

    return ks, Ss
end

# Dispatch without background solution
function source_grid_adaptive(prob::CosmologyProblem, Ss, τs, ks; bgopts = (), kwargs...)
    bgsol = solvebg(prob.bg; bgopts...)
    return source_grid_adaptive(prob, Ss, τs, ks, bgsol; kwargs...)
end


struct ChebyshevInterpolator{T <: Real, F} <: AbstractInterpolator{T}
    xs::Vector{T} # points in input domain: x = f⁻¹(y) (e.g. wavenumbers k)
    ys::Vector{T} # points in interpolation domain: y = f(x)
    f::F
end

Base.minimum(interp::ChebyshevInterpolator) = interp.xs[end] # stored grid is from high-to-low x
Base.maximum(interp::ChebyshevInterpolator) = interp.xs[begin]

function ChebyshevInterpolator(xmin, xmax, order; f = identity, f⁻¹ = nothing)
    xmax > xmin || throw(ArgumentError("Interval $((xmin, xmax)) is not sorted"))
    ymin, ymax = f(xmin), f(xmax)
    ys = chebpoints(order, ymin, ymax)
    issorted(ys; rev = true) || throw(ArgumentError("Domain transformation f(x) is not monotonically increasing"))
    if f == identity && isnothing(f⁻¹)
        f⁻¹ = identity
    end
    if isnothing(f⁻¹)
        # invert numerically
        xs = map(ys) do y
            prob = IntervalNonlinearProblem((x, _) -> f(x) - y, (xmin, xmax))
            sol = solve(prob)
            return sol.u
        end
    else
        # invert analytically
        xs = f⁻¹.(ys)
    end
    xs[end] ≈ xmin && xs[begin] ≈ xmax || throw(ArgumentError("f(x) and f⁻¹(x) are not inverses"))
    xs[end], xs[begin] = xmin, xmax  # prevent floating point bounds errors from f⁻¹(f(k))
    return ChebyshevInterpolator(xs, ys, f)
end

Base.length(interp::AbstractInterpolator) = length(interp.xs)
order(interp::AbstractInterpolator) = length(interp) - 1

function source_kinterp!(out::AbstractVector, Ss_coarse::AbstractVector, kinterp::ChebyshevInterpolator, ys_fine)
    if any(any(!isfinite(Sᵢ) for Sᵢ in S) for S in Ss_coarse)
        out .*= NaN # set to NaN instead of crashing in chebinterp
        return
    end
    ymin, ymax = kinterp.ys[end], kinterp.ys[begin] # stored grid is from high-to-low x
    interp = chebinterp(Ss_coarse, ymin, ymax)
    out .= interp.(ys_fine)
    return out
end

# Special dispatch for returning a vector of interpolation objects (for testing)
function source_grid_interp(prob::CosmologyProblem, S, τs, kinterp::ChebyshevInterpolator, args...; kwargs...)
    Ss = source_grid(prob, S, τs, kinterp.xs, args...; kwargs...)
    ymin, ymax = kinterp.ys[end], kinterp.ys[begin]
    return [chebinterp(Ss[i, :], ymin, ymax) for i in eachindex(τs)]
end


struct PiecewiseChebyshevInterpolator{T <: Real, G <: Tuple} <: AbstractInterpolator{T}
    subgrids::G # NTuple of ChebyshevInterpolator in ascending x-order
    xs::Vector{T} # all unique coarse x-values, in descending x-order
    iranges::Vector{UnitRange{Int}} # index range into xs for each subgrid
end

Base.minimum(interp::PiecewiseChebyshevInterpolator) = interp.xs[end] # stored from high-to-low x
Base.maximum(interp::PiecewiseChebyshevInterpolator) = interp.xs[begin]

function PiecewiseChebyshevInterpolator(xbreaks, orders; f = identity, f⁻¹ = identity)
    N = length(orders) # number of piecewise subgrids
    length(xbreaks) == N + 1 || throw(ArgumentError("Need $(N+1) x-breaks for $N intervals, got $(length(xbreaks))"))
    if !(f isa Tuple)
        f = ntuple(_ -> f, N)
    end
    if !(f⁻¹ isa Tuple)
        f⁻¹ = ntuple(_ -> f⁻¹, N)
    end
    length(f)  == N || throw(ArgumentError("Need $N f, got $(length(f))"))
    length(f⁻¹) == N || throw(ArgumentError("Need $N f⁻¹, got $(length(f⁻¹))"))

    subgrids = ntuple(j -> ChebyshevInterpolator(xbreaks[j], xbreaks[j+1], orders[j]; f = f[j], f⁻¹ = f⁻¹[j]), N)
    xs = reduce(vcat, (subgrids[j].xs[2:end] for j in N-1:-1:1); init = subgrids[end].xs) # combine unique x-points in descending order (boundaries share x points)
    iranges = Vector{UnitRange{Int}}(undef, N)
    i = 1
    for j in N:-1:1
        n = length(subgrids[j].xs)
        iranges[j] = i : i + n - 1 # index range into xs corresponding to subgrid j
        i += n - 1
    end
    return PiecewiseChebyshevInterpolator{eltype(xs), typeof(subgrids)}(subgrids, xs, iranges)
end

function source_kinterp(Ss_coarse::AbstractMatrix, kinterp::PiecewiseChebyshevInterpolator, ks_fine; thread = true)
    Ss_fine = similar(Ss_coarse, size(Ss_coarse, 1), length(ks_fine))
    @inbounds @tasks for j in eachindex(kinterp.subgrids)
        @set scheduler = thread ? :dynamic : :static
        subgrid = kinterp.subgrids[j]
        kmin, kmax = extrema(subgrid)
        in_range = findall(k -> kmin ≤ k ≤ kmax, ks_fine)
        irange = kinterp.iranges[j]
        source_kinterp!(@view(Ss_fine[:, in_range]), @view(Ss_coarse[:, irange]), subgrid, ks_fine[in_range]; thread)
    end
    return Ss_fine
end
