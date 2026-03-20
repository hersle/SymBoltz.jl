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
    Ps = Δs[:, 1, :] .^ 2
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

# TODO: take in ks and N_skip_interp = 1, 2, 3, ... for more efficient? same for τ?
"""
    source_grid(sol::CosmologySolution, Ss::AbstractVector, τs; thread = true)

Evaluate and return source functions ``S(τ,k)`` from the solution `sol`.
The source functions are given by symbolic expressions `Ss`, and evaluated on a grid with conformal times `τs` and the same wavenumbers as `sol` is solved for.
"""
function source_grid(sol::CosmologySolution, Ss::AbstractVector, τs; thread = true)
    # Evaluate integrated perturbations on coarse grid
    ks = sol.ks
    getSs = map(S -> getsym(sol.prob.pt, S), Ss)
    Ss = similar(sol.bg, length(Ss), length(τs), length(ks))
    @tasks for ik in eachindex(ks)
        @set scheduler = thread ? :dynamic : :static
        for iS in eachindex(getSs)
            Ss[:, :, ik] .= permutedims(getSs[iS](sol.pts[ik]))
        end
    end
    return Ss
end

"""
    source_grid(Ss_coarse::AbstractArray, ks_coarse, ks_fine; ktransform = identity, thread = true)

Interpolate values `Ss_coarse` of source functions ``S(τ,k)`` from a coarse wavenumber grid `ks_coarse` to a fine grid `ks_fine`.
The interpolation is linear in `ktransform(k)` (e.g. `identity` for interpolation in ``k`` or `log` for interpolation in ``\\ln k``.
Conformal times are unchanged.
"""
function source_grid(Ss_coarse::AbstractArray, ks_coarse, ks_fine; ktransform = identity, thread = true)
    size_coarse = size(Ss_coarse)
    size_fine = (size_coarse[1], size_coarse[2], length(ks_fine))
    Ns, Nτ, _ = size(Ss_coarse)

    Ss_fine = similar(Ss_coarse, size_fine)
    xs_coarse = ktransform.(ks_coarse) # TODO: user should just pass different ks as input instead
    xs_fine = ktransform.(ks_fine)
    @tasks for iτ in 1:Nτ
        @set scheduler = thread ? :dynamic : :static
        for iS in 1:Ns
            interp = LinearInterpolation(@view(Ss_coarse[iS, iτ, :]), xs_coarse)
            Ss_fine[iS, iτ, :] .= interp.(xs_fine)
        end
    end
    return Ss_fine
end

"""
    source_grid(prob::CosmologyProblem, Ss::AbstractArray, τs, ks; bgopts = (), ptopts = (), thread = true, verbose = false)

Compute and evaluate source functions ``S(τ,k)`` with symbolic expressions `Ss` on a grid with conformal times `τs` and wavenumbers `ks` from the problem `prob`.

The options `bgopts` and `ptopts` are passed to the background and perturbation solves.
"""
function source_grid(prob::CosmologyProblem, Ss::AbstractArray, τs, ks; bgopts = (), ptopts = (), thread = true, verbose = false)
    bgsol = solvebg(prob.bg; bgopts..., verbose)
    getSs = map(S -> getsym(prob.pt, S), Ss)
    Ss = similar(bgsol, length(Ss), length(τs), length(ks))
    minimum(τs) ≥ bgsol.t[begin] && maximum(τs) ≤ bgsol.t[end] || error("input τs and computed background solution have different timespans")
    function output_func(sol, ik)
        for iS in eachindex(getSs)
            Ss[iS, :, ik] .= getSs[iS](sol)
        end
        return nothing
    end
    solvept(prob.pt, bgsol, ks; output_func, saveat = τs, ptopts..., thread, verbose)
    return Ss
end

# TODO: Hermite interpolation
# TODO: create SourceFunction type that does k and τ interpolation?
"""
    source_grid_adaptive(prob::CosmologyProblem, Ss::AbstractVector, τs, ks[, bgsol]; bgopts = (), ptopts = (), refine = true, ktransform = (identity, identity), sort = true, thread = true, verbose = false, kwargs...)

Adaptively compute and evaluate source functions ``S(τ,k)`` with symbolic expressions `Ss` on a grid with fixed conformal times `τs`, but adaptively refined grid of wavenumbers from the problem `prob`.
The source functions are first evaluated on the (coarse) initial grid `ks`.
Each subinterval ``(k₁, k₂)`` of `ks` is then adaptively refined until the linear interpolation ``Sᵢ = (S(k₁)+S(k₂))/2`` to the midpoint ``k=(k₁+k₂)/2`` approximates the actual value ``S(k)`` there within some tolerance.
The comparison ``Sᵢ ≈ S`` is done with `isapprox(Sᵢ, S; kwargs...)`, where `S` and `Sᵢ` are vectors with the (conformal) timeseries of the source function for that wavenumber.
It receives the keyword arguments `kwargs` passed to this function, so `atol`, `rtol` and/or `norm` can be specified to tune the tolerance.

If `τs` is nothing, the source function is evaluated at the final time only (today).

Returns the refined wavenumbers sorted in ascending order and a grid with the corresponding source function values.
If not `sort`, the wavenumbers and source function values are instead left in the order in which they were inserted in the refinement process.

The options `bgopts` and `ptopts` are passed to the background and perturbation solves.
"""
function source_grid_adaptive(prob::CosmologyProblem, Ss::AbstractVector, τs, ks, bgsol::ODESolution; ptopts = (), refine = true, ktransform = (identity, identity), sort = true, thread = true, verbose = false, kwargs...)
    length(ks) ≥ 2 || error("Initial k-grid must have at least 2 values")

    if isnothing(τs)
        Nτs = 1
        ptsaveopts = (save_everystep = false, save_start = false, save_end = true)
    else
        Nτs = length(τs)
        ptsaveopts = (saveat = τs,)
    end

    ptprobgen = setuppt(prob.pt, bgsol)

    getSs = map(S -> getsym(prob.pt, S), Ss)
    function sourcek!(k, ik, Ss)
        ptprob = ptprobgen(k)
        ptsol = solvept(ptprob; ptsaveopts..., ptopts...)
        for iS in eachindex(getSs)
            Ss[iS, :, ik] .= getSs[iS](ptsol)
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
    Ss = similar(bgsol, length(Ss), Nτs, nmaxks)
    ks = resize!(collect(ks), nmaxks)
    iτs = 1 : max(size(Ss, 2) - 1, 1) # exclude χ=0 from refinement comparison (where some CMB sources with 1/χ diverge); except if it is the only point (e.g. matter power spectrum is well-defined)

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

    while true # equivalent to for (i1, i2) in queue, but rely on sentinel value instead of closing queue
        i1, i2 = take!(queue)
        i1 == -1 && break # kill on sentinel value
        task = @spawnif begin
        i = atomic_add!(idx, +1) + 1 # atomic_add! returns old value of idx
        if i > nmaxks
            put!(queue, (-1, -1)) # kill with sentinel value if no more space for ks
            return
        end
        k1 = f(ks[i1]) # e.g. k1 → log(k1)
        k2 = f(ks[i2])
        k = (k1 + k2) / 2
        k = f⁻¹(k) # e.g. log(k1) → exp(log(k1)) = k1
        ks[i] = k
        sourcek!(k, i, Ss)

        verbose && println("Refined k-grid between [$(ks[i1]), $(ks[i2])] on thread $(threadid()) to $i total points")

        # check if interpolation is close enough for all sources
        # (equivalent to finding the source grid of each source separately)
        Sint = (Ss[:, :, i1] .+ Ss[:, :, i2]) ./ 2 # linear interpolation
        if refine && !all(isapprox(Ss[iS, iτs, i], Sint[iS, iτs]; kwargs...) for iS in 1:size(Ss, 1))
            atomic_add!(counter, +2)
            put!(queue, (i, i2))
            put!(queue, (i1, i))
        end

        atomic_add!(counter, -1) == 1 && put!(queue, (-1, -1)) # kill with sentinel value if counter reaches 0 (atomic_add! returns old value)
        end thread
        thread && push!(tasks, task)
    end

    thread && foreach(wait, tasks) # all tasks must finish
    close(queue) # clean up

    # sort according to k
    ks = ks[1:idx[]]
    if sort
        is = sortperm(ks)
        ks = ks[is]
        Ss = Ss[:, :, is]
    else
        Ss = Ss[:, :, 1:length(ks)]
    end

    return ks, Ss
end

# Dispatch without background solution
function source_grid_adaptive(prob::CosmologyProblem, Ss::AbstractVector, τs, ks; bgopts = (), kwargs...)
    bgsol = solvebg(prob.bg; bgopts...)
    return source_grid_adaptive(prob, Ss, τs, ks, bgsol; kwargs...)
end
