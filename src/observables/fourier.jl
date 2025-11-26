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
    return spectrum_primordial(k, sol[M.g.h], sol[M.I.As], sol[M.I.ns])
end
function spectrum_primordial(k, M::System, pars::Dict)
    return spectrum_primordial(k, pars[M.g.h], pars[M.I.As], pars[M.I.ns])
end
function spectrum_primordial(k, prob::CosmologyProblem)
    M = prob.M
    return spectrum_primordial(k, prob.bg.ps[M.g.h], prob.bg.ps[M.I.As], prob.bg.ps[M.I.ns])
end

"""
    spectrum_matter(sol::CosmologySolution, k, τ = sol[τ][end]; species = [:c, :b, :h])

Compute the power spectrum
```math
P(k,τ) = P₀(k) |Δ(k,τ)|²
```
of the total gauge-invariant overdensity
```math
Δ = δ + (3ℋ/k²) θ = (∑ₛδρ)/(∑ₛρₛ) + (3ℋ/k²) (∑ₛ(ρₛ+Pₛ)θₛ) / (∑ₛ(ρₛ+Pₛ))
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
    Δ = δ + 3*M.g.ℋ*θ/M.k^2 # total gauge-independent overdensity

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
    P = P0 .* sol(Δ^2, τ, k) # Baumann (4.4.172)

    return P
end

"""
    spectrum_matter(prob::CosmologyProblem, k::AbstractVector, τ = nothing; species = [:c, :b, :h], kwargs...)

Solve the problem `prob` with exact wavenumber(s) `k`, and then compute the power spectrum with the solution `sol`.
"""
function spectrum_matter(prob::CosmologyProblem, k::AbstractVector, τ = nothing; species = [:c, :b, :h], kwargs...)
    # save only the necessary time(s)
    if isnothing(τ)
        ptextraopts = (save_everystep = false, save_start = false, save_end = true)
    else
        ptextraopts = (saveat = [τ],)
    end
    sol = solve(prob, k; ptextraopts, kwargs...)
    τ = isnothing(τ) ? sol.bg.t[end] : τ
    return spectrum_matter(sol, k, τ; species)
end

"""
    spectrum_matter(prob::CosmologyProblem, k::NTuple{2, Number}, τs = nothing; species = [:c, :b, :h], atol = 4.0, rtol = 4e-3, coarse_length = 7, kwargs...)

Compute the matter power spectrum on the interval ``k`` with adaptively chosen wavenumbers.
Returns wavenumbers and power spectrum values.

The interval is first divided into a grid with `coarse_length` logarithmically spaced wavenumbers.
It is then adaptively refined with tolerances `atol` and `rtol`.
"""
function spectrum_matter(prob::CosmologyProblem, k::NTuple{2, Number}, τs = nothing; species = [:c, :b, :h], atol = 4.0, rtol = 4e-3, coarse_length = 7, kwargs...)
    M = prob.M
    species = getproperty.(M, filter(s -> have(M, s), species))
    δρ = sum(s.ρ*s.δ for s in species)
    ρ = sum(s.ρ for s in species)
    δ = δρ / ρ # total gauge-dependent overdensity
    ρ_plus_P_θ = sum((1+s.w)*s.ρ*s.θ for s in species) # this is the additive perturbation of the energy-momentum tensor
    ρ_plus_P = sum((1+s.w)*s.ρ for s in species) # additive # TODO: meaningful varname for this property?
    θ = ρ_plus_P_θ / ρ_plus_P
    Δ = δ + 3*M.g.ℋ*θ/M.k^2 # total gauge-independent overdensity

    # Initial coarse k-grid
    coarse_length ≥ 2 && isinteger(coarse_length) || error("Coarse must be an integer greater than or equal to 2")
    kmin, kmax = k
    ks = exp.(range(log(kmin), log(kmax), length = coarse_length))
    ks[begin] = kmin # exp(log(k)) ≠ k with floats; ensure ends are exactly what the user passed
    ks[end] = kmax # exp(log(k)) ≠ k with floats; ensure ends are exactly what the user passed

    ktransform = (log, exp)
    ks, Δs = source_grid_adaptive(prob, [Δ], τs, ks; ktransform, atol, rtol, kwargs...)
    Δs = Δs[1, 1, :]
    P0s = spectrum_primordial(ks, prob)
    Ps = P0s .* Δs .^ 2
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
# TODO: source_grid(prob::CosmologyProblem seemed type-stable?
# TODO: return getter for (ks, τs)
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
            Ss[:, :, ik] .= permutedims(getSs[iS](sol.pts[ik])) # TODO: type infer getS?
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
    solvept(prob.pt, bgsol, ks, prob.bgspline; output_func, saveat = τs, ptopts..., thread, verbose)
    return Ss
end

# TODO: Hermite interpolation
# TODO: create SourceFunction type that does k and τ interpolation?
"""
    source_grid_adaptive(prob::CosmologyProblem, Ss::AbstractVector, τs, ks; bgopts = (), ptopts = (), ktransform = (identity, identity), sort = true, thread = true, verbose = false, kwargs...)

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
function source_grid_adaptive(prob::CosmologyProblem, Ss::AbstractVector, τs, ks; bgopts = (), ptopts = (), ktransform = (identity, identity), sort = true, thread = true, verbose = false, kwargs...)
    if isnothing(τs)
        Nτs = 1
        ptsaveopts = (save_everystep = false, save_start = false, save_end = true)
    else
        Nτs = length(τs)
        ptsaveopts = (saveat = τs,)
    end

    bgsol = solvebg(prob.bg; bgopts...)
    ptprob0, ptprobgen = setuppt(prob.pt, bgsol, prob.bgspline)

    getSs = map(S -> getsym(prob.pt, S), Ss)
    function sourcek!(k, ik, Ss)
        ptprob = ptprobgen(ptprob0, k)
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

    length(ks) ≥ 2 || error("Initial k-grid must have at least 2 values")
    ninitks = length(ks)
    Ss = similar(bgsol, length(Ss), Nτs, 1024)
    ks = resize!(collect(ks), 1024)
    iτs = 1 : max(size(Ss, 2) - 1, 1) # exclude χ=0 from refinement comparison (where some CMB sources with 1/χ diverge); except if it is the only point (e.g. matter power spectrum is well-defined)

    @sync begin
    queue = Channel{Tuple{Int, Int}}(1024)
    for i in 1:ninitks
        @spawnif begin
        sourcek!(ks[i], i, Ss)
        i ≥ 2 && put!(queue, (i-1, i))
        verbose && println("Solved k = $(ks[i]) on thread $(threadid()) to $i total points")
        end thread
    end

    for (i1, i2) in queue # equivalent to "while true" with "try take!(queue) catch break end"
        @spawnif begin
        i = atomic_add!(idx, +1) + 1 # atomic_add! returns old value of idx
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
        if !all(isapprox(Ss[iS, iτs, i], Sint[iS, iτs]; kwargs...) for iS in 1:size(Ss, 1))
            atomic_add!(counter, +2)
            put!(queue, (i, i2))
            put!(queue, (i1, i))
        end

        atomic_add!(counter, -1) # finished processing

        counter[] == 0 && close(queue) # close channel when all tasks are done
        end thread
    end
    end

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
