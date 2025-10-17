using DataInterpolations
using MatterPower
using TwoFAST

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
function source_grid(sol::CosmologySolution, Ss::AbstractVector, τs)
    # Evaluate integrated perturbations on coarse grid
    ks = sol.ks
    getSs = map(S -> getsym(sol.prob.pt, S), Ss)
    Ss = similar(sol.bg, length(Ss), length(τs), length(ks))
    @tasks for ik in eachindex(ks)
        for iS in eachindex(getSs)
            Ss[:, :, ik] .= permutedims(getSs[iS](sol.pts[ik])) # TODO: type infer getS?
        end
    end
    return Ss
end

function source_grid(Ss_coarse::AbstractArray, ks_coarse, ks_fine; ktransform = identity)
    size_coarse = size(Ss_coarse)
    size_fine = (size_coarse[1], size_coarse[2], length(ks_fine))
    Ns, Nτ, _ = size(Ss_coarse)

    Ss_fine = similar(Ss_coarse, size_fine)
    xs_coarse = ktransform.(ks_coarse) # TODO: user should just pass different ks as input instead
    xs_fine = ktransform.(ks_fine)
    @tasks for iτ in 1:Nτ
        for iS in 1:Ns
            interp = LinearInterpolation(@view(Ss_coarse[iS, iτ, :]), xs_coarse)
            Ss_fine[iS, iτ, :] .= interp.(xs_fine)
        end
    end
    return Ss_fine
end

# TODO: take in kτ0s and xs
function source_grid(prob::CosmologyProblem, S::AbstractArray, τs, ks; bgopts = (), ptopts = (), thread = true, verbose = false)
    bgsol = solvebg(prob.bg; bgopts..., verbose)
    getSs = map(s -> getsym(prob.pt, s), S)
    Ss = similar(bgsol, length(S), length(τs), length(ks))
    extrema(τs) == extrema(bgsol.t) || error("input τs and computed background solution have different timespans") # TODO: don't rely on
    function output_func(sol, ik)
        for iS in eachindex(getSs)
            Ss[iS, :, ik] .= getSs[iS](sol)
        end
        return nothing, false
    end
    solvept(prob.pt, bgsol, ks, prob.bgspline; output_func, saveat = τs, ptopts..., thread, verbose)
    return Ss
end

# TODO: Hermite interpolation
# TODO: create SourceFunction type that does k and τ interpolation?
function source_grid_adaptive(prob::CosmologyProblem, Ss::AbstractVector, τs, kmin, kmax; bgopts = (), ptopts = (), verbose = false, kwargs...)
    bgsol = solvebg(prob.bg; bgopts...)

    getSs = map(S -> getsym(prob.pt, S), Ss)
    function Sk(k)
        ptsols = solvept(prob.pt, bgsol, [k], prob.bgspline; saveat = τs, ptopts...)
        ptsol = only(ptsols)
        S = transpose(stack(map(getS -> getS(ptsol), getSs)))
        S[:, end] .= 0.0 # avoid NaN/Infs from 1/χ today; fine in LOS integration, since always weighted by 0-valued Bessel function # TODO: avoid
        return S
    end

    Sint = similar(bgsol, length(Ss), length(τs))
    ks = zeros(typeof(kmin), 1024)
    Ss = similar(bgsol, length(Ss), length(τs), length(ks))
    ks[1] = kmin
    ks[2] = kmax
    Ss[:, :, 1] .= Sk(ks[1])
    Ss[:, :, 2] .= Sk(ks[2])

    i = 2
    i1i2s = [(j, j+1) for j in 1:i-1]
    while !isempty(i1i2s)
        i1, i2 = pop!(i1i2s)

        k1 = ks[i1]
        k2 = ks[i2]
        S1 = @view Ss[:, :, i1]
        S2 = @view Ss[:, :, i2]

        verbose && println("Refining k-grid between [$k1, $k2]")

        i += 1
        i > length(ks) && error("fuck")

        k = (k1 + k2) / 2
        Ss[:, :, i] .= Sk(k)
        ks[i] = k
        S = @view Ss[:, :, i]

        Sint .= (S1 .+ S2) ./ 2 # linear interpolation

        # check if interpolation is close enough for all sources
        # (equivalent to finding the source grid of each source separately)
        if !all(isapprox(S[iS, :], Sint[iS, :]; kwargs...) for iS in eachindex(getSs))
            push!(i1i2s, (i, i2), (i1, i)) # refine left and right subintervals
        end
    end

    # sort according to k
    ks = ks[1:i]
    is = sortperm(ks)
    ks = ks[is]
    Ss = Ss[:, :, is]

    return ks, Ss
end
