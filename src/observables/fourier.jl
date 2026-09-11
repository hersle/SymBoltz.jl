using DataInterpolations
using MatterPower
using TwoFAST

"""
    spectrum_primordial(k, h, As, ns=1.0; kp = 0.05/(k0*h))

Compute the primordial power spectrum
```math
P₀(k) = 2π² Aₛ (k/kₚ)^{nₛ-1} / k³
```
with spectral amplitude `As`, spectral index `ns` and pivot scale wavenumber `kp` at the wavenumber(s) `k`.
All wavenumbers are in units of ``H₀/c``, and the default pivot scale is 0.05/Mpc.
"""
function spectrum_primordial(k, h, As, ns=1.0; kp = 0.05/(k0*h)) # 0.05/Mpc in units of H₀/c
    P = 2*π^2 * As ./ k.^3
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
    return spectrum_primordial(k, prob.th.ps[M.g.h], prob.th.ps[M.I.As], prob.th.ps[M.I.ns])
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
- `k` must be a vector of wavenumbers in units of ``H₀/c``.
- `τ` must be a single or a vector of conformal times, or unspecified to use ``τ = τ₀`` today.
- `kwargs...` are keyword arguments that are forwarded to `solve(prob, k; kwargs...)`.
"""
function spectrum_matter(modes::AbstractVector, prob::CosmologyProblem, k, τ::AbstractVector; kwargs...)
    ptextraopts = (saveat = τ,)
    sol = solve(prob, k; ptextraopts, kwargs...)
    return spectrum_matter(modes, sol, k, τ)
end
function spectrum_matter(modes::AbstractVector, prob::CosmologyProblem, k; kwargs...)
    ptextraopts = (save_everystep = false, save_start = false, save_end = true)
    sol = solve(prob, k; ptextraopts, kwargs...)
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
spectrum_matter(modes::AbstractVector, sol::CosmologySolution, k; kwargs...) = spectrum_matter(modes, sol, k, maximum(sol.th.t); kwargs...) # fallback without time (today)
spectrum_matter(modes::AbstractVector, probsol, k, τ::Number; kwargs...) = spectrum_matter(modes, probsol, k, [τ])[:, 1, :] # fallback with single time
spectrum_matter(mode::Symbol, probsol, args...; kwargs...) = selectdim(spectrum_matter([mode], probsol, args...; kwargs...), 1, 1) # fallback with single mode specified
spectrum_matter(probsol::Union{CosmologyProblem, CosmologySolution}, args...; kwargs...) = spectrum_matter(:m, probsol, args...; kwargs...) # fallback with modes unspecified

"""
    spectrum_matter_nonlinear(sol::CosmologySolution, k)

Compute the nonlinear matter power spectrum from the cosmology solution `sol` at wavenumber(s) `k` using halofit implemented in MatterPower.jl.
"""
function spectrum_matter_nonlinear(sol::CosmologySolution, k)
    P = spectrum_matter(sol, k)
    M = sol.prob.M
    h = sol[M.g.h] # halofit searches for the nonlinear scale in Mpc, so convert to 1/Mpc and Mpc³ with H₀/c = k0*h/Mpc, and back
    lgPspl = spline(log.(P ./ (k0*h)^3), log.(k .* (k0*h)))
    Pf(k) = exp(lgPspl(log(k)))
    halofit_params = setup_halofit(Pf)
    Ωm0 = sol[M.m.Ω₀]
    Pf_halofit(k) = MatterPower.halofit(Pf, halofit_params, Ωm0, k)
    return Pf_halofit.(k .* (k0*h)) .* (k0*h)^3
end

# TODO: generalize to arbitrary field?
"""
    variance_matter(sol::CosmologySolution, R)

Compute the variance ``⟨δ²⟩`` of the *linear* matter density field with a top-hat filter with radius `R` in units of ``c/H₀``.
Wraps the implementation in MatterPower.jl.
"""
function variance_matter(sol::CosmologySolution, R)
    M = sol.prob.M
    k = sol.ks
    P = spectrum_matter(sol, k)
    lgPspl = spline(log.(P), log.(k))
    Pf(k) = exp(lgPspl(log(k)))
    return MatterPower.sigma2(Pf, R)
end
"""
    stddev_matter(sol::CosmologySolution, R)

Compute the standard deviation ``√(⟨δ²⟩)`` of the *linear* matter density field with a top-hat filter with radius `R` in units of ``c/H₀``.
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
    source_grid(prob::CosmologyProblem, Ss, τs, ks[, [bgsol, ]thsol]; bgopts = (), thopts = (), ptopts = (), thread = true, verbose = false)

Compute and evaluate source functions ``S(τ,k)`` with symbolic expressions `Ss` on a grid with conformal times `τs` and wavenumbers `ks` from the problem `prob`.
Returns a matrix of size `(Nτ, Nk)`, where each element is a vector of length `NS = length(Ss)` holding all source values at that `(τ, k)` point.

The options `bgopts`/`thopts`/`ptopts` are passed to the `bg`/`th`/`pt` ODE solves.
"""
function source_grid(prob::CosmologyProblem, Ss, τs, ks, bgsol::Union{Nothing, ODESolution}, thsol::ODESolution; ptopts = (), thread = true, verbose = false)
    getSs = getsym(prob.pt, Ss)
    T = source_eltype(Ss, eltype(thsol))
    minimum(τs) ≥ minimum(thsol.t) && maximum(τs) ≤ maximum(thsol.t) || error("input τs and computed background solution have different timespans")

    # Save only the requested source values instead of the full ODE solution with all unknowns
    # Save callback similar to https://github.com/SciML/SciMLBase.jl/blob/97f6d4aff88ab5f2dedc90ef503edabe72f00e93/src/solutions/ode_solutions.jl#L369-L373
    save_func(u, t, integrator) = getSs(SciMLBase.ProblemState(; u, p = SciMLBase.parameter_values(integrator), t))
    savedvalues = [SavedValues(eltype(τs), T) for _ in ks] # one save container per k
    callback(ik) = SavingCallback(save_func, savedvalues[ik]; saveat = τs)
    solvept(prob.pt, bgsol, thsol, ks; callback, save_everystep = false, save_start = false, dense = false, ptopts..., thread, verbose)

    out = Matrix{T}(undef, length(τs), length(ks))
    @inbounds for ik in eachindex(ks), iτ in eachindex(τs)
        out[iτ, ik] = savedvalues[ik].saveval[iτ]
    end
    return out
end
source_grid(prob::CosmologyProblem, Ss, τs, ks, thsol::ODESolution; kwargs...) = source_grid(prob, Ss, τs, ks, nothing, thsol; kwargs...)
function source_grid(prob::CosmologyProblem, Ss, τs, ks; bgopts = (), thopts = (), verbose = false, kwargs...)
    bgsol = solvebg(prob; verbose, bgopts...)
    thsol = solveth(prob, bgsol; verbose, thopts...)
    return source_grid(prob, Ss, τs, ks, bgsol, thsol; verbose, kwargs...)
end

function source_grid(prob::CosmologyProblem, Ss, τs, ks, kinterp::AbstractInterpolator, args...; thread = true, kwargs...)
    Ss = source_grid(prob, Ss, τs, kinterp.xs, args...; thread, kwargs...) # solve perturbation on coarse k-grid
    Ss = source_kinterp(Ss, kinterp, ks; thread) # interpolate to requested fine k-grid
    return Ss
end

function source_kinterp!(out::AbstractVector, Ss_coarse::AbstractVector, kinterp::AbstractInterpolator, ys_fine)
    kinterp(out, Ss_coarse, ys_fine)
    return out
end

# Special dispatch for returning a vector of interpolation objects (for testing)
function source_grid_interp(prob::CosmologyProblem, S, τs, kinterp::ChebyshevInterpolator, args...; kwargs...)
    Ss = source_grid(prob, S, τs, kinterp.xs, args...; kwargs...)
    ymin, ymax = kinterp.ys[end], kinterp.ys[begin]
    return [chebinterp(Ss[i, :], ymin, ymax) for i in eachindex(τs)]
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
