using NumericalIntegration
using Bessels: besselj!, sphericalbesselj
using DataInterpolations
using TwoFAST
using MatterPower

"""
    spectrum_primordial(k, h, As, ns=1.0)

Compute the primordial power spectrum with amplitude `As` at the wavenumber(s) `k`.
"""
function spectrum_primordial(k, h, As, ns=1.0)
    P = @. 2*π^2 / k^3 * As

    # compute k / kpivot with equal wavenumber units
    k = k_dimensionless.(k, h)
    kpivot = k_dimensionless(0.05 / u"Mpc", h)
    @. P *= (k/kpivot)^(ns-1)
    return P
end
function spectrum_primordial(k, sol::CosmologySolution)
    M = sol.prob.M
    return spectrum_primordial(k, sol[M.g.h], sol[M.I.As], sol[M.I.ns])
end
function spectrum_primordial(k, M::ODESystem, pars::Dict)
    return spectrum_primordial(k, pars[M.g.h], pars[M.I.As], pars[M.I.ns])
end

"""
    spectrum_matter(sol::CosmologySolution, k[, t])

Compute the matter power spectrum from the cosmology solution `sol` at wavenumber(s) `k` and conformal time(s) `t` (final, if omitted).
"""
function spectrum_matter(sol::CosmologySolution, k, t = sol[t][end]; species = [:c, :b, :h])
    M = sol.prob.M
    species = getproperty.(M, filter(s -> have(M, s), species))
    ρm = sum(s.ρ for s in species)
    Δm = M.k^2*M.g.Φ / (4π*M.g.a^2*ρm) # TODO: compute sum(s.δ*s.ρ for s in species) / sum(s.ρ for s in species) + 3*M.g.ℰ*θm/k^2, like in https://github.com/lesgourg/class_public/blob/22b49c0af22458a1d8fdf0dd85b5f0840202551b/source/perturbations.c#L6615

    # convert P (through P0) to same units as 1/k^3
    #=
    P0 = sol(k, t, M.I.P)
    kunit = only(unique(unit.(k)))
    if kunit != NoUnits
        H₀ = sol[M.g.h] * H100
        P0 *= (H₀/c / u"m")^(-3)
        Punit = unit(1 / kunit^3)
        P0 = uconvert.(Punit, P0)
    end
    =#

    P0 = spectrum_primordial(k, sol)
    P = P0 .* sol(k, t, Δm^2) # Baumann (4.4.172)

    return P
end

"""
    spectrum_matter(M::ODESystem, pars, k[, t]; species = [:c, :b, :h], solver = KenCarp4(), kwargs...)

Compute the matter power spectrum from the cosmological model `M` with parameter `pars` at wavenumber(s) `k` and conformal time(s) `t` (final, of omitted).
The `solver` and other `kwargs` are passed to `solve`.
"""
function spectrum_matter(prob::CosmologyProblem, k, t = nothing; species = [:c, :b, :h], kwargs...)
    sol = solve(prob, k; kwargs...) # TODO: just save endpoints
    t = isnothing(t) ? sol[prob.M.t][end] : t
    return spectrum_matter(sol, k, t; species)
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

# TODO: create SphericalBesselFunctionMachine-like struct, which can calculate value, derivative, ...
# TODO: contribute back to Bessels.jl
#sphericalbesseljslow(ls::AbstractArray, x) = sphericalbesselj.(ls, x)
#sphericalbesseljfast(ls::AbstractRange, x) = (x == 0.0 ? 1.0 : √(π/(2*x))) * besselj(ls .+ 0.5, x)
function sphericalbesseljfast!(out, ls::AbstractRange, x::Number)
    besselj!(out, ls .+ 0.5, x)
    if x == 0.0 && ls[begin] == 0
        out[begin] = 1.0
    elseif x != 0.0
        @. out *= √(π/(2*x))
    end
    return out
end
function sphericalbesseljslow!(out, ls::AbstractRange, x::Number)
    out .= sphericalbesselj.(ls, x)
end
function sphericalbesseljprime(l, ls, Jls)
    i = 1 + l - ls[begin] # ls[i] == l (assuming stepo of ls is 1)
    return l/(2*l+1)*Jls[i-1] - (l+1)/(2*l+1)*Jls[i+1] # analytical result (see e.g. https://arxiv.org/pdf/astro-ph/9702170 eq. (13)-(15))
end

# TODO: line-of-sight integrate Θl using ODE for evolution of Jl?
# TODO: spline sphericalbesselj for each l, from x=0 to x=kmax*(t0-tini)
# TODO: integrate with ApproxFun? see e.g. https://discourse.julialang.org/t/evaluate-integral-on-many-points-cubature-jl/1723/2
# TODO: RombergEven() works with 513 or 1025 points (do Logging.disable_logging(Logging.Warn) first)
# TODO: gaussian quadrature with weight function? https://juliamath.github.io/QuadGK.jl/stable/weighted-gauss/
# line of sight integration
# TODO: take in symbolic expr?
"""
    los_integrate(S0s::AbstractArray{T}, S1s::AbstractArray{T}, ls::AbstractArray, ks::AbstractRange, ts::AbstractArray, us::AbstractRange, u′s::AbstractArray; integrator = SimpsonEven(), verbose = false) where {T <: Real}

Compute the line-of-sight-integrals ``∫dt S(k,t) jₗ(k(t₀-t)) = ∫dt S₀(k,t) jₗ(k(t₀-t)) + ∫dt S₁(k,t) jₗ′(k(t₀-t))`` over the source function values `S0s` and `S1s` against the spherical kind-1 Bessel functions `jₗ(x)` and their derivatives `jₗ′(x)` for the given `ks` and `ls`.
The element `S0s[i,j]` holds the source function value ``S₀(ks[i], ts[j])`` (and similarly for `S1s`).
An integral substitution `u(t)` can be specified with `us` and `u′s`, so the integral can be performed as ``∫dt f(t) = ∫du f(t(u)) / u′(t)`` on an interval on which the integrand behaves well (e.g. to sample more points closer to the initial time).
"""
function los_integrate(S0s::AbstractArray{T}, S1s::AbstractArray{T}, ls::AbstractArray, ks::AbstractRange, ts::AbstractArray, us::AbstractRange, u′s::AbstractArray; integrator = SimpsonEven(), verbose = false) where {T <: Real}
    @assert size(S0s) == (length(ks), length(us)) # TODO: optimal structure? integration order? @simd?
    verbose && println("LOS integration with $(length(ls)) ls x $(length(ks)) ks x $(length(us)) us")

    lmin, lmax = extrema(ls)
    lmin >= 1 || error("l must be 1 or higher") # TODO: relax?
    lmin, lmax = lmin - 1, lmax + 1 # calculation of l-th bessel function derivatives depend on l-1 and l+1
    ls_all = lmin:1:lmax # range with step 1

    Is = zeros(T, (length(ks), length(ls)))

    @tasks for ik in eachindex(ks)
        @local begin # define task-local values (declared once for all loop iterations)
            Jls_all = zeros(Float64, length(ls_all)) # local task workspace
            ∂I_∂u = zeros(T, (length(ls), length(us))) # TODO: best to do array of arrays without using @view, or to use matrix + @view?
        end

        k = ks[ik]
        for it in eachindex(us)
            t = ts[it]
            S0 = S0s[ik,it]
            S1 = S1s[ik,it]
            kΔt = k * (ts[end]-t)
            sphericalbesseljfast!(Jls_all, ls_all, kΔt)
            for il in eachindex(ls)
                l = ls[il]
                Jl = Jls_all[1+l-lmin]
                Jl′ = sphericalbesseljprime(l, ls_all, Jls_all)
                ∂I_∂u[il,it] = (S0 * Jl + S1 * Jl′) / u′s[it] # ∫dt y(t) = ∫du y(t(u)) / u′(t(u))
            end
        end
        for il in eachindex(ls)
            integrand = @view ∂I_∂u[il,:]
            Is[ik,il] = integrate(us, integrand, integrator) # integrate over t # TODO: add starting I(tini) to fix small l?
        end
    end

    return Is
end

# this one is less elegant, but more numerically stable?
# TODO: saveat = ts
"""
    source_temperature(sol::CosmologySolution, ks::AbstractArray, ts::AbstractArray)

Compute the temperature source function ``Sᵀ(k, t)`` by interpolating in the solution object.
"""
function source_temperature(sol::CosmologySolution, ks::AbstractArray, ts::AbstractArray; sw=true, isw=true, dop=true, pol=true)
    M = sol.prob.M
    out = sol(ks, ts, [M.S0, M.S1])
    S0s, S1s = out[:, :, 1], out[:, :, 2]
    return S0s, S1s
end

"""
    source_polarization(sol::CosmologySolution, ks::AbstractArray, ts::AbstractArray)

Compute the E-mode polarization source function ``Sᴱ(k, t)`` by interpolating in the solution object.
"""
function source_polarization(sol::CosmologySolution, ks::AbstractArray, ts::AbstractArray)
    M = sol.prob.M
    t0 = sol[t][end]
    S0s = sol(ks, ts, 3/4 * M.γ.Π * M.b.rec.v) ./ (ks .* (t0 .- ts)') .^ 2 # TODO: apply integration by parts? # TODO 3/4 -> 3/16?
    S1s = 0.0 .* S0s # == 0
    return S0s, S1s
end

function los_substitution_range(sol::CosmologySolution, u::Function, u⁻¹::Function, u′::Function; kwargs...)
    tmin, tmax = extrema(sol[sol.prob.M.t])
    us = range(u(tmin), u(tmax); kwargs...)
    ts = u⁻¹.(us)
    all(u.(ts) .≈ us) || error("u(u⁻¹(t)) ≠ t")
    u′s = u′.(ts)
    return ts, us, u′s
end

"""
    los_temperature(sol::CosmologySolution, ls::AbstractArray, ks::AbstractArray; u=(t->tanh(t)), u⁻¹=(u->atanh(u)), u′=(t->1/cosh(t)^2), length=500, kwargs...)

Calculate photon temperature multipoles today by line-of-sight integration.
"""
function los_temperature(sol::CosmologySolution, ls::AbstractArray, ks::AbstractArray; u=(t->tanh(t)), u⁻¹=(u->atanh(u)), u′=(t->1/cosh(t)^2), length=500, kwargs...)
    ts, us, u′s = los_substitution_range(sol, u, u⁻¹, u′; length)
    S0s, S1s = source_temperature(sol, ks, ts)
    return los_integrate(S0s, S1s, ls, ks, ts, us, u′s; kwargs...)
end

"""
    los_polarization(sol::CosmologySolution, ls::AbstractArray, ks::AbstractArray; u=(t->tanh(t)), u⁻¹=(u->atanh(u)), u′=(t->1/cosh(t)^2), length=500, kwargs...)

Calculate photon E-mode polarization multipoles today by line-of-sight integration.
"""
function los_polarization(sol::CosmologySolution, ls::AbstractArray, ks::AbstractArray; u=(t->tanh(t)), u⁻¹=(u->atanh(u)), u′=(t->1/cosh(t)^2), length=500, kwargs...)
    ts, us, u′s = los_substitution_range(sol, u, u⁻¹, u′; length) # TODO: do tanh(t/tcmb)?
    S0s, S1s = source_polarization(sol, ks, ts)
    return los_integrate(S0s, S1s, ls, ks, ts, us, u′s; kwargs...) .* transpose(@. √((ls+2)*(ls+1)*(ls+0)*(ls-1)))
end

# TODO: integrate splines instead of trapz! https://discourse.julialang.org/t/how-to-speed-up-the-numerical-integration-with-interpolation/96223/5
"""
    spectrum_cmb(ΘlAs::AbstractArray, ΘlBs::AbstractArray, P0s::AbstractArray, ls::AbstractArray, ks::AbstractRange; integrator = SimpsonEven(), normalization = :Cl)

Compute ``Cₗᴬᴮ = (2/π) ∫dk k^2 P₀(k) Θₗᴬ(k) Θₗᴮ(k)`` for the given `ls`.
If `normaliation == :Dl`, compute ``Dₗ = Cₗ * l * (l+1) / 2π`` instead.
"""
function spectrum_cmb(ΘlAs::AbstractArray, ΘlBs::AbstractArray, P0s::AbstractArray, ls::AbstractArray, ks::AbstractRange; integrator = SimpsonEven(), normalization = :Cl)
    size(ΘlAs) == size(ΘlBs) || error("ΘlAs and ΘlBs have different sizes")
    eltype(ΘlAs) == eltype(ΘlBs) || error("ΘlAs and ΘlBs have different types")

    Cls = similar(ΘlAs, length(ls))
    ks_with0 = [0.0; ks] # add dummy value with k=0 for integration

    # Check that ks have constant spacing if integrator assumes it
    if endswith(string(integrator), "Even()")
        dks_with0 = diff(ks_with0)
        minimum(dks_with0) ≈ maximum(dks_with0) || error("ks_with0 have non-uniform spacing (min=$(minimum(dks_with0)), max=$(maximum(dks_with0)))")
    end

    @tasks for il in eachindex(ls)
        @local dCl_dks_with0 = zeros(eltype(ΘlAs), length(ks_with0)) # local task workspace
        @. dCl_dks_with0[2:end] = 2/π * ks^2 * P0s * ΘlAs[:,il] * ΘlBs[:,il]
        Cls[il] = integrate(ks_with0, dCl_dks_with0, integrator) # integrate over k (_with0 adds one additional point at (0,0))
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
    spectrum_cmb(modes, sol::CosmologySolution, ls::AbstractArray; kwargs...)

Compute the CMB power spectra `modes` (`:TT`, `:EE`, `:TE` or an array thereof) ``C_l^{AB}``'s at angular wavenumbers `ls` from the cosmological solution `sol`.
"""
function spectrum_cmb(modes::AbstractArray, sol::CosmologySolution, ls::AbstractArray; normalization = :Cl, kwargs...)
    _, ks_fine = cmb_ks(ls[end])
    ΘlTs = 'T' in join(modes) ? los_temperature(sol, ls, ks_fine; kwargs...) : nothing
    ΘlPs = 'E' in join(modes) ? los_polarization(sol, ls, ks_fine; kwargs...) : nothing
    P0s = spectrum_primordial(ks_fine, sol) # more accurate
    #P0s = sol(ks_fine, sol[t][begin], sol.M.I.P) # less accurate (requires smaller Δk_S, e.g. Δk_S = 1.0 instead of 10.0)
    spectra = [] # Cls or Dls
    for mode in modes
        mode == :TT && push!(spectra, spectrum_cmb(ΘlTs, ΘlTs, P0s, ls, ks_fine; normalization))
        mode == :EE && push!(spectra, spectrum_cmb(ΘlPs, ΘlPs, P0s, ls, ks_fine; normalization))
        mode == :TE && push!(spectra, spectrum_cmb(ΘlTs, ΘlPs, P0s, ls, ks_fine; normalization))
    end
    return spectra
end

"""
    spectrum_cmb(modes::AbstractArray, prob::CosmologyProblem, ls::AbstractArray; integrator = SimpsonEven(), kwargs...)

Compute the CMB power spectra `modes` (`:TT`, `:EE`, `:TE` or an array thereof) ``C_l^{AB}``'s at angular wavenumbers `ls` from the cosmological problem `prob`.
"""
function spectrum_cmb(modes::AbstractArray, prob::CosmologyProblem, ls::AbstractArray; integrator = SimpsonEven(), normalization = :Cl, kwargs...)
    ks_coarse, _ = cmb_ks(ls[end])
    sol = solve(prob, ks_coarse; kwargs...) # TODO: saveat ts
    return spectrum_cmb(modes, sol, ls; normalization, integrator)
end
spectrum_cmb(mode::Symbol, args...; kwargs...) = only(spectrum_cmb([mode], args...; kwargs...))

function cmb_ks(lmax; Δk = 2π/10, Δk_S = 5.0, kmin = Δk, kmax = 2*lmax)
    # Assumes t0 = 1 (e.g. t0 = 1/H0 = 1) # TODO: don't assume t0 = 1
    ks_fine = range(kmin, kmax, step=Δk) # use integer multiple so endpoints are the same
    ks_coarse = range(ks_fine[begin], ks_fine[end], length = Int(floor((kmax-kmin)/Δk_S+1)))
    ks_coarse[begin] == ks_fine[begin] && ks_coarse[end] == ks_fine[end] || error("different wavenumber endpoints")
    # kmax is a guideline and may be slightly different from ks_fine[end] and ks_coarse[end]
    return ks_coarse, ks_fine
end

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

# TODO: not really a spectrum...
"""
    distance_luminosity(sol::CosmologySolution, t = sol[sol.prob.M.t], t0 = sol[sol.prob.M.t][end])

Compute luminosity distances at the time(s) `t` relative to the (present) time `t0`.
"""
function distance_luminosity(sol::CosmologySolution, t = sol[sol.prob.M.t], t0 = sol[sol.prob.M.t][end])
    M = sol.prob.M
    χ = t0 .- t
    Ωk0 = have(M, :K) ? sol[M.K.Ω₀] : 0.0
    r = sinc.(√(-Ωk0+0im)*χ/π) .* χ |> real # Julia's sinc(x) = sin(π*x) / (π*x)
    H0 = H100 * sol[M.g.h]
    a = sol(t, M.g.a)
    return @. r / a * SymBoltz.c / H0 # to meters
end
