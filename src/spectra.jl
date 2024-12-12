using NumericalIntegration
using Bessels: besselj!, sphericalbesselj
using DataInterpolations
using ForwardDiff
using TwoFAST
using MatterPower

"""
    spectrum_primordial(k; As=2e-9)

Compute the primordial power spectrum with amplitude `As` at the wavenumber(s) `k`.
"""
spectrum_primordial(k; As=2e-9) = @. 2*π^2 / k^3 * As # TODO: add kpivot and ns # TODO: make separate InflationModel with these parameters

"""
    spectrum_matter(sol::CosmologySolution, k)

Compute the matter power spectrum from the cosmology solution `sol` at wavenumber(s) `k`.
"""
function spectrum_matter(sol::CosmologySolution, k; species = [:c, :b, :h])
    tend = sol[t][end]
    M = sol.M
    species = getproperty.(M, filter(s -> have(M.sys, s), species))
    ρm = sum(s.ρ for s in species)
    Δm = M.k^2*M.g.Φ / (4π*M.g.a^2*ρm) # TODO: compute sum(s.δ*s.ρ for s in species) / sum(s.ρ for s in species) + 3*M.g.ℰ*θm/k^2, like in https://github.com/lesgourg/class_public/blob/22b49c0af22458a1d8fdf0dd85b5f0840202551b/source/perturbations.c#L6615
    P0 = spectrum_primordial(k)
    return P0 .* sol(k, tend, [Δm])[:, 1, 1] .^ 2 # Baumann (4.4.172)
end

"""
    spectrum_matter(M::CosmologyModel, pars, k; species = [:c, :b, :h], solver = KenCarp4(), kwargs...)

Compute the matter power spectrum from the cosmological model `M` with parameter `pars` at wavenumber(s) `k`.
The `solver` and other `kwargs` are passed to `solve`.
"""
function spectrum_matter(M::CosmologyModel, pars, k; species = [:c, :b, :h], solver = KenCarp4(), kwargs...)
    sol = solve(M, pars, k; save_everystep=false, solver, kwargs...) # just save endpoints
    return spectrum_matter(sol, k; species)
end

"""
    spectrum_matter_nonlinear(sol::CosmologySolution, k)

Compute the nonlinear matter power spectrum from the cosmology solution `sol` at wavenumber(s) `k` using halofit implemented in MatterPower.jl.
"""
function spectrum_matter_nonlinear(sol::CosmologySolution, k)
    P = spectrum_matter(sol, k)
    lgPspl = CubicSpline(log.(ustrip(P)), log.(ustrip(k)); extrapolate=true)
    Pf(k) = exp(lgPspl(log(k)))
    halofit_params = setup_halofit(Pf)
    Ωm0 = sol[sol.M.c.Ω₀ + sol.M.b.Ω₀] # TODO: generalize to massive neutrinos, redshift etc.
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
    k = sol.ks
    P = spectrum_matter(sol, k)
    lgPspl = CubicSpline(log.(P), log.(k); extrapolate=true)
    Pf(k) = exp(lgPspl(log(k)))
    R = 1 / k_dimensionless(1 / R, sol[sol.M.g.h]) # make dimensionless
    return MatterPower.sigma2(Pf, R)
end
"""
    stddev_matter(sol::CosmologySolution, R)

Compute the standard deviation ``√(⟨δ²⟩)`` of the *linear* matter density field with a top-hat filter with radius `R`.
"""
stddev_matter(sol::CosmologySolution, R) = √(variance_matter(sol, R))

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

# TODO: line-of-sight integrate Θl using ODE for evolution of Jl?
# TODO: spline sphericalbesselj for each l, from x=0 to x=kmax*(t0-tini)
# TODO: integrate with ApproxFun? see e.g. https://discourse.julialang.org/t/evaluate-integral-on-many-points-cubature-jl/1723/2
# TODO: RombergEven() works with 513 or 1025 points (do Logging.disable_logging(Logging.Warn) first)
# TODO: gaussian quadrature with weight function? https://juliamath.github.io/QuadGK.jl/stable/weighted-gauss/
# line of sight integration
# TODO: take in symbolic expr?
"""
    los_integrate(Ss::AbstractArray, ls::AbstractArray, ks::AbstractRange, lnts::AbstractRange; integrator = SimpsonEven(), verbose = true)

Compute the line-of-sight-integrals ``∫dt S(k,t) jₗ(k(t₀-t))`` over the source function values `Ss` against the spherical kind-1 Bessel functions `jₗ` for the given `ks` and `ls`.
The element `Ss[i,j]` holds the source function value ``S(k[i], exp(lnts[j]))``.
"""
function los_integrate(Ss::AbstractArray, ls::AbstractArray, ks::AbstractRange, lnts::AbstractRange; t0 = exp(lnts[end]), integrator = SimpsonEven(), verbose = true)
    @assert size(Ss) == (length(ks), length(lnts)) # TODO: optimal structure?
    verbose && println("LOS integration with $(length(ls)) ls x $(length(ks)) ks x $(length(lnts)) lnts")

    lmin, lmax = extrema(ls)
    ls_all = lmin:1:lmax # range with step 1

    T = eltype(Ss)
    Is = zeros(T, (length(ks), length(ls)))

    @tasks for ik in eachindex(ks)
        @local begin # define task-local values (declared once for all loop iterations)
            Jls_all = zeros(Float64, length(ls_all)) # local task workspace
            ∂I_∂lnt = zeros(T, (length(ls), length(lnts))) # TODO: best to do array of arrays without using @view, or to use matrix + @view?
        end

        k = ks[ik]
        for it in eachindex(lnts)
            t = exp(lnts[it])
            St = Ss[ik,it] * t # multiply S(k,t) to get (S(k, log(t)))
            kΔt = k * (t0-t)
            sphericalbesseljfast!(Jls_all, ls_all, kΔt) # TODO: reuse ∂Θ_∂lnt's memory?
            for il in eachindex(ls) # TODO: @simd if I can make Jl access with unit stride? also need @inbounds?
                l = ls[il]
                Jl = Jls_all[1+l-lmin]
                ∂I_∂lnt[il,it] = St * Jl
                # TODO: integrate in this loop instead?
            end
        end
        for il in eachindex(ls)
            integrand = @view ∂I_∂lnt[il,:]
            Is[ik,il] = integrate(lnts, integrand, integrator) # integrate over t # TODO: add starting I(tini) # TODO: calculate ∂Θ_∂logΘ and use Even() methods
        end
    end

    return Is
end

# this one is less elegant, but more numerically stable?
# TODO: saveat = ts
# TODO: restore sol(ks, exp.(lnts), sol.M.S)
"""
    source_temperature(sol::CosmologySolution, ks::AbstractArray, ts::AbstractArray)

Compute the temperature source function ``Sᵀ(k, t)`` by interpolating in the solution object.
"""
function source_temperature(sol::CosmologySolution, ks::AbstractArray, ts::AbstractArray)
    M = sol.M
    Ss = zeros((length(ks), length(ts))) # TODO: change order to get DenseArray during integrations?

    τ = sol(ts, M.b.rec.τ) # TODO: assume ts[end] is today
    v = -D_spline(τ, ts) .* exp.(-τ)
    idxs = [M.γ.δ, M.g.Ψ, M.γ.Π, M.g.Φ, M.b.u]
    @tasks for ik in eachindex(ks)
        k = ks[ik]
        out = sol(k, ts, idxs)
        δ, Ψ, Π, Φ, ub = selectdim.(Ref(out), 2, eachindex(idxs))
        Ss[ik,:] .= v .* (δ/4 + Ψ + Π/4) + exp.(-τ) .* D_spline(Ψ + Φ, ts) + D_spline(v .* ub, ts) / k + 3/(4*k^2) * D_spline(v .* Π, ts; order = 2) # Dodelson (9.57) with Φ → -Φ and polarization
        #Ss[ik,:] .= exp.(-τ) .* (D_spline(Φ, ts) - τ̇/4 .* (δ + Π)) + D_spline(exp.(-τ) .* (Ψ - ub.*τ̇/k), ts) + 3/(4*k^2) * D_spline(v .* Π, ts; order = 2) # Dodelson (9.55) with Φ → -Φ
        #Ss[ik,:] .= v .* (δ/4 + Ψ + Π/4) + v .* (Φ-Ψ) + 2 * exp.(-τ) .* D_spline(Φ, ts) + D_spline(v .* ub, ts) / k + exp.(-τ) * k .* (Ψ - Φ) + 3/(4*k^2) * D_spline(v .* Π, ts; order = 2) # CLASS' expression with added polarization
    end

    return Ss
end

"""
    source_polarization(sol::CosmologySolution, ks::AbstractArray, ts::AbstractArray)

Compute the E-mode polarization source function ``Sᴱ(k, t)`` by interpolating in the solution object.
"""
function source_polarization(sol::CosmologySolution, ks::AbstractArray, ts::AbstractArray)
    M = sol.M
    t0 = sol[t][end]
    return sol(ks, ts, 3/4 * M.γ.Π * M.b.rec.v) ./ (ks .* (t0 .- ts)') .^ 2
end

# TODO: increase Δlnt. should use same Δlnt in T and E when cross-correlating. TE seems to need ≲ 0.01
"""
    los_temperature(sol::CosmologySolution, ls::AbstractArray, ks::AbstractArray; Δlnt = 0.01, kwargs...)

Calculate photon temperature multipoles today by line-of-sight integration.
"""
function los_temperature(sol::CosmologySolution, ls::AbstractArray, ks::AbstractArray; Δlnt = 0.01, kwargs...)
    tmin, tmax = extrema(sol[sol.M.t])
    lnts = range(log(tmin), log(tmax), step=Δlnt)
    STs = source_temperature(sol, ks, exp.(lnts))
    return los_integrate(STs, ls, ks, lnts; kwargs...)
end

"""
    los_polarization(sol::CosmologySolution, ls::AbstractArray, ks::AbstractArray; Δlnt = 0.01, kwargs...)

Calculate photon E-mode polarization multipoles today by line-of-sight integration.
"""
function los_polarization(sol::CosmologySolution, ls::AbstractArray, ks::AbstractArray; Δlnt = 0.01, kwargs...)
    tmin, tmax = extrema(sol[sol.M.t])
    lnts = range(log(tmin), log(tmax), step=Δlnt)
    SPs = source_polarization(sol, ks, exp.(lnts))
    return los_integrate(SPs, ls, ks, lnts; t0=sol[sol.M.t][end], kwargs...) .* transpose(@. √((ls+2)*(ls+1)*(ls+0)*(ls-1)))
end

# TODO: integrate splines instead of trapz! https://discourse.julialang.org/t/how-to-speed-up-the-numerical-integration-with-interpolation/96223/5
# TODO: take scalar l as input
"""
    spectrum_cmb(ΘlAs::AbstractArray, ΘlBs::AbstractArray, P0s::AbstractArray, ls::AbstractArray, ks::AbstractRange; integrator = SimpsonEven())

Compute ``Cₗᴬᴮ = (2/π) ∫dk k^2 P0(k) Θₗᴬ(k) Θₗᴮ(k)`` for the given `ls`.
"""
function spectrum_cmb(ΘlAs::AbstractArray, ΘlBs::AbstractArray, P0s::AbstractArray, ls::AbstractArray, ks::AbstractRange; integrator = SimpsonEven())
    size(ΘlAs) == size(ΘlBs) || error("ΘlAs and ΘlBs have different size")
    eltype(ΘlAs) == eltype(ΘlBs) || error("ΘlAs and ΘlBs have different types")

    Cls = similar(ΘlAs, length(ls))
    ks_with0 = [0.0; ks] # add dummy value with k=0 for integration

    @tasks for il in eachindex(ls)
        @local dCl_dks_with0 = zeros(eltype(ΘlAs), length(ks_with0)) # local task workspace
        @. dCl_dks_with0[2:end] = 2/π * ks^2 * P0s * ΘlAs[:,il] * ΘlBs[:,il]
        Cls[il] = integrate(ks_with0, dCl_dks_with0, integrator) # integrate over k (_with0 adds one additional point at (0,0))
    end

    return Cls
end

"""
    spectrum_cmb(modes, sol::CosmologySolution, ls::AbstractArray; kwargs...)

Compute the CMB power spectra `modes` (`:TT`, `:EE`, `:TE` or an array thereof) ``C_l^{AB}``'s at angular wavenumbers `ls` from the cosmological solution `sol`.
"""
function spectrum_cmb(modes::AbstractArray, sol::CosmologySolution, ls::AbstractArray; kwargs...)
    _, ks_fine = cmb_ks(ls[end])
    ΘlTs = 'T' in join(modes) ? los_temperature(sol, ls, ks_fine; kwargs...) : nothing
    ΘlPs = 'E' in join(modes) ? los_polarization(sol, ls, ks_fine; kwargs...) : nothing
    P0s = spectrum_primordial(ks_fine)
    Cls = []
    for mode in modes
        mode == :TT && push!(Cls, spectrum_cmb(ΘlTs, ΘlTs, P0s, ls, ks_fine))
        mode == :EE && push!(Cls, spectrum_cmb(ΘlPs, ΘlPs, P0s, ls, ks_fine))
        mode == :TE && push!(Cls, spectrum_cmb(ΘlTs, ΘlPs, P0s, ls, ks_fine))
    end
    return Cls
end

"""
    spectrum_cmb(modes, M::CosmologyModel, pars::Dict, ls::AbstractArray; integrator = SimpsonEven(), kwargs...)

Compute the CMB power spectra `modes` (`:TT`, `:EE`, `:TE` or an array thereof) ``C_l^{AB}``'s at angular wavenumbers `ls` from the cosmological model `M` with parameters `pars`.
"""
function spectrum_cmb(modes::AbstractArray, M::CosmologyModel, pars::Dict, ls::AbstractArray; integrator = SimpsonEven(), kwargs...)
    ks_coarse, _ = cmb_ks(ls[end])
    sol = solve(M, pars, ks_coarse; kwargs...) # TODO: saveat ts
    return spectrum_cmb(modes, sol, ls; integrator)
end

spectrum_cmb(mode::Symbol, args...; kwargs...) = only(spectrum_cmb([mode], args...; kwargs...))

function cmb_ks(lmax; Δk = 2π/24, Δk_S = 10.0, kmin = Δk, kmax = lmax)
    # Assumes t0 = 1 (e.g. t0 = 1/H0 = 1) # TODO: don't assume t0 = 1
    ks_coarse = range(kmin, kmax, length = Int(floor((kmax-kmin)/Δk_S+1)))
    ks_fine = range(ks_coarse[begin], ks_coarse[end], length = length(ks_coarse) * Int(floor(Δk_S/Δk))) # use integer multiple so endpoints are the same
    return ks_coarse, ks_fine
end

Dl(Cl, l) = @. Cl * l * (l+1) / 2π

"""
    correlation_function(sol::CosmologySolution; N = 2048, spline = true)

Compute the two-point correlation function in real space by Fourier transforming the matter power spectrum of `sol` with `N` points the FFTLog algorithm implemented in TwoFAST.
Returns `N` radii and correlation function values (e.g. `r`, `ξ`).
"""
function correlation_function(sol::CosmologySolution; N = 2048, spline = true)
    ks = sol.ks
    if spline
        P = CubicSpline(spectrum_matter(sol, ks), ks) # create spline interpolation (fast)
    else
        P(k) = only(spectrum_matter(sol, k)) # use solution's built-in interpolation (elegant)
    end
    kmin, kmax = extrema(ks)
    rmin = 2π / kmax
    return xicalc(P, 0, 0; N, kmin, kmax, r0=rmin)
end
