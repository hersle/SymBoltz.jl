using Integrals
using Bessels: besselj!, sphericalbesselj
using DataInterpolations
using TwoFAST
using MatterPower
using ForwardDiff # TODO: get rid of; to convert to Float64 inside spectrum_cmb

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
function spectrum_primordial(k, M::ODESystem, pars::Dict)
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
Δ = δ + (3ℰ/k²) θ = (∑ₛδρ)/(∑ₛρₛ) + (3ℰ/k²) (∑ₛ(ρₛ+Pₛ)θₛ) / (∑ₛ(ρₛ+Pₛ))
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
    Δ = δ + 3*M.g.ℰ*θ/M.k^2 # total gauge-independent overdensity

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
    P = P0 .* sol(k, τ, Δ^2) # Baumann (4.4.172)

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

# TODO: create SphericalBesselFunctionMachine-like struct, which can calculate value, derivative, ...
# TODO: contribute back to Bessels.jl
#sphericalbesseljslow(ls::AbstractVector, x) = sphericalbesselj.(ls, x)
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
    i = 1 + l - ls[begin] # ls[i] == l (assuming step of ls is 1)
    return l/(2*l+1)*Jls[i-1] - (l+1)/(2*l+1)*Jls[i+1] # analytical result (see e.g. https://arxiv.org/pdf/astro-ph/9702170 eq. (13)-(15))
end
function sphericalbesselj_over_x2!(out, ls::AbstractRange, x::Number)
    x = max(x, 1e-20) # avoid division by zero
    sphericalbesseljfast!(out, ls, x)
    out ./= x .^2
    return out
end

function integrate(xs, ys; integrator = TrapezoidalRule())
    prob = SampledIntegralProblem(ys, xs)
    sol = solve(prob, integrator)
    return sol.u
end

# TODO: line-of-sight integrate Θl using ODE for evolution of Jl?
# TODO: spline sphericalbesselj for each l, from x=0 to x=kmax*(τ0-τini)
# TODO: integrate with ApproxFun? see e.g. https://discourse.julialang.org/t/evaluate-integral-on-many-points-cubature-jl/1723/2
# TODO: RombergEven() works with 513 or 1025 points (do Logging.disable_logging(Logging.Warn) first)
# TODO: gaussian quadrature with weight function? https://juliamath.github.io/QuadGK.jl/stable/weighted-gauss/
# line of sight integration
# TODO: take in symbolic expr?
# TODO: implement chain rule using fundamental theorem of calculus, to make faster and remove need of differentiating Bessel functions? https://juliadiff.org/ChainRulesCore.jl/stable/rule_author/which_functions_need_rules.html
# TODO: use u = k*χ as integration variable, so oscillations of Bessel functions are the same for every k?
"""
    los_integrate(S0s::AbstractMatrix{T}, S1s::AbstractMatrix{T}, ls::AbstractVector, ks::AbstractVector, τs::AbstractVector, Rl0!::Function; integrator = TrapezoidalRule(), verbose = false) where {T <: Real}

Compute the line-of-sight-integrals ``∫dτ S(k,τ) jₗ(k(τ₀-τ)) = ∫dτ S₀(k,τ) jₗ(k(τ₀-τ)) + ∫dτ S₁(k,τ) jₗ′(k(τ₀-τ))`` over the source function values `S0s` and `S1s` against the spherical kind-1 Bessel functions `jₗ(x)` and their derivatives `jₗ′(x)` for the given `ks` and `ls`.
The element `S0s[i,j]` holds the source function value ``S₀(ks[i], τs[j])`` (and similarly for `S1s`).
An integral substitution `u(τ)` can be specified with `us` and `u′s`, so the integral can be performed as ``∫dτ f(τ) = ∫du f(τ(u)) / u′(τ)`` on an interval on which the integrand behaves well (e.g. to sample more points closer to the initial time).
"""
function los_integrate(S0s::AbstractMatrix{T}, S1s::AbstractMatrix{T}, ls::AbstractVector, ks::AbstractVector, τs::AbstractVector, Rl0!::Function; integrator = TrapezoidalRule(), verbose = false) where {T <: Real}
    @assert size(S0s) == (length(ks), length(τs)) # TODO: optimal structure? integration order? @simd?
    verbose && println("LOS integration with $(length(ls)) ls x $(length(ks)) ks x $(length(τs)) τs")

    # TODO: only do to total τ0 - τ? otherwise Dual errors can propagate if one is dual and the other isn't
    τs = ForwardDiff.value.(τs) # convert from Duals to e.g. Float64s; it is the independent parameter and is only used as an integration range, so should not affect the final result # TODO: safe?
    χs = τs[end] .- τs

    lmin, lmax = extrema(ls)
    lmin >= 1 || error("l must be 1 or higher") # TODO: relax?
    lmin, lmax = lmin - 1, lmax + 1 # calculation of l-th bessel function derivatives depend on l-1 and l+1
    ls_all = lmin:1:lmax # range with step 1

    Is = zeros(T, (length(ks), length(ls)))

    @tasks for ik in eachindex(ks)
        @local begin # define task-local values (declared once for all loop iterations)
            lmin = lmin
            Jls_all = zeros(Float64, length(ls_all)) # local task workspace
            ∂I_∂τ = zeros(T, (length(ls), length(τs))) # TODO: best to do array of arrays without using @view, or to use matrix + @view?
        end

        k = ks[ik]
        for iτ in eachindex(τs)
            χ = χs[iτ]
            S0 = S0s[ik,iτ]
            S1 = S1s[ik,iτ] # TODO: don't need?
            kχ = k * χ
            Rl0!(Jls_all, ls_all, kχ)
            for il in eachindex(ls)
                # TODO: skip and set to zero if l ≳ kτ0?
                l = ls[il]
                Jl = Jls_all[1+l-lmin]
                Jl′ = sphericalbesseljprime(l, ls_all, Jls_all) # TODO: only valid if Rl0 = jₗ
                ∂I_∂τ[il,iτ] = S0 * Jl + S1 * Jl′
            end
        end
        for il in eachindex(ls)
            # TODO: skip and set to zero if l ≳ kτ0?
            integrand = @view ∂I_∂τ[il,:]
            Is[ik,il] = integrate(τs, integrand; integrator) # integrate over τ # TODO: add starting I(τini) to fix small l?
        end
    end

    return Is
end

"""
    los_temperature(sol::CosmologySolution, ls::AbstractVector, ks::AbstractVector, τs::AbstractVector; ktransform = identity, kwargs...)

Calculate photon temperature multipoles today by line-of-sight integration.
"""
function los_temperature(sol::CosmologySolution, ls::AbstractVector, ks::AbstractVector, τs::AbstractVector; ktransform = identity, kwargs...)
    M = sol.prob.M
    out = sol(ks, τs, [M.ST0, M.ST1]; ktransform)
    S0s, S1s = out[:, :, 1], out[:, :, 2]
    return los_integrate(S0s, S1s, ls, ks, τs, sphericalbesseljfast!; kwargs...)
end

# TODO: make a function for calculating jl(x)/x^2 for use in polarization
"""
    los_polarization(sol::CosmologySolution, ls::AbstractVector, ks::AbstractVector, τs::AbstractVector; ktransform = identity, kwargs...)

Calculate photon E-mode polarization multipoles today by line-of-sight integration.
"""
function los_polarization(sol::CosmologySolution, ls::AbstractVector, ks::AbstractVector, τs::AbstractVector; ktransform = identity, kwargs...)
    M = sol.prob.M
    S0s = sol(ks, τs, 3/16*M.γ.Π*M.b.rec.v; ktransform) # TODO: apply integration by parts?
    S1s = 0.0 .* S0s # == 0
    return los_integrate(S0s, S1s, ls, ks, τs, sphericalbesselj_over_x2!; kwargs...) .* transpose(@. √((ls+2)*(ls+1)*(ls+0)*(ls-1)))
end

# TODO: integrate splines instead of trapz! https://discourse.julialang.org/t/how-to-speed-up-the-numerical-integration-with-interpolation/96223/5
raw"""
    spectrum_cmb(ΘlAs::AbstractMatrix, ΘlBs::AbstractMatrix, P0s::AbstractVector, ls::AbstractVector, ks::AbstractVector; integrator = TrapezoidalRule(), normalization = :Cl)

Compute the angular power spectrum
```math
Cₗᴬᴮ = (2/π) ∫\mathrm{d}k \, k² P₀(k) Θₗᴬ(k,τ₀) Θₗᴮ(k,τ₀)
```
for the given `ls`.
If `normaliation == :Dl`, compute ``Dₗ = Cₗ l (l+1) / 2π`` instead.
"""
function spectrum_cmb(ΘlAs::AbstractMatrix, ΘlBs::AbstractMatrix, P0s::AbstractVector, ls::AbstractVector, ks::AbstractVector; integrator = TrapezoidalRule(), normalization = :Cl)
    size(ΘlAs) == size(ΘlBs) || error("ΘlAs and ΘlBs have different sizes")
    eltype(ΘlAs) == eltype(ΘlBs) || error("ΘlAs and ΘlBs have different types")

    Cls = similar(ΘlAs, length(ls))
    ks_with0 = [0.0; ks] # add dummy value with k=0 for integration

    @tasks for il in eachindex(ls)
        # TODO: skip kτ0 ≲ l
        @local dCl_dks_with0 = zeros(eltype(ΘlAs), length(ks_with0)) # local task workspace
        @. dCl_dks_with0[2:end] = 2/π * ks^2 * P0s * ΘlAs[:,il] * ΘlBs[:,il]
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
    spectrum_cmb(modes::AbstractVector, sol::CosmologySolution, ls::AbstractVector, ks::AbstractVector, τs::AbstractVector; normalization = :Cl, unit = nothing, integrator = TrapezoidalRule(), kwargs...)

Compute the CMB power spectra `modes` (`:TT`, `:EE`, `:TE` or an array thereof) ``C_l^{AB}``'s at angular wavenumbers `ls` from the cosmological solution `sol`.
If `unit` is `nothing` the spectra are of dimensionless temperature fluctuations relative to the present photon temperature; while if `unit` is a temperature unit the spectra are of dimensionful temperature fluctuations.
"""
function spectrum_cmb(modes::AbstractVector, sol::CosmologySolution, ls::AbstractVector, ks::AbstractVector, τs::AbstractVector; normalization = :Cl, unit = nothing, integrator = TrapezoidalRule(), kwargs...)
    ΘlTs = 'T' in join(modes) ? los_temperature(sol, ls, ks, τs; integrator, kwargs...) : nothing
    ΘlPs = 'E' in join(modes) ? los_polarization(sol, ls, ks, τs; integrator, kwargs...) : nothing
    P0s = spectrum_primordial(ks, sol) # more accurate
    #P0s = sol(ks, sol[τ][begin], sol.M.I.P) # less accurate (requires smaller Δk_S, e.g. Δk_S = 1.0 instead of 10.0)

    if isnothing(unit)
        factor = 1 # keep dimensionless
    elseif dimension(unit) == dimension(u"K")
        factor = uconvert(unit, sol[sol.prob.M.γ.T₀] * u"K") # convert to a temperature unit
    else
        error("Requested unit $unit is not a temperature unit")
    end

    spectra = [] # Cls or Dls
    for mode in modes
        if mode == :TT
            spectrum = spectrum_cmb(ΘlTs, ΘlTs, P0s, ls, ks; integrator, normalization)
        elseif mode == :EE
            spectrum = spectrum_cmb(ΘlPs, ΘlPs, P0s, ls, ks; integrator, normalization)
        elseif mode == :TE
            spectrum = spectrum_cmb(ΘlTs, ΘlPs, P0s, ls, ks; integrator, normalization)
        else
            error("Unknown CMB power spectrum mode $mode")
        end
        spectrum *= factor^2 # possibly make dimensionful
        push!(spectra, spectrum)
    end

    return spectra
end

"""
    spectrum_cmb(modes::AbstractVector, prob::CosmologyProblem, ls::AbstractVector; normalization = :Cl, unit = nothing, Δkτ0 = 2π/2, Δkτ0_S = 8.0, kτ0min = 0.1*ls[begin], kτ0max = 3*ls[end], u = (τ->tanh(τ)), u⁻¹ = (u->atanh(u)), Nlos = 768, integrator = TrapezoidalRule(), bgopts = (alg = Rodas4P(), reltol = 1e-8), ptopts = (alg = KenCarp4(), reltol = 1e-8), kwargs...)

Compute the CMB power spectra `modes` (`:TT`, `:EE`, `:TE` or an array thereof) ``C_l^{AB}``'s at angular wavenumbers `ls` from the cosmological problem `prob`.
"""
function spectrum_cmb(modes::AbstractVector, prob::CosmologyProblem, ls::AbstractVector; normalization = :Cl, unit = nothing, Δkτ0 = 2π/2, Δkτ0_S = 8.0, kτ0min = 0.1*ls[begin], kτ0max = 3*ls[end], u = (τ->tanh(τ)), u⁻¹ = (u->atanh(u)), Nlos = 768, integrator = TrapezoidalRule(), bgopts = (alg = Rodas4P(), reltol = 1e-8), ptopts = (alg = KenCarp4(), reltol = 1e-8), kwargs...)
    kτ0s_coarse, kτ0s_fine = cmb_kτ0s(ls[begin], ls[end]; Δkτ0, Δkτ0_S, kτ0min, kτ0max)
    sol = solve(prob; bgopts)
    τ0 = getp(sol, prob.M.τ0)(sol) |> ForwardDiff.value
    ks_coarse = kτ0s_coarse / τ0
    τs = sol.bg.t # by default, use background (thermodynamics) time points for line of sight integration
    if Nlos != 0 # instead choose Nlos time points τ = τ(u) corresponding to uniformly spaced u
        τmin, τmax = extrema(τs)
        umin, umax = u(τmin), u(τmax)
        us = range(umin, umax, length = Nlos)
        τs = u⁻¹.(us)
    end
    #ptopts = merge(ptopts, (saveat = τs,)) # TODO: use, but would need to redo k integration
    sol = solve(prob, ks_coarse; ptopts)
    ks_fine = kτ0s_fine / τ0
    return spectrum_cmb(modes, sol, ls, ks_fine, τs; normalization, unit, integrator, kwargs...)
end
spectrum_cmb(mode::Symbol, args...; kwargs...) = only(spectrum_cmb([mode], args...; kwargs...))

function cmb_kτ0s(lmin, lmax; Δkτ0 = 2π/2, Δkτ0_S = 8.0, kτ0min = 0.1*lmin, kτ0max = 3*lmax)
    kτ0s_fine = range(kτ0min, kτ0max, step=Δkτ0) # use integer multiple so endpoints are the same
    kτ0s_coarse = range(kτ0s_fine[begin], kτ0s_fine[end], length = Int(floor((kτ0max-kτ0min)/Δkτ0_S+1)))
    kτ0s_coarse[begin] == kτ0s_fine[begin] && kτ0s_coarse[end] == kτ0s_fine[end] || error("different wavenumber endpoints")
    return kτ0s_coarse, kτ0s_fine
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
# TODO: add formula
"""
    distance_luminosity(sol::CosmologySolution, ivs = sol.bg.t, τ0 = sol[sol.prob.M.τ0])

Compute luminosity distances at the independent variable values `ivs` relative to the (present) time `τ0`.
"""
function distance_luminosity(sol::CosmologySolution, ivs = sol.bg.t, τ0 = sol[sol.prob.M.τ0])
    M = sol.prob.M
    χ = sol(ivs, M.χ)
    Ωk0 = have(M, :K) ? sol[M.K.Ω₀] : 0.0
    r = sinc.(√(-Ωk0+0im)*χ/π) .* χ |> real # Julia's sinc(x) = sin(π*x) / (π*x)
    H0 = H100 * sol[M.g.h]
    a = sol(ivs, M.g.a)
    return @. r / a * SymBoltz.c / H0 # to meters
end

# TODO: test @inferred
function distance_luminosity_function(M::ODESystem, pars_fixed, pars_varying, zs; bgopts = (alg = SymBoltz.Tsit5(), reltol = 1e-5, maxiters = 1e3))
    isequal(independent_variable(M), M.g.a) || error("Independent variable must be $(M.g.a)")

    pars = merge(pars_fixed, Dict(pars_varying .=> NaN))
    as = @. 1 / (zs + 1)
    prob = CosmologyProblem(M, pars; pt = false, ivspan = (minimum(as), 1.0))
    probgen = SymBoltz.parameter_updater(prob, pars_varying; build_initializeprob = Val{false})

    geta = getsym(prob.bg, M.g.a)
    getτ = getsym(prob.bg, M.τ)
    geth = getsym(prob.bg, M.g.h)
    getΩk0 = getsym(prob.bg, M.K.Ω₀)

    return p -> begin
        prob = probgen(p)
        sol = solve(prob; bgopts, saveat = as, save_end = true)
        a = geta(sol.bg)
        τ = getτ(sol.bg)
        h = geth(sol.bg)
        Ωk0 = getΩk0(sol.bg)
        τ0 = τ[end] # time today
        χ = τ0 .- τ # TODO: use M.χ
        r = @. real(sinc(√(-Ωk0+0im)*χ/π) * χ) # Julia's sinc(x) = sin(π*x) / (π*x)
        H0 = SymBoltz.H100 * h
        return @. r / a * SymBoltz.c / H0 # luminosity distance in meters
    end
end
