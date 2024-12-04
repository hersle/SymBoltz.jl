using NumericalIntegration
using Bessels: besselj!, sphericalbesselj
using DataInterpolations
using ForwardDiff
using Base.Threads
using TwoFAST

"""
    P0(k; As=2e-9)

Compute the primordial power spectrum with amplitude `As` at the wavenumber(s) `k`.
"""
P0(k; As=2e-9) = @. 2*π^2 / k^3 * As # TODO: add kpivot and ns # TODO: make separate InflationModel with these parameters

"""
    power_spectrum(sol::CosmologySolution, k)

Compute the matter power spectrum from the cosmology solution `sol` at wavenumber(s) `k`.
"""
function power_spectrum(sol::CosmologySolution, k)
    tend = sol[t][end]
    M = sol.pts[1].prob.f.sys
    ρm = M.c.ρ + M.b.ρ # TODO: massive neutrinos
    return P0(k) .* sol(k, tend, [M.k^2*M.g.Φ / (4*Num(π)*M.g.a^2*ρm)])[:, 1, 1] .^ 2 # Baumann (4.4.172)
end

"""
    power_spectrum(M::CosmologyModel, pars, k; solver = KenCarp4(), kwargs...)

Compute the matter power spectrum from the cosmological model `M` with parameter `pars` at wavenumber(s) `k`.
The `solver` and other `kwargs` are passed to `solve`.
"""
function power_spectrum(M::CosmologyModel, pars, k; solver = KenCarp4(), kwargs...)
    sol = solve(M, pars, k; save_everystep=false, solver, kwargs...) # just save endpoints
    return power_spectrum(sol, k)
end

# this one is less elegant, but more numerically stable?
# TODO: saveat = ts
function S_splined(sol::CosmologySolution, ks::AbstractArray, ts::AbstractArray)
    M = sol.M
    Ss = zeros((length(ks), length(ts))) # TODO: change order to get DenseArray during integrations?

    τ = sol(ts, M.b.rec.τ) # TODO: assume ts[end] is today
    v = -D_spline(τ, ts) .* exp.(-τ)
    @threads for ik in eachindex(ks)
        k = ks[ik]
        idxs = [M.γ.δ, M.g.Ψ, M.γ.Π, M.g.Φ, M.b.u]
        out = sol(k, ts, idxs)
        δ, Ψ, Π, Φ, ub = selectdim.(Ref(out), 2, eachindex(idxs))
        Ss[ik,:] .= v .* (δ/4 + Ψ + Π/4) + exp.(-τ) .* D_spline(Ψ + Φ, ts) + D_spline(v .* ub, ts) / k + 3/(4*k^2) * D_spline(v .* Π, ts; order = 2) # Dodelson (9.57) with Φ → -Φ and polarization
        #Ss[ik,:] .= exp.(-τ) .* (D_spline(Φ, ts) - τ̇/4 .* (δ + Π)) + D_spline(exp.(-τ) .* (Ψ - ub.*τ̇/k), ts) + 3/(4*k^2) * D_spline(v .* Π, ts; order = 2) # Dodelson (9.55) with Φ → -Φ
        #Ss[ik,:] .= v .* (δ/4 + Ψ + Π/4) + v .* (Φ-Ψ) + 2 * exp.(-τ) .* D_spline(Φ, ts) + D_spline(v .* ub, ts) / k + exp.(-τ) * k .* (Ψ - Φ) + 3/(4*k^2) * D_spline(v .* Π, ts; order = 2) # CLASS' expression with added polarization
    end

    return Ss
end

# TODO: contribute back to Bessels.jl
#sphericalbesseljslow(ls::AbstractArray, x) = sphericalbesselj.(ls, x)
#sphericalbesseljfast(ls::AbstractRange, x) = (x == 0.0 ? 1.0 : √(π/(2*x))) * besselj(ls .+ 0.5, x)
function sphericalbesseljfast!(out, ls::AbstractRange, x)
    besselj!(out, ls .+ 0.5, x)
    if x != 0.0 # TODO: well-behaved?
        @. out *= √(π/(2*x))
    end
    return out
end
# TODO: test and compare to fast!
function sphericalbesseljslow!(out, ls::AbstractRange, x)
    out .= sphericalbesselj.(ls, x)
end

# TODO: line-of-sight integrate Θl using ODE for evolution of Jl?
# TODO: spline sphericalbesselj for each l, from x=0 to x=kmax*(t0-tini)
# TODO: integrate with ApproxFun? see e.g. https://discourse.julialang.org/t/evaluate-integral-on-many-points-cubature-jl/1723/2
# TODO: RombergEven() works with 513 or 1025 points (do Logging.disable_logging(Logging.Warn) first)
# TODO: gaussian quadrature with weight function? https://juliamath.github.io/QuadGK.jl/stable/weighted-gauss/
# line of sight integration
# TODO: take in symbolic expr?
# TODO: use tasks, not threads!
function Θl(Ss::AbstractArray, ls::AbstractArray, ks::AbstractRange, lnts::AbstractRange; integrator = SimpsonEven(), verbose = true)
    @assert size(Ss) == (length(ks), length(lnts)) # TODO: optimal structure?
    verbose && println("LOS integration with $(length(ls)) ls x $(length(ks)) ks x $(length(lnts)) lnts")

    ts = exp.(lnts)
    t0 = ts[end]

    Θls = similar(Ss, (length(ks), length(ls)))
    ∂Θ_∂lnt = [similar(Ss, (length(ls), length(ts))) for _ in 1:nthreads()] # TODO: best to do array of arrays without using @view, or to use matrix + @view?
    lmin, lmax = extrema(ls)
    ls_all = lmin:1:lmax # range with step 1
    Jls_all = [zeros(length(ls_all)) for _ in 1:nthreads()] # separate workspace per thread

    @threads for ik in eachindex(ks)
        k = ks[ik]
        for it in eachindex(ts)
            t = ts[it]
            St = Ss[ik,it] * t # multiply S(k,t) to get (S(k, log(t)))
            kΔt = k * (t0-t)
            sphericalbesseljfast!(Jls_all[threadid()], ls_all, kΔt) # TODO: reuse ∂Θ_∂lnt's memory?
            for il in eachindex(ls) # TODO: @simd if I can make Jl access with unit stride? also need @inbounds?
                l = ls[il]
                Jl = Jls_all[threadid()][1+l-lmin]
                ∂Θ_∂lnt[threadid()][il,it] = St * Jl
                # TODO: integrate in this loop instead?
            end
        end
        for il in eachindex(ls)
            integrand = @view ∂Θ_∂lnt[threadid()][il,:]
            Θls[ik,il] = integrate(lnts, integrand, integrator) # integrate over t # TODO: add starting Θl(tini) # TODO: calculate ∂Θ_∂logΘ and use Even() methods
        end
    end

    return Θls
end

function Θl(sol::CosmologySolution, ls::AbstractArray, ks::AbstractRange, lnts::AbstractRange, pars, ks_S; kwargs...)
    Ss = S_splined(sol, ks, exp.(lnts)) # sol(ks, exp.(lnts), M.S)
    return Θl(Ss, ls, ks, lnts; kwargs...)
end

# TODO: integrate splines instead of trapz! https://discourse.julialang.org/t/how-to-speed-up-the-numerical-integration-with-interpolation/96223/5
function Cl(Θls::AbstractArray, P0s::AbstractArray, ls::AbstractArray, ks::AbstractRange; integrator = SimpsonEven())
    Cls = similar(Θls, length(ls))
    ks_with0 = [0.0; ks] # add dummy value with k=0 for integration
    dCl_dks_with0 = [zeros(eltype(Θls), length(ks_with0)) for _ in 1:nthreads()] # separate workspace per thread

    @threads for il in eachindex(ls)
        integrand = @view Θls[:,il]
        @. dCl_dks_with0[threadid()][2:end] = 2/π * ks^2 * P0s * integrand^2
        Cls[il] = integrate(ks_with0, dCl_dks_with0[threadid()], integrator) # integrate over k (_with0 adds one additional point at (0,0))
    end

    return Cls
end

# TODO: try to disable CLASS' want_lcmb_full_limber?
function Cl(sol::CosmologySolution, ls::AbstractArray, ks::AbstractArray, lnts::AbstractArray; kwargs...)
    Ss = S_splined(sol, ks, exp.(lnts)) # TODO: restore sol(ks, exp.(lnts), sol.M.S)
    Θls = Θl(Ss, ls, ks, lnts; kwargs...)
    P0s = P0(ks)
    return Cl(Θls, P0s, ls, ks)
end


"""
    Cl(M::CosmologyModel, pars::Dict, ls::AbstractArray; kwargs...)

Compute the ``C_l``'s of the CMB power spectrum from the cosmological model `M` with parameters `pars` at angular wavenumbers `ls`.
"""
function Cl(M::CosmologyModel, pars::Dict, ls::AbstractArray; integrator = SimpsonEven(), kwargs...) # TODO: Δlnt shifts Cls <->, Δkt0 seems fine, should test interpolation with Δkt0_S! kt0max_lmax?
    @assert issorted(ls)
    sol, ks, lnts = solve_for_Cl(M, pars, ls[end]; kwargs...)
    return Cl(sol, ls, ks, lnts; integrator)
end

function solve_for_Cl(M::CosmologyModel, pars::Dict, lmax; Δk = 2π/24, Δk_S = 10.0, kmax = 1.0 * lmax, Δlnt=0.03, kwargs...)
    # Assumes t0 = 1 (e.g. t0 = 1/H0 = 1) # TODO: don't assume t0 = 1
    kmin = Δk
    ks = range(kmin, kmax, length = Int(floor((kmax-kmin)/Δk_S+1)))
    sol = solve(M, pars, ks; kwargs...)

    ti, t0 = sol[t][begin] + 1e-10, sol[t][end] # add tiny number to ti; otherwise the lengths of ts and ODESolution(... ; saveat = ts) differs by 1
    ti, t0 = ForwardDiff.value.([ti, t0]) # TODO: do I lose some gradient information here?! no? ti/t0 is just a shift of the integration interval?
    lnts = range(log(ti), log(t0), step=Δlnt) # logarithmic spread to capture early-time oscillations # TODO: dynamic/adaptive spacing!

    ks = range(ks[begin], ks[end], length = length(ks) * Int(floor(Δk_S/Δk))) # use integer multiple so endpoints are the same

    return sol, ks, lnts
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
        P = CubicSpline(power_spectrum(sol, ks), ks) # create spline interpolation (fast)
    else
        P(k) = only(power_spectrum(sol, k)) # use solution's built-in interpolation (elegant)
    end
    kmin, kmax = extrema(ks)
    rmin = 2π / kmax
    return xicalc(P, 0, 0; N, kmin, kmax, r0=rmin)
end
