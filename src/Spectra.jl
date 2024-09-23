using NumericalIntegration
using Bessels: besselj!, sphericalbesselj
using DataInterpolations
using ForwardDiff
using Base.Threads

# primordial power spectrum
P0(k, As=2e-9) = @. 2*π^2 / k^3 * As # TODO: add kpivot and ns

# power spectrum
function power_spectrum(sol::CosmologySolution, species::ODESystem, k)
    tend = sol[t][end]
    return P0(k) .* sol(k, tend, [species.Δ^2])[:, 1, 1]
end

function power_spectrum(M::CosmologyModel, species::ODESystem, pars, k; kwargs...)
    sol = solve(M, pars, k; save_everystep=false, kwargs...) # just save endpoints
    return power_spectrum(sol, species, k)
end


#= # TODO: make work again?
# source function
# this one is more elegant, but a little numerically unstable (would really like to use this one)
function S_observed(pt::PerturbationsSystem, ts::AbstractArray, ks::AbstractArray, Ωγ0, Ων0, Ωc0, Ωb0, h, Yp)
    pt_sols = solve(pt, ks, Ωγ0, Ων0, Ωc0, Ωb0, h, Yp; saveat = ts)
    Ss = zeros(eltype([Ωγ0, Ων0, Ωc0, Ωb0, h, Yp]), (length(ts), length(ks))) # TODO: change order to get DenseArray during integrations?
    @threads for ik in eachindex(ks)
        pt_sol = pt_sols[ik]
        Ss[:,ik] .= pt_sol[pt.ssys.S] # whether this gives a accurate CMB spectrum depends on the perturbation ODE solver (e.g. KenCarp{4,47,5,58}) and its reltol
    end
    return Ss
end
=#

# this one is less elegant, but more numerically stable
function S_splined(M::CosmologyModel, ts::AbstractArray, ks::AbstractArray, pars; kwargs...)
    sol = solve(M, pars, ks; saveat = ts, kwargs...)
    τ = sol[M.b.rec.τ] .- sol[M.b.rec.τ][end] # make τ = 0 today # TODO: assume ts[end] is today
    τ′ = D_spline(τ, ts)
    g = @. -τ′ * exp(-τ)
    
    # TODO: add source functions as observed perturbation functions?
    Ss = zeros(eltype([par[2] for par in pars]), (length(ts), length(ks))) # TODO: change order to get DenseArray during integrations?
    @threads for ik in eachindex(ks)
        k = ks[ik]
        Θ0 = sol[ik, M.γ.F0]
        Ψ = sol[ik, M.g.Ψ]
        Φ = sol[ik, M.g.Φ]
        Π = sol[ik, M.γ.Π]
        ub = sol[ik, M.b.u]
        gub′ = D_spline(g .* ub, ts) # TODO: why is the observed version messy?
        Ψ_minus_Φ′ = D_spline(Ψ .- Φ, ts) # TODO: use pt_sol(..., Val{1}) when this is fixed: https://github.com/SciML/ModelingToolkit.jl/issues/2697 and https://github.com/SciML/ModelingToolkit.jl/pull/2574 # TODO: why is the observed version so messy?
        gΠ″ = D_spline(g .* Π, ts; order = 2) # TODO: why is the observed version messy?
        @. Ss[:,ik] = g*(Θ0+Ψ+Π/4) + gub′/k + exp(-τ)*Ψ_minus_Φ′ + 3/(4*k^2)*gΠ″ # SW + Doppler + ISW + polarization
    end

    return Ss
end

function S(M::CosmologyModel, ts::AbstractArray, ksfine::AbstractArray, pars, kscoarse::AbstractArray; kwargs...)
    Sscoarse = S_splined(M::CosmologyModel, ts, kscoarse, pars) # TODO: restore S_observed
    Ssfine = similar(Sscoarse, (length(ts), length(ksfine)))
    for it in eachindex(ts)
        Sscoarset = @view Sscoarse[it,:]
        Ssfine[it,:] .= spline(Sscoarset, kscoarse)(ksfine)
    end
    return Ssfine
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

# TODO: line-of-sight integrate Θl using ODE for evolution of Jl?
# TODO: spline sphericalbesselj for each l, from x=0 to x=kmax*(t0-tini)
# TODO: integrate with ApproxFun? see e.g. https://discourse.julialang.org/t/evaluate-integral-on-many-points-cubature-jl/1723/2
# TODO: RombergEven() works with 513 or 1025 points (do Logging.disable_logging(Logging.Warn) first)
# TODO: gaussian quadrature with weight function? https://juliamath.github.io/QuadGK.jl/stable/weighted-gauss/
# line of sight integration
function Θl(ls::AbstractArray, ks::AbstractRange, lnts::AbstractRange, Ss::AbstractArray; integrator = SimpsonEven())
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
            St = Ss[it,ik] * t
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

function Θl(M::CosmologyModel, ls::AbstractArray, ks::AbstractRange, lnts::AbstractRange, pars, args...; kwargs...)
    ts = exp.(lnts)
    Ss = S(M, ts, ks, pars, args...; kwargs...)
    return Θl(ls, ks, lnts, Ss)
end

# TODO: integrate splines instead of trapz! https://discourse.julialang.org/t/how-to-speed-up-the-numerical-integration-with-interpolation/96223/5
function Cl(ls::AbstractArray, ks::AbstractRange, Θls::AbstractArray, P0s::AbstractArray; integrator = SimpsonEven())
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

function Cl(M::CosmologyModel, pars, ls::AbstractArray, ks::AbstractRange, lnts::AbstractRange, ks_S::AbstractArray; kwargs...)
    Θls = Θl(M, ls, ks, lnts, pars, ks_S; kwargs...)
    P0s = P0(ks)
    return Cl(ls, ks, Θls, P0s)
end

function Cl(M::CosmologyModel, pars, ls::AbstractArray; Δlnt = 0.03, Δkt0 = 2π/4, Δkt0_S = 50.0, observe = false)
    bg_sol = solve(M, pars)

    ti, t0 = 1e-4, bg_sol[t][end] # add tiny number to ti; otherwise the lengths of ts and ODESolution(... ; saveat = ts) differs by 1
    ti, t0 = ForwardDiff.value.([ti, t0]) # TODO: do I lose some gradient information here?! no? ti/t0 is just a shift of the integration interval?
    lnts = range(log(ti), log(t0), step=Δlnt) # logarithmic spread to capture early-time oscillations # TODO: dynamic/adaptive spacing!

    kt0max = 2.0 * ls[end]
    ks_Cl = range_until(0, kt0max, Δkt0; skip_start=true) ./ t0
    ks_S = range_until(0, kt0max, Δkt0_S; skip_start=true) ./ t0 # Δk = 50/t0

    return Cl(M, pars, ls, ks_Cl, lnts, ks_S; observe)
end
