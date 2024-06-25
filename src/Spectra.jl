using NumericalIntegration
using Bessels: besselj!, sphericalbesselj
using ForwardDiff
using Base.Threads

# primordial power spectrum
P0(k, As) = @. 2*π^2 / k^3 * As # TODO: add kpivot and ns

# total matter power spectrum
function P(pt::ODESystem, ks, Ωγ0, Ων0, Ωc0, Ωb0, h, As, Yp; solver = Rodas5P(), reltol=1e-7, verbose=true)
    pts = structural_simplify(pt)
    th = pts.th
    bg = th.bg

    probs = EnsembleProblem(; safetycopy = false, prob = nothing, prob_func = (_, i, _) -> begin
        verbose && println("$i/$(length(ks)) k = $(ks[i]*k0) Mpc/h")
        return ODEProblem(pts, [], (1e-5, 4.0), [bg.ph.Ω0 => Ωγ0, bg.neu.Ω0 => Ων0, bg.cdm.Ω0 => Ωc0, bg.bar.Ω0 => Ωb0, bg.g.h => h, th.Yp => Yp, k => ks[i]])
    end)
    sols = solve(probs, solver, EnsembleThreads(), trajectories = length(ks); reltol) # TODO: test GPU parallellization
    return P0(ks, As) .* sols(4.0, idxs=pt.Δm) .^ 2
end

#=
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

# this one is less elegant, but more numerically stable
function S_splined(pt::PerturbationsSystem, ts::AbstractArray, ks::AbstractArray, Ωγ0, Ων0, Ωc0, Ωb0, h, Yp)
    th = pt.th
    th_sol = solve(th, Ωγ0, Ων0, Ωc0, Ωb0, h, Yp; saveat = ts)
    τ = th_sol[th.sys.τ] .- th_sol[th.sys.τ][end] # make τ = 0 today # TODO: assume ts[end] is today
    τ′ = D_spline(τ, ts)
    τ″ = D_spline(τ′, ts)
    g = @. -τ′ * exp(-τ)
    g′ = @. (τ′^2 - τ″) * exp(-τ)
    
    # TODO: add source functions as observed perturbation functions? but difficult with cumulative τ(t)? must anyway wait for this to be fixed: https://github.com/SciML/ModelingToolkit.jl/issues/2697
    pt_sols = solve(pt, ks, Ωγ0, Ων0, Ωc0, Ωb0, h, Yp; saveat = ts)
    Ss = zeros(eltype([Ωγ0, Ων0, Ωc0, Ωb0, h, Yp]), (length(ts), length(ks))) # TODO: change order to get DenseArray during integrations?
    @threads for ik in eachindex(ks)
        pt_sol = pt_sols[ik]
        k = ks[ik]
        Θ0 = pt_sol[pt.sys.ph.Θ0]
        Ψ = pt_sol[pt.sys.gravpt.Ψ]
        Φ = pt_sol[pt.sys.gravpt.Φ]
        Π = pt_sol[pt.sys.ph.Π]
        ub = pt_sol[pt.sys.bar.u]
        Ψ′ = D_spline(Ψ, ts) # TODO: use pt_sol(..., Val{1}) when this is fixed: https://github.com/SciML/ModelingToolkit.jl/issues/2697 and https://github.com/SciML/ModelingToolkit.jl/pull/2574
        Φ′ = D_spline(Φ, ts)
        ub′ = D_spline(ub, ts)
        gΠ″ = D_spline(g .* Π, ts; order = 2)
        @. Ss[:,ik] = g*(Θ0+Ψ+Π/4) + (g′*ub+g*ub′)/k + exp(-τ)*(Ψ′-Φ′) + 3/(4*k^2)*gΠ″ # SW + Doppler + ISW + polarization
    end

    return Ss
end

function S(pt::PerturbationsSystem, ts::AbstractArray, ksfine::AbstractArray, Ωγ0, Ων0, Ωc0, Ωb0, h, Yp, kscoarse::AbstractArray; Spline = CubicSpline, kwargs...) 
    Sscoarse = S_observed(pt, ts, kscoarse, Ωγ0, Ων0, Ωc0, Ωb0, h, Yp) # TODO: restore observed
    Ssfine = similar(Sscoarse, (length(ts), length(ksfine)))
    for it in eachindex(ts)
        Sscoarset = @view Sscoarse[it,:]
        Ssfine[it,:] .= Spline(Sscoarset, kscoarse)(ksfine)
    end
    return Ssfine
end

#Ss = S(ts, ks, Ωγ0, Ων0, Ωc0, Ωb0, h, Yp)
#plot(ts, asinh.(Ss[:,[1,9]]))

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

function Θl(pt::PerturbationsSystem, ls::AbstractArray, ks::AbstractRange, lnts::AbstractRange, Ωγ0, Ων0, Ωc0, Ωb0, h, Yp, args...; kwargs...)
    ts = exp.(lnts)
    Ss = S(pt, ts, ks, Ωγ0, Ων0, Ωc0, Ωb0, h, Yp, args...; kwargs...)
    return Θl(ls, ks, lnts, Ss)
end

# TODO: integrate CubicSplines instead of trapz! https://discourse.julialang.org/t/how-to-speed-up-the-numerical-integration-with-interpolation/96223/5
function Cl(ls::AbstractArray, ks::AbstractRange, Θls::AbstractArray, P0s::AbstractArray; integrator = SimpsonEven())
    Cls = similar(Θls, length(ls))
    ks_with0 = [0.0; ks]
    dCl_dks_with0 = [zeros(eltype(Θls), length(ks_with0)) for _ in 1:nthreads()] # separate workspace per thread

    @threads for il in eachindex(ls)
        integrand = @view Θls[:,il]
        @. dCl_dks_with0[threadid()][2:end] = 2/π * ks^2 * P0s * integrand^2
        Cls[il] = integrate(ks_with0, dCl_dks_with0[threadid()], integrator) # integrate over k (_with0 adds one additional point at (0,0))
    end

    return Cls
end

function Cl(pt::PerturbationsSystem, ls::AbstractArray, ks::AbstractRange, lnts::AbstractRange, Ωγ0, Ων0, Ωc0, Ωb0, h, As, Yp, ks_S::AbstractArray; kwargs...)
    Θls = Θl(pt, ls, ks, lnts, Ωγ0, Ων0, Ωc0, Ωb0, h, Yp, ks_S; kwargs...)
    P0s = P0(ks, As)
    return Cl(ls, ks, Θls, P0s)
end

function Cl(pt::PerturbationsSystem, ls::AbstractArray, Ωγ0, Ων0, Ωc0, Ωb0, h, As, Yp; Δlnt = 0.03, Δkt0 = 2π/4, Δkt0_S = 50.0, observe = false)
    bg_sol = solve(pt.bg, Ωγ0, Ων0, Ωc0, Ωb0)

    ti, t0 = max(bg_sol[t][begin], 1e-4), bg_sol[t][end]
    ti, t0 = ForwardDiff.value.([ti, t0]) # TODO: do I lose some gradient information here?! no? ti/t0 is just a shift of the integration interval?
    lnts = range(log(ti), log(t0), step=Δlnt) # logarithmic spread to capture early-time oscillations # TODO: dynamic/adaptive spacing!

    kt0max = 2.0 * ls[end]
    ks_Cl = range_until(0, kt0max, Δkt0; skip_start=true) ./ t0
    ks_S = range_until(0, kt0max, Δkt0_S; skip_start=true) ./ t0 # Δk = 50/t0

    return Cl(pt, ls, ks_Cl, lnts, Ωγ0, Ων0, Ωc0, Ωb0, h, As, Yp, ks_S; observe)
end

Dl(pt::PerturbationsSystem, ls::AbstractArray, args...; kwargs...) = Cl(pt, ls, args...; kwargs...) .* ls .* (ls .+ 1) ./ 2π

#Dls = Dl(ls, ks, ts, Ωγ0, Ων0, Ωc0, Ωb0, h, As, Yp; Sspline_ks)
#plot(ls, Dls; xlabel="l", ylabel="Dl = l (l+1) Cl / 2π")

function D_spline(y, x; Spline = CubicSpline, order = 1)
    y_spline = Spline(y, x; extrapolate=true)
    y′ = DataInterpolations.derivative.(Ref(y_spline), x, order)
    return y′
end

range_until(start, stop, step; skip_start=false) = range(skip_start ? start+step : start, step=step, length=Int(ceil((stop-start)/step+1)))
=#