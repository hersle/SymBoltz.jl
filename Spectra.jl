using NumericalIntegration
using Bessels: besselj!, sphericalbesselj
using ForwardDiff
using Base.Threads

# primordial power spectrum
P0(k, As) = @. 2*π^2 / k^3 * As # TODO: add kpivot and ns

# total matter power spectrum
function P(pt::PerturbationsSystem, k, Ωr0, Ωm0, Ωb0, h, As, Yp)
    pt_sols = solve(pt, k, Ωr0, Ωm0, Ωb0, h, Yp)
    ηtoday = pt_sols[1].prob.tspan[end] # TODO: something more robust?
    return P0(k, As) .* pt_sols(ηtoday, idxs=pt.sys.Δm) .^ 2
end

# source function
# this one is more elegant, but a little numerically unstable (would really like to use this one)
function S_observed(pt::PerturbationsSystem, ηs::AbstractArray, ks::AbstractArray, Ωr0, Ωm0, Ωb0, h, Yp)
    pt_sols = solve(pt, ks, Ωr0, Ωm0, Ωb0, h, Yp; saveat = ηs)
    Ss = zeros(eltype([Ωr0, Ωm0, Ωb0, h, Yp]), (length(ηs), length(ks))) # TODO: change order to get DenseArray during integrations?
    @threads for ik in eachindex(ks)
        pt_sol = pt_sols[ik]
        Ss[:,ik] .= pt_sol[pt.ssys.S] # whether this gives a accurate CMB spectrum depends on the perturbation ODE solver (e.g. KenCarp{4,47,5,58}) and its reltol
    end
    return Ss
end

# this one is less elegant, but more numerically stable
function S_splined(pt::PerturbationsSystem, ηs::AbstractArray, ks::AbstractArray, Ωr0, Ωm0, Ωb0, h, Yp)
    th = pt.th
    th_sol = solve(th, Ωr0, Ωm0, Ωb0, h, Yp; saveat = ηs)
    τ = th_sol[th.sys.τ] .- th_sol[th.sys.τ][end] # make τ = 0 today # TODO: assume ηs[end] is today
    τ′ = D_spline(τ, ηs)
    τ″ = D_spline(τ′, ηs)
    g = @. -τ′ * exp(-τ)
    g′ = @. (τ′^2 - τ″) * exp(-τ)
    
    # TODO: add source functions as observed perturbation functions? but difficult with cumulative τ(η)? must anyway wait for this to be fixed: https://github.com/SciML/ModelingToolkit.jl/issues/2697
    pt_sols = solve(pt, ks, Ωr0, Ωm0, Ωb0, h, Yp; saveat = ηs)
    Ss = zeros(eltype([Ωr0, Ωm0, Ωb0, h, Yp]), (length(ηs), length(ks))) # TODO: change order to get DenseArray during integrations?
    @threads for ik in eachindex(ks)
        pt_sol = pt_sols[ik]
        k = ks[ik]
        Θ0 = pt_sol[pt.sys.ph.Θ0]
        Ψ = pt_sol[pt.sys.gravpt.Ψ]
        Φ = pt_sol[pt.sys.gravpt.Φ]
        Π = pt_sol[pt.sys.ph.Π]
        ub = pt_sol[pt.sys.bar.u]
        Ψ′ = D_spline(Ψ, ηs) # TODO: use pt_sol(..., Val{1}) when this is fixed: https://github.com/SciML/ModelingToolkit.jl/issues/2697 and https://github.com/SciML/ModelingToolkit.jl/pull/2574
        Φ′ = D_spline(Φ, ηs)
        ub′ = D_spline(ub, ηs)
        gΠ″ = D_spline(g .* Π, ηs; order = 2)
        @. Ss[:,ik] = g*(Θ0+Ψ+Π/4) + (g′*ub+g*ub′)/k + exp(-τ)*(Ψ′-Φ′) + 3/(4*k^2)*gΠ″ # SW + Doppler + ISW + polarization
    end

    return Ss
end

function S(pt::PerturbationsSystem, ηs::AbstractArray, ksfine::AbstractArray, Ωr0, Ωm0, Ωb0, h, Yp, kscoarse::AbstractArray; Spline = CubicSpline, kwargs...) 
    Sscoarse = S_splined(pt, ηs, kscoarse, Ωr0, Ωm0, Ωb0, h, Yp)
    Ssfine = similar(Sscoarse, (length(ηs), length(ksfine)))
    for iη in eachindex(ηs)
        Sscoarseη = @view Sscoarse[iη,:]
        Ssfine[iη,:] .= Spline(Sscoarseη, kscoarse)(ksfine)
    end
    return Ssfine
end

#Ss = S(ηs, ks, Ωr0, Ωm0, Ωb0, h, Yp)
#plot(ηs, asinh.(Ss[:,[1,9]]))

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
# TODO: spline sphericalbesselj for each l, from x=0 to x=kmax*(η0-ηini)
# TODO: integrate with ApproxFun? see e.g. https://discourse.julialang.org/t/evaluate-integral-on-many-points-cubature-jl/1723/2
# TODO: RombergEven() works with 513 or 1025 points (do Logging.disable_logging(Logging.Warn) first)
# line of sight integration
function Θl(ls::AbstractArray, ks::AbstractRange, lnηs::AbstractRange, Ss::AbstractArray; integrator = SimpsonEven())
    ηs = exp.(lnηs)
    η0 = ηs[end]

    Θls = similar(Ss, (length(ks), length(ls)))
    ∂Θ_∂lnη = [similar(Ss, (length(ls), length(ηs))) for _ in 1:nthreads()] # TODO: best to do array of arrays without using @view, or to use matrix + @view?
    lmin, lmax = extrema(ls)
    ls_all = lmin:1:lmax # range with step 1
    Jls_all = [zeros(length(ls_all)) for _ in 1:nthreads()] # separate workspace per thread

    @threads for ik in eachindex(ks)
        k = ks[ik]
        for iη in eachindex(ηs)
            η = ηs[iη]
            Sη = Ss[iη,ik] * η
            kΔη = k * (η0-η)
            sphericalbesseljfast!(Jls_all[threadid()], ls_all, kΔη) # TODO: reuse ∂Θ_∂lnη's memory?
            for il in eachindex(ls) # TODO: @simd if I can make Jl access with unit stride? also need @inbounds?
                l = ls[il]
                Jl = Jls_all[threadid()][1+l-lmin]
                ∂Θ_∂lnη[threadid()][il,iη] = Sη * Jl
                # TODO: integrate in this loop instead?
            end
        end
        for il in eachindex(ls)
            integrand = @view ∂Θ_∂lnη[threadid()][il,:]
            Θls[ik,il] = integrate(lnηs, integrand, integrator) # integrate over η # TODO: add starting Θl(ηini) # TODO: calculate ∂Θ_∂logΘ and use Even() methods
        end
    end

    return Θls
end

function Θl(pt::PerturbationsSystem, ls::AbstractArray, ks::AbstractRange, lnηs::AbstractRange, Ωr0, Ωm0, Ωb0, h, Yp, args...; kwargs...)
    ηs = exp.(lnηs)
    Ss = S(pt, ηs, ks, Ωr0, Ωm0, Ωb0, h, Yp, args...; kwargs...)
    return Θl(ls, ks, lnηs, Ss)
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

function Cl(pt::PerturbationsSystem, ls::AbstractArray, ks::AbstractRange, lnηs::AbstractRange, Ωr0, Ωm0, Ωb0, h, As, Yp, ks_S::AbstractArray; kwargs...)
    Θls = Θl(pt, ls, ks, lnηs, Ωr0, Ωm0, Ωb0, h, Yp, ks_S; kwargs...)
    P0s = P0(ks, As)
    return Cl(ls, ks, Θls, P0s)
end

function Cl(pt::PerturbationsSystem, ls::AbstractArray, Ωr0, Ωm0, Ωb0, h, As, Yp; Δlnη = 0.03, Δkη0 = 2π/4, Δkη0_S = 50.0, observe = false)
    bg_sol = solve(pt.bg, Ωr0, Ωm0)

    ηi, η0 = max(bg_sol[η][begin], 1e-4), bg_sol[η][end]
    ηi, η0 = ForwardDiff.value.([ηi, η0]) # TODO: do I lose some gradient information here?! no? ηi/η0 is just a shift of the integration interval?
    lnηs = range(log(ηi), log(η0), step=Δlnη) # logarithmic spread to capture early-time oscillations # TODO: dynamic/adaptive spacing!

    kη0max = 2.0 * ls[end]
    ks_Cl = range_until(0, kη0max, Δkη0; skip_start=true) ./ η0
    ks_S = range_until(0, kη0max, Δkη0_S; skip_start=true) ./ η0 # Δk = 50/η0

    return Cl(pt, ls, ks_Cl, lnηs, Ωr0, Ωm0, Ωb0, h, As, Yp, ks_S; observe)
end

Dl(pt::PerturbationsSystem, ls::AbstractArray, args...; kwargs...) = Cl(pt, ls, args...; kwargs...) .* ls .* (ls .+ 1) ./ 2π

#Dls = Dl(ls, ks, ηs, Ωr0, Ωm0, Ωb0, h, As, Yp; Sspline_ks)
#plot(ls, Dls; xlabel="l", ylabel="Dl = l (l+1) Cl / 2π")

function D_spline(y, x; Spline = CubicSpline, order = 1)
    y_spline = Spline(y, x; extrapolate=true)
    y′ = DataInterpolations.derivative.(Ref(y_spline), x, order)
    return y′
end

range_until(start, stop, step; skip_start=false) = range(skip_start ? start+step : start, step=step, length=Int(ceil((stop-start)/step+1)))