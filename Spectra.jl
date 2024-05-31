using NumericalIntegration
using Bessels: besselj, sphericalbesselj
using ForwardDiff
using Base.Threads

# matter power spectrum
P0(k, As) = @. 2*π^2 / k ^ 3 * As # TODO: add kpivot and ns
function P(pt::PerturbationsSystem, k, Ωr0, Ωm0, Ωb0, h, As, Yp)
    pt_sols = solve(pt, k, Ωr0, Ωm0, Ωb0, h, Yp)
    ηtoday = pt_sols[1].prob.tspan[end] # TODO: something more robust?
    return P0(k, As) .* pt_sols(ηtoday, idxs=pt.sys.Δm) .^ 2
end

range_until(start, stop, step; skip_start=false) = range(skip_start ? start+step : start, step=step, length=Int(ceil((stop-start)/step+1)))
# TODO: only need as from a = 1e-4 till today
function S(pt::PerturbationsSystem, ηs::AbstractArray, ks::AbstractArray, Ωr0, Ωm0, Ωb0, h, Yp; Sspline_ks=nothing)    
    Ss = zeros(eltype([Ωr0, Ωm0, Ωb0, h, Yp]), (length(ηs), length(ks)))

    if !isnothing(Sspline_ks)
        Ssplinedata = S(pt, ηs, Sspline_ks, Ωr0, Ωm0, Ωb0, h, Yp)
        for i_η in eachindex(ηs)
            Sspline = CubicSpline(Ssplinedata[i_η,:], Sspline_ks)
            Ss[i_η,:] .= Sspline(ks)
        end
        return Ss
    end

    # TODO: use saveat for ηs
    pt_sols = solve(pt, ks, Ωr0, Ωm0, Ωb0, h, Yp)
    for (i_k, pt_sol) in enumerate(pt_sols) # TODO: parallellize over threads?
        # observe Sachs-Wolfe and Doppler contributions from perturbations
        Ss[:,i_k] .+= pt_sol(ηs, idxs=pt.ssys.SSW)
        Ss[:,i_k] .+= pt_sol(ηs, idxs=pt.ssys.SD)
        
        # spline ISW contribution from perturbations (gives wiggles when observed)
        # TODO: use pt_sol(..., Val{1}) when this is fixed: https://github.com/SciML/ModelingToolkit.jl/issues/2697 and https://github.com/SciML/ModelingToolkit.jl/pull/2574
        τs = pt_sol(ηs, idxs=pt.sys.τ).u
        Ψs = pt_sol(ηs, idxs=pt.sys.gravpt.Ψ).u
        Φs = pt_sol(ηs, idxs=pt.sys.gravpt.Φ).u
        Ψ′s = (Ψspl = CubicSpline(Ψs, ηs; extrapolate=true); DataInterpolations.derivative.(Ref(Ψspl), ηs))
        Φ′s = (Φspl = CubicSpline(Φs, ηs; extrapolate=true); DataInterpolations.derivative.(Ref(Φspl), ηs))
        Ss[:,i_k] .+= exp.(.-τs) .* (Ψ′s - Φ′s) # ISW
    end

    return Ss
end

#Ss = S(ηs, ks, Ωr0, Ωm0, Ωb0, h, Yp)
#plot(ηs, asinh.(Ss[:,[1,9]]))

# TODO: contribute back to Bessels.jl
sphericalbesseljslow(ls::AbstractRange, x) = sphericalbesselj.(ls, x)
sphericalbesseljfast(ls::AbstractRange, x) = (x == 0.0 ? 1.0 : √(π/(2*x))) * besselj(ls .+ 0.5, x)

# TODO: integrate CubicSplines instead of trapz! https://discourse.julialang.org/t/how-to-speed-up-the-numerical-integration-with-interpolation/96223/5
function Cl(pt::PerturbationsSystem, ls::AbstractArray, ks::AbstractRange, lnηs::AbstractRange, Ωr0, Ωm0, Ωb0, h, As, Yp; kwargs...)
    ηs = exp.(lnηs)
    Ss = S(pt, ηs, ks, Ωr0, Ωm0, Ωb0, h, Yp; kwargs...) # TODO: reduce memory allocation!
    η0 = ηs[end]

    T = eltype([Ωr0, Ωm0, Ωb0, h, As, Yp]) # TODO: handle with/without As differently?
    ∂Θ_∂lnη = zeros(T, (length(ls), length(ηs)))
    Θls = zeros(T, (length(ls), length(ks)))
    ls_all = minimum(ls):1:maximum(ls) # range with step 1
    ls_indices = 1 .+ ls .- ls_all[1] # indices such that ls = ls_all[ls_indices]
    Jls = zeros((length(ls)))

    # TODO: parallellize some of this over threads
    for (i_k, k) in enumerate(ks)
        for (i_η, η) in enumerate(ηs)
            # TODO: try to spline sphericalbesselj for each l, from x=0 to x=kmax*(η0-ηini)
            Jls[:] .= sphericalbesseljfast(ls_all, k*(η0-η))[ls_indices] # @btime reference (~1000 ls): 5.4 s # TODO: use a grid of constant k*(η0-η)=y; use η with logarithmic spacing; use resulting k = y/(η0-η)
            #Jls .= sphericalbesselj.(ls, k*(η0-η)) # @btime reference (~1000 ls): 32 s
            ∂Θ_∂lnη[:, i_η] .= Ss[i_η, i_k] * Jls[:] * η # TODO: sphericalbesselj takes a lot of time!
        end
        @threads for i_l in eachindex(ls)
            Θls[i_l, i_k] = integrate(lnηs, ∂Θ_∂lnη[i_l, :], SimpsonEven()) # integrate over η # TODO: add starting Θl(ηini) # TODO: calculate ∂Θ_∂logΘ and use Even() methods
        end
    end

    ks_with0 = [0.0; ks]
    dCl_dks_with0 = zeros(T, (length(ks_with0), nthreads()))
    Cls = zeros(T, length(ls))

    @threads for i_l in eachindex(ls)
        dCl_dks_with0[2:end, threadid()] .= @. 2/π * ks^2 * P0(ks, As) * Θls[i_l, :]^2
        Cls[i_l] = integrate(ks_with0, dCl_dks_with0[:, threadid()], SimpsonEven()) # integrate over k (_with0 adds one additional point at (0,0))
    end

    return Cls
end

function Cl(pt::PerturbationsSystem, ls::AbstractArray, Ωr0, Ωm0, Ωb0, h, As, Yp; splineS = true)
    bg_sol = solve(pt.bg, Ωr0, Ωm0)

    ηi, η0 = bg_sol[η][begin], bg_sol[η][end]
    ηi, η0 = ForwardDiff.value.([ηi, η0]) # TODO: do I lose some gradient information here?!
    lnηs = range(log(ηi), log(η0), length=400) # logarithmic spread to capture early-time oscillations # TODO: dynamic/adaptive spacing!

    kmax = 2.0 * ls[end]
    ks = range_until(0, kmax, 2π/4; skip_start=true) ./ η0
    if splineS
        Sspline_ks = range_until(0, kmax, 50; skip_start=true) ./ η0 # Δk = 50/η0
    else
        Sspline_ks = nothing
    end

    return Cl(pt, ls, ks, lnηs, Ωr0, Ωm0, Ωb0, h, As, Yp; Sspline_ks)
end

Dl(pt::PerturbationsSystem, ls::AbstractArray, args...; kwargs...) = Cl(pt, ls, args...; kwargs...) .* ls .* (ls .+ 1) / (2*π)

#Dls = Dl(ls, ks, ηs, Ωr0, Ωm0, Ωb0, h, As, Yp; Sspline_ks)
#plot(ls, Dls; xlabel="l", ylabel="Dl = l (l+1) Cl / 2π")