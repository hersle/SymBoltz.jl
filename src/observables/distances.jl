@doc raw"""
    distance_luminosity(χ, a, h, Ωk0 = 0)

Compute luminosity distances (in meters)
```math
d_L = \frac{c}{H_0} \frac{r}{a}, \quad \mathrm{where} \quad r = \chi \, \frac{\sin\left(\sqrt{-Ω_{k0}} \, \chi\right)}{\sqrt{-Ω_{k0}} \, \chi},
```
from conformal lookback times `χ`, scale factors `a`, Hubble parameter `h` and curvature density `Ωk0`.

!!! warning
    This function will be removed in the next version. Use [`distance`](@ref) with the mode `:L` instead.
"""
function distance_luminosity(χ, a, h, Ωk0 = 0)
    Base.depwarn("distance_luminosity is deprecated and will be removed in the next version. Use distance(:L, sol, t) instead.", :distance_luminosity; force = true)
    H0 = H100 * h
    r = @. real(sinc(√(-Ωk0+0im)*χ/π) * χ) # Julia's sinc(x) = sin(π*x) / (π*x)
    return @. r / a * c / H0 # to meters
end

@doc raw"""
    distance(modes, sol::CosmologySolution, t)

Get the distances `modes` at the values `t` of the independent variable from the solution `sol`.
The times `t` can be a vector of values of the independent variable, or a pair like `M.g.z => zs` with values of another variable.
The `modes` can be one or a vector of:
- `:χ`: line-of-sight comoving distance ``χ = τ₀ - τ``,
- `:M`: transverse comoving distance ``D_M = \sin(\sqrt{-Ω_{k0}} χ) / \sqrt{-Ω_{k0}}``,
- `:A`: angular diameter distance ``D_A = a D_M``,
- `:L`: luminosity distance ``D_L = D_M / a``,
- `:H`: Hubble distance ``D_H = 1 / H``,
- `:V`: volume-averaged distance ``D_V = (z D_M^2 D_H)^{1/3}``.
Returned distances are dimensionless in units of ``c/H₀``.
The model must have the variables `χ`, `a` and `H`.
"""
function distance(modes::AbstractVector, sol::CosmologySolution, t::Union{AbstractVector, Pair})
    bg = sol.bg[end]
    Ωk0 = SymbolicIndexingInterface.is_parameter(bg, :Ωk0) ? bg.ps[:Ωk0] : 0.0 # flat without curvature parameter # TODO: curvature(sol) function that returns K?
    χ, a, H = eachrow(sol([:χ, :a, :H], t))
    DM = @. real(sinc(√(-Ωk0+0im)*χ/π)) * χ # Julia's sinc(x) = sin(π*x) / (π*x)
    return [distance(mode, χ[i], DM[i], a[i], H[i]) for mode in modes, i in eachindex(χ)]
end
function distance(mode::Symbol, χ, DM, a, H)
    mode == :χ && return χ
    mode == :M && return DM
    mode == :A && return DM * a
    mode == :L && return DM / a
    mode == :H && return 1 / H
    mode == :V && return cbrt((1/a-1) * DM^2 / H)
    error("Unknown distance mode :$mode. Use :χ, :M, :A, :L, :H or :V.")
end
distance(mode::Symbol, sol::CosmologySolution, t) = distance([mode], sol, t)[1, :] # fallback with single mode

@doc raw"""
    sound_horizon(sol::CosmologySolution[, t])

Get the photon-baryon sound horizon
```math
    rₛ(τ) = ∫_0^τ dτ cₛ = ∫_0^τ \frac{dτ}{√(3(1+3ρ_b/4ρ_γ))}
```
at the time steps of the solution `sol`, or at the values `t` of the independent variable.
It is integrated with the background as the variable `rₛ`.
"""
sound_horizon(sol::CosmologySolution) = sol[sol.prob.M.rₛ]
sound_horizon(sol::CosmologySolution, t) = sol(sol.prob.M.rₛ, t)

@doc raw"""
    time_drag(sol::CosmologySolution)

Get the value of the independent time variable ``t`` at the baryon drag epoch, when the drag optical depth satisfies
```math
    κ_d(t) = -∫_t^{t_0} \frac{κ'}{3ρ_b/4ρ_γ} dt' = 1,
```
where ``κ`` is the Thomson optical depth and ``κ' = dκ/dt``.
The sound horizon at the drag epoch ``r_d`` is then `sound_horizon(sol, time_drag(sol))`.
"""
time_drag(sol::CosmologySolution) = timeseries(sol, sol.prob.M.κd, 1.0)
