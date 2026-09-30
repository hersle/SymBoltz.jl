@doc raw"""
    distance_luminosity(χ, a, h, Ωk0 = 0)

Compute luminosity distances (in meters)
```math
d_L = \frac{c}{H_0} \frac{r}{a}, \quad \mathrm{where} \quad r = \chi \, \frac{\sin\left(\sqrt{-Ω_{k0}} \, \chi\right)}{\sqrt{-Ω_{k0}} \, \chi},
```
from conformal lookback times `χ`, scale factors `a`, Hubble parameter `h` and curvature density `Ωk0`.
"""
function distance_luminosity(χ, a, h, Ωk0 = 0)
    H0 = H100 * h
    r = @. real(sinc(√(-Ωk0+0im)*χ/π) * χ) # Julia's sinc(x) = sin(π*x) / (π*x)
    return @. r / a * c / H0 # to meters
end

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
