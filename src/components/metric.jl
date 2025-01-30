"""
    metric(; name = :g, kwargs...)

Create a symbolic component for the perturbed FLRW spacetime metric in the conformal Newtonian gauge with sign signature diag(-1, +1, +1, +1).
"""
function metric(; name = :g, kwargs...)
    vars = a, ℰ, E, H, ℋ, Φ, Ψ, g₁₁, g₂₂, z = GlobalScope.(@variables a(t) ℰ(t) E(t) H(t) ℋ(t) Φ(t) Ψ(t) g₁₁(t) g₂₂(t) z(t))
    pars = h, = GlobalScope.(@parameters h)
    guesses = [
        a => 1e-5
    ]
    description = "Spacetime FLRW metric in Newtonian gauge"
    return ODESystem([
        z ~ 1/a - 1
        ℰ ~ D(a) / a # ℰ = ℋ/ℋ0 = ℋ/H₀
        E ~ ℰ / a # E = H/H₀
        ℋ ~ ℰ * H100 * h
        H ~ E * H100 * h
        g₁₁ ~ a^2 * (-1 - 2*ϵ*Ψ)
        g₂₂ ~ a^2 * (+1 - 2*ϵ*Φ)
    ], t, vars, pars; guesses, name, description, kwargs...)
end
