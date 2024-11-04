"""
    metric(; name = :g, kwargs...)

Create a symbolic component for the perturbed FLRW spacetime metric in the conformal Newtonian gauge with signature

    a^2 * diag(-1-2*ϵ*Ψ, +1 - 2*ϵ*Φ, +1 - 2*ϵ*Φ, +1 - 2*ϵ*Φ)
"""
function metric(; name = :g, kwargs...)
    vars = a, ℰ, E, H, ℋ, Φ, Ψ, Φ′, Ψ′, g11, g22, h11, h22, z = GlobalScope.(@variables a(t) ℰ(t) E(t) H(t) ℋ(t) Φ(t) Ψ(t) Φ′(t) Ψ′(t) g11(t) g22(t) h11(t) h22(t) z(t))
    pars = H0, h = GlobalScope.(@parameters H0 h)
    defs = [
        H0 => H100 * h
        h => H0 / H100
    ]
    description = "Spacetime FLRW metric in Newtonian gauge"
    return ODESystem([
        z ~ 1/a - 1
        ℰ ~ D(a) / a # ℰ = ℋ/ℋ0 = ℋ/H0
        E ~ ℰ / a # E = H/H0
        ℋ ~ ℰ * H0
        H ~ E * H0
        g11 ~ -a^2
        g22 ~ +a^2
        h11 * ϵ ~ -2 * a^2 * Ψ * ϵ
        h22 * ϵ ~ -2 * a^2 * Φ * ϵ
        Φ′ * ϵ ~ D(Φ) * ϵ # TODO: why unstable at early times?
        Ψ′ * ϵ ~ D(Ψ) * ϵ # TODO: why unstable at early times?
    ], t, vars, pars; defaults = defs, name, description, kwargs...)
end
