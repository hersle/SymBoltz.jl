"""
    metric(; name = :g, kwargs...)

Create a symbolic component for the perturbed FLRW spacetime metric in the conformal Newtonian gauge with sign signature diag(-1, +1, +1, +1).
"""
function metric(; name = :g, kwargs...)
    vars = @variables begin
        a(τ), [description = "Scale factor"]
        ℋ(τ), [description = "Conformal Hubble function (in units of H₀)"]
        H(τ), [description = "Cosmic Hubble function (in units of H₀)"]
        Ψ(τ, k), [description = "Gravitational metric potential in gₜₜ = -a²(1+2Ψ)"]
        Φ(τ, k), [description = "Curvature metric potential in gᵢⱼ = a²(1-2Φ)δᵢⱼ"]
        Ψ̇(τ, k), [description = "Derivative of Φ wrt. conformal time"]
        Φ̇(τ, k), [description = "Derivative of Ψ wrt. conformal time"]
        z(τ), [description = "Redshift"]
        ż(τ), [description = "Redshift derivative"]
    end
    vars = a, ℋ, H, Ψ, Φ, Ψ̇, Φ̇, z, ż = GlobalScope.(vars)

    pars = h, = GlobalScope.(@parameters h [description = "Dimensionless Hubble parameter today (H₀/(100km/s/Mpc))"])
    description = "Spacetime FLRW metric in Newtonian gauge"
    return System([
        z ~ 1/a - 1
        ż ~ D(z)
        ℋ ~ D(a) / a
        H ~ ℋ / a
    ], τ, vars, pars; name, description, kwargs...)
end
