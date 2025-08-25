# TODO: incorporate in metric/gravity; don't define as effective species
"""
    curvature(g; name = :K, kwargs...)

Create a species that effectively accounts for curvature in the spacetime with metric `g`.
"""
function curvature(g; name = :K, kwargs...)
    description = "Curvature"
    vars = @variables begin
        ρ(τ), [description = "Effective background density"]
        P(τ), [description = "Effective background pressure"]
        w(τ), [description = "Effective equation of state"]
        δ(τ, k),
        θ(τ, k),
        cₛ²(τ),
        σ(τ, k)
    end
    pars = @parameters begin
        Ω₀, [description = "Effective reduced background density today"] # dimless K is physical K*c²/H₀²
    end
    eqs = [
        w ~ -1//3
        ρ ~ 3/(8*Num(π))*Ω₀ / g.a^2
        P ~ w*ρ
        δ ~ 0
        θ ~ 0
        cₛ² ~ 0
        σ ~ 0
    ]
    return System(eqs, τ, vars, pars; name, description, kwargs...)
end
