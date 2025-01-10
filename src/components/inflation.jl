function harrison_zeldovich(g; kwargs...)
    @parameters As ns kpivot
    @variables P(t)
    eqs = [
        P ~ 2*Num(π)^2 / k^3 * As * (k/kpivot)^(ns-1)
    ]
    defaults = [
        kpivot => 0.05 / (g.H₀/c) / Mpc
    ]
    description = "Harrison-Zel'dovich inflation"
    return ODESystem(eqs, t; defaults, description, kwargs...)
end
