function harrison_zeldovich(g; kwargs...)
    @parameters begin
        As, [description = "Spectral amplitude"]
        ns, [description = "Spectral index"]
        kpivot, [description = "Pivot scale wavenumber"]
    end
    @variables P(τ)
    eqs = [
        P ~ 2*Num(π)^2 / k^3 * As * (k/kpivot)^(ns-1)
    ]
    defaults = [
        kpivot => 0.05 / (H100*g.h/c) / Mpc # k = 0.05/Mpc ≠ 0.05/(Mpc/h)
        As => NaN # keep uninitialized if not needed
        ns => NaN
    ]
    description = "Harrison-Zel'dovich inflation"
    return ODESystem(eqs, τ; defaults, description, kwargs...)
end
