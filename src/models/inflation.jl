function harrison_zeldovich(g; kwargs...)
    @parameters begin
        ln_As1e10, [description = "ln(10¹⁰ As)"]
        As, [description = "Spectral amplitude"]
        ns, [description = "Spectral index"]
        kpivot, [description = "Pivot scale wavenumber"]
    end
    @variables P(τ, k)
    eqs = [
        As ~ exp(ln_As1e10) / 1e10
        P ~ 2*Num(π)^2 / k^3 * As * (k/kpivot)^(ns-1)
    ]
    defaults = [
        kpivot => 0.05 / (H100*g.h/c) / Mpc # k = 0.05/Mpc ≠ 0.05/(Mpc/h)
        ln_As1e10 => NaN # keep uninitialized if not needed
        ns => NaN
    ]
    description = "Harrison-Zel'dovich inflation"
    return System(eqs, τ; defaults, description, kwargs...)
end
