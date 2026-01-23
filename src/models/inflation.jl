function harrison_zeldovich(g; kwargs...)
    pars = @parameters begin
        ln_As1e10 = NaN, [description = "ln(10¹⁰ As)"] # keep uninitialized if not needed
        As = exp(ln_As1e10) / 1e10, [description = "Spectral amplitude"]
        ns = NaN, [description = "Spectral index"]
        kpivot = 0.05 / Mpc / (H100/c) / g.h, [description = "Pivot scale wavenumber"] # k = 0.05/Mpc ≠ 0.05/(Mpc/h)
    end
    vars = @variables P(τ, k)
    eqs = [
        P ~ 2*Num(π)^2 / k^3 * As * (k/kpivot)^(ns-1)
    ]
    description = "Harrison-Zel'dovich inflation"
    return System(eqs, τ, vars, pars; description, kwargs...)
end
