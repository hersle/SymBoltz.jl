using NumericalIntegration
using Bessels: besselj, sphericalbesselj

# matter power spectrum
P0(k, As) = @. As / k ^ 3
function P(pt::PerturbationsSystem, k, Ωr0, Ωm0, Ωb0, h, As, Yp)
    pt_sols = solve(pt, k, Ωr0, Ωm0, Ωb0, h, Yp)
    ηtoday = pt_sols[1].prob.tspan[end] # TODO: something more robust?
    return P0(k, As) .* pt_sols(ηtoday, idxs=pt.sys.gravpt.Δm) .^ 2
end