using ModelingToolkit
using Plots
using Printf

@named M = SymBoltz.ΛCDM()
prob = SymBoltz.CosmologyProblem(M)
pars = [
    M.γ.Ω0 => 5.5e-5
    M.c.Ω0 => 0.267
    M.b.Ω0 => 0.05
    M.ν.Neff => 3.046
    M.g.h => 0.67
    M.b.rec.Yp => 0.245
]

ks = 10 .^ range(-1, 1, length=3) ./ SymBoltz.k0
sol = SymBoltz.solve(prob, pars, ks)

p = plot(layout=(3,3), size=(1600, 1200), margin=5*Plots.mm)

# plot background
plot!(p[1,1], sol, SymBoltz.t, M.g.a; ylims=(0, 1))
plot!(p[1,2], sol, log10(M.g.a), [M.γ.ρ, M.ν.ρ, M.h.ρ, M.c.ρ, M.b.ρ, M.Λ.ρ] ./ M.G.ρcrit; legend=:left, ylims=(0, 1))
plot!(p[2,1], sol, log10(M.g.a), [M.b.rec.Xe, M.b.rec.XH⁺, M.b.rec.XHe⁺, M.b.rec.XHe⁺⁺]; ylims=(0, 1.5))
plot!(p[2,2], sol, log10(M.g.a), log10.([M.b.rec.Tγ, M.b.rec.Tb]))

# plot perturbations
# TODO: color wavelengths like EM spectrum?
plot!(p[3,1], sol, ks, log10(M.g.a), [M.g.Φ, M.g.Ψ])
plot!(p[3,2], sol, ks, log10(M.g.a), log10.(abs.([M.γ.δ, M.c.δ, M.b.δ, M.ν.δ, M.h.δ#=, M.h.θ=#])); ylabel="log10(abs(δ(t))")
plot!(p[3,3], sol, ks, log10(M.g.a), [M.γ.θ, M.c.θ, M.b.θ, M.ν.θ] ./ M.k; ylabel="θ(t) / k")