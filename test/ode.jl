using ModelingToolkit
using Plots; Plots.default(label=nothing)
using Printf

# TODO: make plot recipes

@named M = Symboltz.ΛCDM()
prob = Symboltz.CosmologyProblem(M)
pars = [
    M.γ.Ω0 => 5.5e-5
    M.c.Ω0 => 0.267
    M.b.Ω0 => 0.05
    M.ν.Neff => 3.046
    M.g.h => 0.67
    M.b.rec.Yp => 0.245
]

p = plot(layout=(3,3), size=(1600, 1200), margin=5*Plots.mm)

if true
    bg_sol = Symboltz.solve(prob, pars)
    plot!(p[1,1], bg_sol[Symboltz.t], bg_sol[M.g.a]; xlabel="t / (1/H0)", ylabel="a", ylims=(0, 1))
    plot!(p[1,2], log10.(bg_sol[M.g.a]), stack(bg_sol[[M.γ.ρ, M.ν.ρ, M.h.ρ, M.c.ρ, M.b.ρ, M.Λ.ρ]] ./ bg_sol[M.G.ρcrit])'; xlabel="lg(a)", ylabel="Ω", label=["Ω = Ωγ" "Ω = Ων" "Ω = Ωh" "Ω = Ωc" "Ω = Ωb" "Ω = ΩΛ"], legend=:left)
    plot!(p[2,1], log10.(bg_sol[M.g.a]), stack(bg_sol[[M.b.rec.Xe, M.b.rec.XH⁺, M.b.rec.XHe⁺, M.b.rec.XHe⁺⁺]])'; xlabel="lg(a)", ylabel="X", ylims=(0, 1.5), label=["X = Xe" "X = XH⁺" "X = XHe⁺" "X = XHe⁺⁺"])
    plot!(p[2,2], log10.(bg_sol[M.g.a]), log10.(stack(bg_sol[[M.b.rec.Tγ, M.b.rec.Tb]])'); xlabel = "lg(a)", ylabel = "lg(T/K)", labels = ["T = Tγ" "T = Tb"])
    display(p)
end

# TODO: color wavelengths like EM spectrum
if true
    ks = 10 .^ range(-1, 1, length=3) ./ Symboltz.k0
    for (i, k) in enumerate(ks)
        pt_sol = Symboltz.solve(prob, pars, k)
        color = i
        plot!(p[3,1], log10.(pt_sol.pt[M.g.a]), pt_sol[M.g.Φ]; linestyle=:solid, xlabel="lg(a)", color)
        plot!(p[3,1], log10.(pt_sol.pt[M.g.a]), pt_sol[M.g.Ψ]; linestyle=:dash,  xlabel="lg(a)", color)
        plot!(p[3,2], log10.(pt_sol.pt[M.g.a]), log10.(abs.(pt_sol[M.γ.δ]));  linestyle=:solid, color)
        plot!(p[3,2], log10.(pt_sol.pt[M.g.a]), log10.(abs.(pt_sol[M.c.δ])); linestyle=:dash, xlabel="lg(a)", ylabel="lg(|δ|)", color)
        plot!(p[3,2], log10.(pt_sol.pt[M.g.a]), log10.(abs.(pt_sol[M.b.δ])); linestyle=:dot, color)
        plot!(p[3,2], log10.(pt_sol.pt[M.g.a]), log10.(abs.(pt_sol[M.ν.δ])); linestyle=:dashdot, color)
        plot!(p[3,2], log10.(pt_sol.pt[M.g.a]), log10.(abs.(pt_sol[M.h.δ])); linestyle=:dashdotdot, color)
        plot!(p[3,3], log10.(pt_sol.pt[M.g.a]), pt_sol[M.γ.θ] ./ k; color, xlabel="lg(a)", ylabel="θ / k")
        plot!(p[3,3], log10.(pt_sol.pt[M.g.a]), pt_sol[M.c.θ] ./ k; color)
        plot!(p[3,3], log10.(pt_sol.pt[M.g.a]), pt_sol[M.b.θ] ./ k; color)
        plot!(p[3,3], log10.(pt_sol.pt[M.g.a]), pt_sol[M.ν.θ] ./ k; color)
        #plot!(p[3,3], log10.(pt_sol[M.g.a]), pt_sol[M.h.θ] ./ k; color) # TODO: compute
    end
    hline!(p[3,1], [NaN NaN], linestyle=[:solid :dash], label=["Φ" "Ψ"], color=:black, legend_position=:topright)
    hline!(p[3,1], fill(NaN, 1, length(ks)), color=permutedims(eachindex(ks)), label=permutedims([(@sprintf "k = %f h/Mpc" k * Symboltz.k0) for k in ks]))
    hline!(p[3,2], [NaN NaN NaN NaN NaN], linestyle=[:solid :dash :dot :dashdot :dashdotdot], label="δ = δ" .* ["γ" "c" "b" "ν" "h"], color=:black, legend_position=:topleft)
    hline!(p[3,3], [NaN NaN NaN NaN #=NaN=#], linestyle=[:solid :dash :dot :dashdot :dashdotdot], label="θ = θ" .* ["γ" "c" "b" "ν" #="h"=#], color=:black, legend_position=:topleft)
    display(p)
end