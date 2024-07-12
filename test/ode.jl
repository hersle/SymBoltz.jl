using ModelingToolkit
using Plots; Plots.default(label=nothing)
using Printf

# TODO: make plot recipes

@named M = Symboltz.ΛCDM()
M = Symboltz.CosmologicalModel(M)
par = Symboltz.CosmologicalParameters()

p = plot(layout=(3,3), size=(1600, 1200), margin=5*Plots.mm)

if true
    bg_sol = Symboltz.solve_background(M, par)
    plot!(p[1,1], bg_sol[Symboltz.t], bg_sol[M.bg_sim.g.a]; xlabel="t / (1/H0)", ylabel="a", ylims=(0, 1))
    plot!(p[1,2], log10.(bg_sol[M.bg_sim.g.a]), stack(bg_sol[[M.bg_sim.γ.ρ, M.bg_sim.ν.ρ, #=M.bg_sim.mneu.ρ,=# M.bg_sim.c.ρ, M.bg_sim.b.ρ, M.bg_sim.Λ.ρ]] ./ bg_sol[M.bg_sim.G.ρcrit])'; xlabel="lg(a)", ylabel="Ω", label=["Ω = Ωγ" "Ω = Ων" #="Ω = Ωmν"=# "Ω = Ωc" "Ω = Ωb" "Ω = ΩΛ"], legend=:left)
    plot!(p[2,1], log10.(bg_sol[M.bg_sim.g.a]), stack(bg_sol[[M.bg_sim.b.rec.Xe, M.bg_sim.b.rec.XH⁺, M.bg_sim.b.rec.XHe⁺, M.bg_sim.b.rec.XHe⁺⁺]])'; xlabel="lg(a)", ylabel="X", ylims=(0, 1.5), label=["X = Xe" "X = XH⁺" "X = XHe⁺" "X = XHe⁺⁺"])
    plot!(p[2,2], log10.(bg_sol[M.bg_sim.g.a]), log10.(stack(bg_sol[[M.bg_sim.b.rec.Tγ, M.bg_sim.b.rec.Tb]])'); xlabel = "lg(a)", ylabel = "lg(T/K)", labels = ["T = Tγ" "T = Tb"])
    display(p)
end

# TODO: color wavelengths like EM spectrum
if true
    ks = 10 .^ range(-1, 1, length=3) ./ Symboltz.k0
    pt_sols = Symboltz.solve_perturbations(M, ks, par)
    for (i, (k, pt_sol)) in enumerate(zip(ks, pt_sols))
        color = i
        plot!(p[3,1], log10.(pt_sol[M.bg.g.a]), pt_sol[M.pt.g.Φ]; linestyle=:solid, xlabel="lg(a)", color)
        plot!(p[3,1], log10.(pt_sol[M.bg.g.a]), pt_sol[M.pt.g.Ψ]; linestyle=:dash,  xlabel="lg(a)", color)
        plot!(p[3,2], log10.(pt_sol[M.bg.g.a]), log10.(abs.(pt_sol[M.pt.γ.δ]));  linestyle=:solid, color)
        plot!(p[3,2], log10.(pt_sol[M.bg.g.a]), log10.(abs.(pt_sol[M.pt.c.δ])); linestyle=:dash, xlabel="lg(a)", ylabel="lg(|δ|)", color)
        plot!(p[3,2], log10.(pt_sol[M.bg.g.a]), log10.(abs.(pt_sol[M.pt.b.δ])); linestyle=:dot, color)
        plot!(p[3,2], log10.(pt_sol[M.bg.g.a]), log10.(abs.(pt_sol[M.pt.ν.δ])); linestyle=:dashdot, color)
        #plot!(p[3,2], log10.(pt_sol[M.bg.g.a]), log10.(abs.(pt_sol[M.pt.mneu.δ])); linestyle=:dashdotdot, color)
        plot!(p[3,3], log10.(pt_sol[M.bg.g.a]), pt_sol[M.pt.γ.θ] ./ k; color, xlabel="lg(a)", ylabel="θ / k")
        plot!(p[3,3], log10.(pt_sol[M.bg.g.a]), pt_sol[M.pt.c.θ] ./ k; color)
        plot!(p[3,3], log10.(pt_sol[M.bg.g.a]), pt_sol[M.pt.b.θ] ./ k; color)
        plot!(p[3,3], log10.(pt_sol[M.bg.g.a]), pt_sol[M.pt.ν.θ] ./ k; color)
        #plot!(p[3,3], log10.(pt_sol[M.bg.g.a]), pt_sol[M.pt.mneu.θ]; label="θ = θmν")
    end
    hline!(p[3,1], [NaN NaN], linestyle=[:solid :dash], label=["Φ" "Ψ"], color=:black, legend_position=:topright)
    hline!(p[3,1], fill(NaN, 1, length(ks)), color=permutedims(eachindex(ks)), label=permutedims([(@sprintf "k = %f h/Mpc" k * Symboltz.k0) for k in ks]))
    hline!(p[3,2], [NaN NaN NaN NaN #=NaN=#], linestyle=[:solid :dash :dot :dashdot :dashdotdot], label="δ = δ" .* ["γ" "c" "b" "ν" #="mν"=#], color=:black, legend_position=:topleft)
    hline!(p[3,3], [NaN NaN NaN NaN #=NaN=#], linestyle=[:solid :dash :dot :dashdot :dashdotdot], label="θ = θ" .* ["γ" "c" "b" "ν" #="mν"=#], color=:black, legend_position=:topleft)
    display(p)
end