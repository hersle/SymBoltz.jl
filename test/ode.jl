using ModelingToolkit
using Plots; Plots.default(label=nothing)
using Printf

# TODO: make plot recipes

model = Symboltz.ΛCDM()
par = Symboltz.CosmologicalParameters()

p = plot(layout=(3,3), size=(1600, 1200), margin=5*Plots.mm)

if true
    bg_sol = Symboltz.solve_background(model, par)
    plot!(p[1,1], bg_sol[Symboltz.t], bg_sol[model.bg_sim.g.a]; xlabel="t / (1/H0)", ylabel="a", ylims=(0, 1))
    plot!(p[1,2], log10.(bg_sol[model.bg_sim.g.a]), stack(bg_sol[[model.bg_sim.ph.ρ, model.bg_sim.neu.ρ, model.bg_sim.mneu.ρ, model.bg_sim.cdm.ρ, model.bg_sim.bar.ρ, model.bg_sim.de.ρ]] ./ bg_sol[model.bg_sim.grav.ρcrit])'; xlabel="lg(a)", ylabel="Ω", label=["Ω = Ωγ" "Ω = Ων" "Ω = Ωmν" "Ω = Ωc" "Ω = Ωb" "Ω = ΩΛ"], legend=:left)
    display(p)
end

if true
    th_sol = Symboltz.solve_thermodynamics(model, par)
    plot!(p[2,1], log10.(th_sol[model.th_sim.bg.g.a]), log10.(abs.(stack(th_sol[[model.th_sim.rec.Xe, model.th_sim.rec.XH⁺, model.th_sim.rec.XHe⁺, model.th_sim.rec.XHe⁺⁺]])')); xlabel="lg(a)", ylabel="X", ylims=(-5, 1), label=["X = Xe" "X = XH⁺" "X = XHe⁺" "X = XHe⁺⁺"], legend=:bottomleft)
    plot!(p[2,2], log10.(th_sol[model.th_sim.bg.g.a]), log10.(stack(th_sol[[model.th_sim.rec.Tγ, model.th_sim.rec.Tb]])'); xlabel = "lg(a)", ylabel = "lg(T/K)", labels = ["T = Tγ" "T = Tb"])
    display(p)
end

if true
    ks = 10 .^ range(-1, 1, length=3) ./ Symboltz.k0
    pt_sols = Symboltz.solve_perturbations(model, ks, par)
    for (i, (k, pt_sol)) in enumerate(zip(ks, pt_sols))
        color = i
        plot!(p[3,1], log10.(pt_sol[model.bg.g.a]), pt_sol[model.pt.g1.Φ]; linestyle=:solid, xlabel="lg(a)", color)
        plot!(p[3,1], log10.(pt_sol[model.bg.g.a]), pt_sol[model.pt.g1.Ψ]; linestyle=:dash,  xlabel="lg(a)", color)
        plot!(p[3,2], log10.(pt_sol[model.bg.g.a]), log10.(abs.(pt_sol[model.pt.ph.δ]));  linestyle=:solid, color)
        plot!(p[3,2], log10.(pt_sol[model.bg.g.a]), log10.(abs.(pt_sol[model.pt.cdm.δ])); linestyle=:dash, xlabel="lg(a)", ylabel="lg(|δ|)", color)
        plot!(p[3,2], log10.(pt_sol[model.bg.g.a]), log10.(abs.(pt_sol[model.pt.bar.δ])); linestyle=:dot, color)
        plot!(p[3,2], log10.(pt_sol[model.bg.g.a]), log10.(abs.(pt_sol[model.pt.neu.δ])); linestyle=:dashdot, color)
        plot!(p[3,2], log10.(pt_sol[model.bg.g.a]), log10.(abs.(pt_sol[model.pt.mneu.δ])); linestyle=:dashdotdot, color)
        plot!(p[3,3], log10.(pt_sol[model.bg.g.a]), pt_sol[model.pt.cdm.θ] ./ k; color, xlabel="lg(a)", ylabel="θ / k")
        plot!(p[3,3], log10.(pt_sol[model.bg.g.a]), pt_sol[model.pt.bar.θ] ./ k; color)
        plot!(p[3,3], log10.(pt_sol[model.bg.g.a]), pt_sol[model.pt.ph.θ] ./ k; color)
        plot!(p[3,3], log10.(pt_sol[model.bg.g.a]), pt_sol[model.pt.neu.θ] ./ k; color)
        #plot!(p[3,3], log10.(pt_sol[model.bg.g.a]), pt_sol[model.pt.mneu.θ]; label="θ = θmν")
    end
    hline!(p[3,1], [NaN NaN], linestyle=[:solid :dash], label=["Φ" "Ψ"], color=:black, legend_position=:topright)
    hline!(p[3,1], fill(NaN, 1, length(ks)), color=permutedims(eachindex(ks)), label=permutedims([(@sprintf "k = %f h/Mpc" k * Symboltz.k0) for k in ks]))
    hline!(p[3,2], [NaN NaN NaN NaN NaN], linestyle=[:solid :dash :dot :dashdot :dashdotdot], label="δ = δ" .* ["γ" "c" "b" "ν" "mν"], color=:black, legend_position=:topleft)
    hline!(p[3,3], [NaN NaN NaN NaN NaN], linestyle=[:solid :dash :dot :dashdot :dashdotdot], label="θ = θ" .* ["γ" "c" "b" "ν" "mν"], color=:black, legend_position=:topleft)
    display(p)
end