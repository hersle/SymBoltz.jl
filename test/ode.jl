using ModelingToolkit
using Plots; Plots.default(label=nothing)

# TODO: make plot recipes

model = Symboltz.ΛCDM()
par = Symboltz.CosmologicalParameters()

p = plot(layout=(3,2), size=(1000, 1200), margin=5*Plots.mm)

if true
    bg_sol = Symboltz.solve_background(model, par)
    plot!(p[1,1], bg_sol[Symboltz.t], bg_sol[model.bg_sim.g.a]; xlabel="t / (1/H0)", ylabel="a", ylims=(0, 1))
    plot!(p[1,2], log10.(bg_sol[model.bg_sim.g.a]), stack(bg_sol[[model.bg_sim.ph.ρ, model.bg_sim.neu.ρ, model.bg_sim.cdm.ρ, model.bg_sim.bar.ρ, model.bg_sim.de.ρ]] ./ bg_sol[model.bg_sim.grav.ρcrit])'; xlabel="lg(a)", ylabel="Ω", label=["Ωγ" "Ων" "Ωc" "Ωb" "ΩΛ"], legend=:left)
    display(p)
end

if true
    th_sol = Symboltz.solve_thermodynamics(model, par)
    plot!(p[2,1], log10.(th_sol[model.bg.g.a]), log10.(abs.(stack(th_sol[[model.th_sim.Xe, model.th_sim.XH⁺, model.th_sim.XHe⁺, model.th_sim.XHe⁺⁺]])')); xlabel="lg(a)", ylabel="X", ylims=(-5, 1), label=["Xe" "XH⁺" "XHe⁺" "XHe⁺⁺"], legend=:bottomleft)
    plot!(p[2,2], log10.(th_sol[model.bg.g.a]), log10.(stack(th_sol[[model.th_sim.Tγ, model.th_sim.Tb]])'); xlabel = "lg(a)", ylabel = "lg(T/K)", labels = ["Tγ" "Tb"])
    display(p)
end

if true
    pt_sol = Symboltz.solve_perturbations(model, 1.0 / Symboltz.k0, par)
    i = 1 #for (i, pt_sol) in enumerate(pt_sols)
        plot!(p[3,1], log10.(pt_sol[model.bg.g.a]), pt_sol[model.pt.g1.Φ] / pt_sol[model.pt.g1.Φ][begin]; xlabel="lg(a)", ylabel="Φ/Φᵢ")
        plot!(p[3,2], log10.(pt_sol[model.bg.g.a]), log10.(abs.(pt_sol[model.pt.cdm.δ])); color=i, xlabel="lg(a)", ylabel="lg(|δb|), lg(|δc|), lg(|δγ|), lg(|δν|)")
        plot!(p[3,2], log10.(pt_sol[model.bg.g.a]), log10.(abs.(pt_sol[model.pt.bar.δ])); color=i, xlabel="lg(a)")
        plot!(p[3,2], log10.(pt_sol[model.bg.g.a]), log10.(abs.(pt_sol[model.pt.ph.δ]));  color=i, xlabel="lg(a)")
        plot!(p[3,2], log10.(pt_sol[model.bg.g.a]), log10.(abs.(pt_sol[model.pt.neu.δ])); color=i, xlabel="lg(a)")
    #end
    display(p)
end