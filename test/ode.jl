using ModelingToolkit
using Plots; Plots.default(label=nothing)

# TODO: make plot recipes

model = Symboltz.ΛCDM()
par = Symboltz.CosmologicalParameters()

p = plot(layout=(3,2), size=(1000, 1200), margin=5*Plots.mm)

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
    pt_sol = Symboltz.solve_perturbations(model, 1.0 / Symboltz.k0, par)
    plot!(p[3,1], log10.(pt_sol[model.bg.g.a]), pt_sol[model.pt.g1.Φ] / pt_sol[model.pt.g1.Φ][begin]; xlabel="lg(a)", ylabel="Φ/Φᵢ")
    plot!(p[3,2], log10.(pt_sol[model.bg.g.a]), log10.(abs.(pt_sol[model.pt.cdm.δ])); label="δ = δc", xlabel="lg(a)", ylabel="lg(|δ|)")
    plot!(p[3,2], log10.(pt_sol[model.bg.g.a]), log10.(abs.(pt_sol[model.pt.bar.δ])); label="δ = δb", xlabel="lg(a)")
    plot!(p[3,2], log10.(pt_sol[model.bg.g.a]), log10.(abs.(pt_sol[model.pt.ph.δ]));  label="δ = δγ", xlabel="lg(a)")
    plot!(p[3,2], log10.(pt_sol[model.bg.g.a]), log10.(abs.(pt_sol[model.pt.neu.δ])); label="δ = δν", xlabel="lg(a)")
    plot!(p[3,2], log10.(pt_sol[model.bg.g.a]), log10.(abs.(pt_sol[model.pt.mneu.δ])); label="δ = δmν", xlabel="lg(a)")
    display(p)
end