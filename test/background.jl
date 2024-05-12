import .Symboltz
using ModelingToolkit
using Plots; Plots.default(label=nothing)

@kwdef struct Parameters
    Ωr0 = 5.5e-5
    Ωm0 = 0.317
    Ωb0 = 0.05
    h = 0.67
    As = 2.1e-9
    Yp = 0.245
end
par = Parameters()

@named bg = Symboltz.background_ΛCDM()

if true
    bg_sol = Symboltz.solve(bg, par.Ωr0, par.Ωm0)

    # TODO: define plot() on BackgroundSolution
    p = plot(layout=(1,2), size=(1000, 400), margin=5*Plots.mm); display(p) # TODO: add plot recipe!
    plot!(p[1], bg_sol[Symboltz.η], bg_sol[bg.ssys.g.a]; xlabel="η / (1/H0)", ylabel="a", ylims=(0, 1)); display(p)
    plot!(p[2], log10.(bg_sol[bg.ssys.g.a]), stack(bg_sol[[bg.ssys.rad.ρ, bg.ssys.mat.ρ, bg.ssys.de.ρ]] ./ bg_sol[bg.sys.grav.ρcrit])'; xlabel="lg(a)", ylabel="Ω", label=["Ωr" "Ωm" "ΩΛ"], legend=:left); display(p)
end

@named th = Symboltz.thermodynamics_ΛCDM(bg)

if true
    th_sol = Symboltz.solve(th, par.Ωr0, par.Ωm0, par.Ωb0, par.h, par.Yp)

    # TODO: thermodynamics plot recipe
    p = plot(layout=(1,2), size=(1000, 400), margin=5*Plots.mm); display(p) # TODO: add plot recipe!
    plot!(p[1], log10.(th_sol[bg.ssys.g.a]), log10.(abs.(stack(th_sol[[Herec.Xe, Hrec.Xe, reion1.Xe, reion2.Xe, th.ssys.Xe]])')); xlabel="lg(a)", ylabel="lg(Xe)", ylims=(-5, +1), label=["XeS" "XeP" "XeRE1" "XeRE2" "Xe"], legend=:bottomleft); display(p)
    plot!(p[2], log10.(th_sol[bg.ssys.g.a]), log10.(stack(th_sol[[th.ssys.temp.Tγ, th.ssys.temp.Tb]])'); xlabel = "lg(a)", ylabel = "lg(T/K)", labels = ["Tγ" "Tb"]); display(p)
end

@named pt = Symboltz.perturbations_ΛCDM(th, 6)

if true
    ks = 10 .^ range(-4, +2, length=100) / Symboltz.k0 # in code units of k0 = H0/c
    pt_sols = Symboltz.solve(pt, ks, par.Ωr0, par.Ωm0, par.Ωb0, par.h, par.Yp; reltol=1e-8)

    p = plot(layout=(1,2), size=(1000, 400), margin=5*Plots.mm); display(p) # TODO: add plot recipe!
    for (i, pt_sol) in enumerate(pt_sols)
        plot!(p[1], log10.(pt_sol[pt.sys.bg.g.a]), pt_sol[pt.sys.gpt.Φ] / pt_sol[pt.sys.gpt.Φ][begin]; xlabel="lg(a)", ylabel="Φ/Φᵢ")
        plot!(p[2], log10.(pt_sol[pt.sys.bg.g.a]), log10.(abs.(pt_sol[pt.sys.cdm.δ])); color=i, xlabel="lg(a)", ylabel="lg(|δb|), lg(δc)")
        plot!(p[2], log10.(pt_sol[pt.sys.bg.g.a]), log10.(abs.(pt_sol[pt.sys.bar.δ])); color=i, xlabel="lg(a)", ylabel="lg(|δb|), lg(δc)")
    end
    display(p)
end