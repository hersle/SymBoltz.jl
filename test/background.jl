using .Symboltz
using .Symboltz: η
using ModelingToolkit
using Plots; Plots.default(label=nothing, markershape=:pixel)
using ForwardDiff, DiffResults, FiniteDiff

@kwdef struct Parameters
    Ωr0 = 5.5e-5
    Ωm0 = 0.317
    Ωb0 = 0.05
    h = 0.67
    As = 2.1e-9
    Yp = 0.245
end
par = Parameters()

@named g = Symboltz.background_metric()
@named grav = Symboltz.background_gravity_GR(g)
@named rad = Symboltz.background_radiation(g)
@named mat = Symboltz.background_matter(g)
@named de = Symboltz.background_cosmological_constant(g)
@named bg = BackgroundSystem(g, grav, [rad, mat, de])

if false
    bg_sol = solve(bg, par.Ωr0, par.Ωm0)
    println("ηi = $(bg_sol[η][begin]), η0 = $(bg_sol[η][end])")

    # TODO: define plot() on BackgroundSolution
    p = plot(layout=(1,2), size=(1000, 400), margin=5*Plots.mm); display(p) # TODO: add plot recipe!
    plot!(p[1], bg_sol[η], bg_sol[bg.ssys.g.a]; xlabel="η / (1/H0)", ylabel="a", ylims=(0, 1)); display(p)
    plot!(p[2], log10.(bg_sol[bg.ssys.g.a]), stack(bg_sol[[bg.ssys.rad.Ω, bg.ssys.mat.Ω, bg.ssys.de.Ω]])'; xlabel="lg(a)", ylabel="Ω", label=["Ωr" "Ωm" "ΩΛ"], legend=:left); display(p)
end

@named temp = Symboltz.thermodynamics_temperature(g)
@named Herec = Symboltz.recombination_helium_saha()
@named Hrec = Symboltz.recombination_hydrogen_peebles(g)
@named reion1 = Symboltz.reionization_smooth_step(g, 8.0, 0.5, 1 + Herec.Yp/(4*(1-Herec.Yp))) # TODO: separate reionH₊, reionHe₊, reionHe₊₊
@named reion2 = Symboltz.reionization_smooth_step(g, 3.5, 0.5, Herec.Yp/(4*(1-Herec.Yp)))
@named th = ThermodynamicsSystem(bg, Herec, Hrec, temp, [reion1, reion2])

if false
    th_sol = solve(th, par.Ωr0, par.Ωm0, par.Ωb0, par.h, par.Yp)

    # TODO: thermodynamics plot recipe
    p = plot(layout=(1,2), size=(1000, 400), margin=5*Plots.mm); display(p) # TODO: add plot recipe!
    plot!(p[1], log10.(th_sol[bg.ssys.g.a]), log10.(abs.(stack(th_sol[[Herec.Xe, Hrec.Xe, reion1.Xe, reion2.Xe, th.ssys.Xe]])')); xlabel="lg(a)", ylabel="lg(Xe)", ylims=(-5, +1), label=["XeS" "XeP" "XeRE1" "XeRE2" "Xe"], legend=:bottomleft); display(p)
    plot!(p[2], log10.(th_sol[bg.ssys.g.a]), log10.(stack(th_sol[[th.ssys.temp.Tγ, th.ssys.temp.Tb]])'); xlabel = "lg(a)", ylabel = "lg(T/K)", labels = ["Tγ" "Tb"]); display(p)
end

lmax = 6
@named gpt = Symboltz.perturbations_metric()
@named ph = Symboltz.perturbations_photon_hierarchy(gpt, lmax, true)
@named pol = Symboltz.perturbations_polarization_hierarchy(gpt, lmax)
@named cdm = Symboltz.perturbations_matter(g, gpt, false)
@named bar = Symboltz.perturbations_matter(g, gpt, true)
@named gravpt = Symboltz.perturbations_gravity(g, gpt, )
@named pt = PerturbationsSystem(bg, th, gpt, gravpt, ph, pol, cdm, bar)

if true
    ks = 10 .^ range(-4, +2, length=100) / Symboltz.k0 # in code units of k0 = H0/c
    pt_sols = solve(pt, ks, par.Ωr0, par.Ωm0, par.Ωb0, par.h, par.Yp; reltol=1e-8)

    p = plot(layout=(1,2), size=(1000, 400), margin=5*Plots.mm); display(p) # TODO: add plot recipe!
    for (i, pt_sol) in enumerate(pt_sols)
        plot!(p[1], log10.(pt_sol[pt.sys.bg.g.a]), pt_sol[pt.sys.gpt.Φ] / pt_sol[pt.sys.gpt.Φ][begin]; xlabel="lg(a)", ylabel="Φ/Φᵢ")
        plot!(p[2], log10.(pt_sol[pt.sys.bg.g.a]), log10.(abs.(pt_sol[pt.sys.cdm.δ])); color=i, xlabel="lg(a)", ylabel="lg(|δb|), lg(δc)")
        plot!(p[2], log10.(pt_sol[pt.sys.bg.g.a]), log10.(abs.(pt_sol[pt.sys.bar.δ])); color=i, xlabel="lg(a)", ylabel="lg(|δb|), lg(δc)")
    end
    display(p)
end