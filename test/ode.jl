include("../Symboltz.jl")
import .Symboltz
using ModelingToolkit
using Plots; Plots.default(label=nothing)

# TODO: make plot recipes

@kwdef struct Parameters
    Ωγ0 = 5.5e-5
    Ων0 = 3.046 * 7/8 * (4/11)^(4/3) * 5.5e-5 # TODO: handle more elegantly with Neff/Tν0
    Ωc0 = 0.267
    Ωb0 = 0.05
    h = 0.67
    As = 2.1e-9
    Yp = 0.245
end

par = Parameters()
@named bg = Symboltz.background_ΛCDM()
@named th = Symboltz.thermodynamics_ΛCDM(bg)
@named pt = Symboltz.perturbations_ΛCDM(th, 6)

p = plot(layout=(3,2), size=(1000, 1200), margin=5*Plots.mm)

if true
    bg_sol = Symboltz.solve(bg, par.Ωγ0, par.Ων0, par.Ωc0, par.Ωb0)
    plot!(p[1,1], bg_sol[Symboltz.η], bg_sol[bg.ssys.g.a]; xlabel="η / (1/H0)", ylabel="a", ylims=(0, 1))
    plot!(p[1,2], log10.(bg_sol[bg.ssys.g.a]), stack(bg_sol[[bg.ssys.ph.ρ, bg.ssys.neu.ρ, bg.ssys.cdm.ρ, bg.ssys.bar.ρ, bg.ssys.de.ρ]] ./ bg_sol[bg.sys.grav.ρcrit])'; xlabel="lg(a)", ylabel="Ω", label=["Ωγ" "Ων" "Ωc" "Ωb" "ΩΛ"], legend=:left)
    display(p)
end

if true
    th_sol = Symboltz.solve(th, par.Ωγ0, par.Ων0, par.Ωc0, par.Ωb0, par.h, par.Yp)
    plot!(p[2,1], log10.(th_sol[bg.ssys.g.a]), log10.(abs.(stack(th_sol[[th.ssys.Xe, th.ssys.Xp]])')); xlabel="lg(a)", ylabel="X", ylims=(-5, 1), label=["Xe" "Xp"], legend=:bottomleft)
    plot!(p[2,2], log10.(th_sol[bg.ssys.g.a]), log10.(stack(th_sol[[th.sys.Tγ, th.sys.Tb]])'); xlabel = "lg(a)", ylabel = "lg(T/K)", labels = ["Tγ" "Tb"])
    display(p)
end

if true
    ks = 10 .^ range(-4, +2, length=100) / Symboltz.k0 # in code units of k0 = H0/c
    pt_sols = Symboltz.solve(pt, ks, par.Ωγ0, par.Ων0, par.Ωc0, par.Ωb0, par.h, par.Yp; reltol=1e-8)
    for (i, pt_sol) in enumerate(pt_sols)
        plot!(p[3,1], log10.(pt_sol[pt.sys.bg.g.a]), pt_sol[pt.sys.g1.Φ] / pt_sol[pt.sys.g1.Φ][begin]; xlabel="lg(a)", ylabel="Φ/Φᵢ")
        plot!(p[3,2], log10.(pt_sol[pt.sys.bg.g.a]), log10.(abs.(pt_sol[pt.sys.cdm.δ])); color=i, xlabel="lg(a)", ylabel="lg(|δb|), lg(|δc|), lg(|δγ|), lg(|δν|)")
        plot!(p[3,2], log10.(pt_sol[pt.sys.bg.g.a]), log10.(abs.(pt_sol[pt.sys.bar.δ])); color=i, xlabel="lg(a)")
        plot!(p[3,2], log10.(pt_sol[pt.sys.bg.g.a]), log10.(abs.(pt_sol[pt.sys.ph.δ]));  color=i, xlabel="lg(a)")
        plot!(p[3,2], log10.(pt_sol[pt.sys.bg.g.a]), log10.(abs.(pt_sol[pt.sys.neu.δ])); color=i, xlabel="lg(a)")
    end
    display(p)
end