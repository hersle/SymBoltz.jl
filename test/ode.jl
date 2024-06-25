using ModelingToolkit
using DifferentialEquations
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

bgs = structural_simplify(bg)
ths = structural_simplify(th)
pts = structural_simplify(pt)

bg_prob = ODEProblem(bgs, [], (1e-5, 4.0), [bgs.ph.Ω0 => par.Ωγ0, bgs.neu.Ω0 => par.Ων0, bgs.cdm.Ω0 => par.Ωc0, bgs.bar.Ω0 => par.Ωb0, bgs.g.H0 => NaN])
th_prob = ODEProblem(ths, [], (1e-5, 4.0), [bg.ph.Ω0 => par.Ωγ0, bg.neu.Ω0 => par.Ων0, bg.cdm.Ω0 => par.Ωc0, bg.bar.Ω0 => par.Ωb0, bg.g.H0 => Symboltz.H100 * par.h, ths.Yp => par.Yp])
pt_prob = ODEProblem(pts, [], (1e-5, 4.0), [th.bg.ph.Ω0 => par.Ωγ0, th.bg.neu.Ω0 => par.Ων0, th.bg.cdm.Ω0 => par.Ωc0, th.bg.bar.Ω0 => par.Ωb0, th.bg.g.H0 => Symboltz.H100 * par.h, th.Yp => par.Yp, Symboltz.k => 1.0 / Symboltz.k0])

p = plot(layout=(3,2), size=(1000, 1200), margin=5*Plots.mm)

if true
    bg_sol = solve(bg_prob, Tsit5())
    plot!(p[1,1], bg_sol[Symboltz.t], bg_sol[bgs.g.a]; xlabel="t / (1/H0)", ylabel="a", ylims=(0, 1))
    plot!(p[1,2], log10.(bg_sol[bgs.g.a]), stack(bg_sol[[bgs.ph.ρ, bgs.neu.ρ, bgs.cdm.ρ, bgs.bar.ρ, bgs.de.ρ]] ./ bg_sol[bgs.grav.ρcrit])'; xlabel="lg(a)", ylabel="Ω", label=["Ωγ" "Ων" "Ωc" "Ωb" "ΩΛ"], legend=:left)
    display(p)
end

if true
    th_sol = solve(th_prob, Rodas5P(); reltol=1e-6)
    plot!(p[2,1], log10.(th_sol[bg.g.a]), log10.(abs.(stack(th_sol[[ths.Xe, ths.XH⁺, ths.XHe⁺, ths.XHe⁺⁺]])')); xlabel="lg(a)", ylabel="X", ylims=(-5, 1), label=["Xe" "XH⁺" "XHe⁺" "XHe⁺⁺"], legend=:bottomleft)
    plot!(p[2,2], log10.(th_sol[bg.g.a]), log10.(stack(th_sol[[ths.Tγ, ths.Tb]])'); xlabel = "lg(a)", ylabel = "lg(T/K)", labels = ["Tγ" "Tb"])
    display(p)
end

if true
    pt_sol = solve(pt_prob, Rodas5P(); reltol=1e-6)
    i = 1 #for (i, pt_sol) in enumerate(pt_sols)
        plot!(p[3,1], log10.(pt_sol[pt.th.bg.g.a]), pt_sol[pt.g1.Φ] / pt_sol[pt.g1.Φ][begin]; xlabel="lg(a)", ylabel="Φ/Φᵢ")
        plot!(p[3,2], log10.(pt_sol[pt.th.bg.g.a]), log10.(abs.(pt_sol[pt.cdm.δ])); color=i, xlabel="lg(a)", ylabel="lg(|δb|), lg(|δc|), lg(|δγ|), lg(|δν|)")
        plot!(p[3,2], log10.(pt_sol[pt.th.bg.g.a]), log10.(abs.(pt_sol[pt.bar.δ])); color=i, xlabel="lg(a)")
        plot!(p[3,2], log10.(pt_sol[pt.th.bg.g.a]), log10.(abs.(pt_sol[pt.ph.δ]));  color=i, xlabel="lg(a)")
        plot!(p[3,2], log10.(pt_sol[pt.th.bg.g.a]), log10.(abs.(pt_sol[pt.neu.δ])); color=i, xlabel="lg(a)")
    #end
    display(p)
end