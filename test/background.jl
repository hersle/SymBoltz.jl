using .Symboltz
using .Symboltz: η
using ModelingToolkit
using Plots; Plots.default(label=nothing, markershape=:pixel)
using ForwardDiff, DiffResults, FiniteDiff

@named g = Symboltz.background_metric()
@named grav = Symboltz.background_gravity_GR(g)
@named rad = Symboltz.background_radiation(g)
@named mat = Symboltz.background_matter(g)
@named de = Symboltz.background_cosmological_constant(g)
@named bg = BackgroundSystem(g, grav, [rad, mat, de])

Ωr0 = 5.5e-5
Ωm0 = 0.317
bg_sol = solve(bg, Ωr0, Ωm0)
println("ηi = $(bg_sol[η][begin]), η0 = $(bg_sol[η][end])")

# TODO: define plot() on BackgroundSolution
p = plot(layout=(1,2), size=(1200, 400), left_margin=bottom_margin=30*Plots.px); display(p) # TODO: add plot recipe!
plot!(p[1], bg_sol[η], bg_sol[bg.ssys.g.a]; xlabel="η / (1/H0)", ylabel="a", ylims=(0, 1)); display(p)
plot!(p[2], log10.(bg_sol[bg.ssys.g.a]), stack(bg_sol[[bg.ssys.rad.Ω,bg.ssys.mat.Ω,bg.ssys.de.Ω]])'; xlabel="lg(a)", ylabel="Ω", label=["Ωr" "Ωm" "ΩΛ"], legend=:left); display(p)