using ModelingToolkit
using Plots; Plots.default(label=nothing)
using ForwardDiff, DiffResults, FiniteDiff

M = SymBoltz.ΛCDM()
prob = SymBoltz.CosmologyProblem(M)
ks = 10 .^ range(-3, 0, length=150) / SymBoltz.k0
θ0 = [5.5e-5, 0.267, 0.05, 3.046, 0.67, 0.245]

P(θ) = SymBoltz.P(prob, M.c, [M.γ.Ω0, M.c.Ω0, M.b.Ω0, M.ν.Neff, M.g.h, M.b.rec.Yp] .=> θ, ks) / SymBoltz.k0^3
lgP(lgθ) = log10.(P(10 .^ lgθ))

# computer power spectrum and derivatives wrt. input parameters using autodiff in one go
Pres = DiffResults.JacobianResult(ks, θ0)
ForwardDiff.jacobian!(Pres, lgP, log10.(θ0))
lgPs, dlgP_dθs_ad = DiffResults.value(Pres), DiffResults.jacobian(Pres)

p = plot(layout=(2,1), size=(800, 1000), left_margin=bottom_margin=30*Plots.px); display(p)
plot!(p[1], log10.(ks*SymBoltz.k0), lgPs; xlabel="lg(k/(h/Mpc))", ylabel="lg(P/(Mpc/h)³)"); display(p)
plot!(p[2], log10.(ks*SymBoltz.k0), dlgP_dθs_ad; xlabel = "lg(k/(h/Mpc))", ylabel = "∂ lg(P) / ∂ lg(θᵢ)", labels = "θᵢ=" .* ["Ωγ0" "Ωm0" "Ωb0" "h" "Neff" "As" "Yp"]); display(p)

# TODO: fix finite differences
#dlgP_dθs_fd = FiniteDiff.finite_difference_jacobian(log10Ph3, log10.(θ0); relstep=1e-4) # relstep is important for finite difference accuracy!
#plot!(p[2], log10.(ks*SymBoltz.k0), dlgP_dθs_fd; xlabel = "lg(k/(h/Mpc))", ylabel = "∂ lg(P) / ∂ lg(θᵢ)", color = :black); display(p)