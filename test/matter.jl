using ModelingToolkit
using Plots; Plots.default(label=nothing)
using ForwardDiff, DiffResults, FiniteDiff

@named M = Symboltz.ΛCDM()
prob = Symboltz.CosmologyProblem(M)
ks = 10 .^ range(-3, 0, length=150) / Symboltz.k0
θ0 = [5.5e-5, 0.267, 0.05, 3.046, 0.67, 0.245]

function Pc(Ωγ0, Ωc0, Ωb0, Neff, h, Yp)
    pars = [M.γ.Ω0 => Ωγ0, M.c.Ω0 => Ωc0, M.b.Ω0 => Ωb0, M.ν.Neff => Neff, M.g.h => h, M.b.rec.Yp => Yp]
    return Symboltz.Pc(prob, pars, ks)
end

# computer power spectrum and derivatives wrt. input parameters using autodiff in one go
Pres = DiffResults.JacobianResult(ks, θ0)
log10Ph3(log10θ) = log10.(Pc((10 .^ log10θ)...) / Symboltz.k0^3)
ForwardDiff.jacobian!(Pres, log10Ph3, log10.(θ0))
lgPs, dlgP_dθs_ad = DiffResults.value(Pres), DiffResults.jacobian(Pres)

p = plot(layout=(2,1), size=(800, 1000), left_margin=bottom_margin=30*Plots.px); display(p)
plot!(p[1], log10.(ks*Symboltz.k0), lgPs; xlabel="lg(k/(h/Mpc))", ylabel="lg(P/(Mpc/h)³)"); display(p)
plot!(p[2], log10.(ks*Symboltz.k0), dlgP_dθs_ad; xlabel = "lg(k/(h/Mpc))", ylabel = "∂ lg(P) / ∂ lg(θᵢ)", labels = "θᵢ=" .* ["Ωγ0" "Ωm0" "Ωb0" "h" "Neff" "As" "Yp"]); display(p)

# TODO: fix finite differences
#dlgP_dθs_fd = FiniteDiff.finite_difference_jacobian(log10Ph3, log10.(θ0); relstep=1e-4) # relstep is important for finite difference accuracy!
#plot!(p[2], log10.(ks*Symboltz.k0), dlgP_dθs_fd; xlabel = "lg(k/(h/Mpc))", ylabel = "∂ lg(P) / ∂ lg(θᵢ)", color = :black); display(p)