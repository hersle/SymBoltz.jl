import .Symboltz
using ModelingToolkit
using Plots; Plots.default(label=nothing)
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
@named bg = Symboltz.background_ΛCDM()
@named th = Symboltz.thermodynamics_ΛCDM(bg)
@named pt = Symboltz.perturbations_ΛCDM(th, 6)

ks = 10 .^ range(-4, +2, length=300) / Symboltz.k0
θ0 = [par.Ωr0, par.Ωm0, par.Ωb0, par.h, par.As, par.Yp]

# computer power spectrum and derivatives wrt. input parameters using autodiff in one go
Pres = DiffResults.JacobianResult(ks, θ0)
log10Ph3(log10θ) = log10.(Symboltz.P(pt, ks, (10 .^ log10θ)...) / Symboltz.k0^3)
ForwardDiff.jacobian!(Pres, log10Ph3, log10.(θ0))
lgPs, dlgP_dθs_ad = DiffResults.value(Pres), DiffResults.jacobian(Pres)

p = plot(layout=(2,1), size=(800, 1000), left_margin=bottom_margin=30*Plots.px); display(p)
plot!(p[1], log10.(ks*Symboltz.k0), lgPs; xlabel="lg(k/(h/Mpc))", label="lg(P/(Mpc/h)³)"); display(p)
plot!(p[2], log10.(ks*Symboltz.k0), dlgP_dθs_ad; xlabel = "lg(k/(h/Mpc))", ylabel = "∂ lg(P) / ∂ lg(θᵢ)", labels = "θᵢ=" .* ["Ωr0" "Ωm0" "Ωb0" "h" "As" "Yp"]); display(p)

# TODO: fix finite differences
#dlgP_dθs_fd = FiniteDiff.finite_difference_jacobian(log10Ph3, log10.(θ0); relstep=1e-4) # relstep is important for finite difference accuracy!
#plot!(p[2], log10.(ks*Symboltz.k0), dlgP_dθs_fd; xlabel = "lg(k/(h/Mpc))", ylabel = "∂ lg(P) / ∂ lg(θᵢ)", color = :black); display(p)