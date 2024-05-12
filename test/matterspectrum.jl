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

@named g = Symboltz.background_metric()
@named grav = Symboltz.background_gravity_GR(g)
@named rad = Symboltz.background_radiation(g)
@named mat = Symboltz.background_matter(g)
@named de = Symboltz.background_cosmological_constant(g)
@named bg = Symboltz.BackgroundSystem(g, grav, [rad, mat, de])

@named temp = Symboltz.thermodynamics_temperature(g)
@named Herec = Symboltz.recombination_helium_saha()
@named Hrec = Symboltz.recombination_hydrogen_peebles(g)
@named reion1 = Symboltz.reionization_smooth_step(g, 8.0, 0.5, 1 + Herec.Yp/(4*(1-Herec.Yp))) # TODO: separate reionH₊, reionHe₊, reionHe₊₊
@named reion2 = Symboltz.reionization_smooth_step(g, 3.5, 0.5, Herec.Yp/(4*(1-Herec.Yp)))
@named th = Symboltz.ThermodynamicsSystem(bg, Herec, Hrec, temp, [reion1, reion2])

lmax = 6
@named gpt = Symboltz.perturbations_metric()
@named ph = Symboltz.perturbations_photon_hierarchy(gpt, lmax, true)
@named pol = Symboltz.perturbations_polarization_hierarchy(gpt, lmax)
@named cdm = Symboltz.perturbations_matter(g, gpt, false)
@named bar = Symboltz.perturbations_matter(g, gpt, true)
@named gravpt = Symboltz.perturbations_gravity(g, gpt)
@named pt = Symboltz.PerturbationsSystem(bg, th, gpt, gravpt, ph, pol, cdm, bar)

ks = 10 .^ range(-4, +2, length=300) / Symboltz.k0
θ0 = [par.Ωr0, par.Ωm0, par.Ωb0, par.h, par.As, par.Yp]

# computer power spectrum and derivatives wrt. input parameters using autodiff in one go
Pres = DiffResults.JacobianResult(ks, θ0)
log10Ph3(log10θ) = log10.(Symboltz.P(pt, ks, (10 .^ log10θ)...) / Symboltz.k0^3)
ForwardDiff.jacobian!(Pres, log10Ph3, log10.(θ0))
lgPs, dlgP_dθs_ad = DiffResults.value(Pres), DiffResults.jacobian(Pres)
#lgPs = log10.(Symboltz.P(pt, ks, θ0...))

p = plot(layout=(2,1), size=(800, 1000), left_margin=bottom_margin=30*Plots.px); display(p)
plot!(p[1], log10.(ks*Symboltz.k0), lgPs; xlabel="lg(k/(h/Mpc))", label="lg(P/(Mpc/h)³)"); display(p)
plot!(p[2], log10.(ks*Symboltz.k0), dlgP_dθs_ad; xlabel = "lg(k/(h/Mpc))", ylabel = "∂ lg(P) / ∂ lg(θᵢ)", labels = "θᵢ=" .* ["Ωr0" "Ωm0" "Ωb0" "h" "As" "Yp"]); display(p)

# TODO: fix finite differences
#dlgP_dθs_fd = FiniteDiff.finite_difference_jacobian(log10Ph3, log10.(θ0); relstep=1e-4) # relstep is important for finite difference accuracy!
#plot!(p[2], log10.(ks*Symboltz.k0), dlgP_dθs_fd; xlabel = "lg(k/(h/Mpc))", ylabel = "∂ lg(P) / ∂ lg(θᵢ)", color = :black); display(p)