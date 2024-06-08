include("../Symboltz.jl")
using .Symboltz
using ModelingToolkit
using Plots; Plots.default(label=nothing)
using ForwardDiff, DiffResults, FiniteDiff

@kwdef struct Parameters
    Ωr0 = 5.5e-5
    Ωc0 = 0.267
    Ωb0 = 0.05
    h = 0.67
    As = 2.1e-9
    Yp = 0.245 # TODO: fix, make work with 0.245
end

par = Parameters()
@named bg = background_ΛCDM()
@named th = thermodynamics_ΛCDM(bg)
@named pt = perturbations_ΛCDM(th, 6)

ls = [2:1:8; 10; 12; 16; 22; 30:15:3000]
θ0 = [par.Ωr0, par.Ωc0, par.Ωb0, par.h, par.As, par.Yp]

# differentiated CMB power spectrum
lgDl(lgθ) = log10.(Dl(pt, ls, (10 .^ lgθ)...))
lgDlres = DiffResults.JacobianResult(Float64.(ls), θ0)
ForwardDiff.jacobian!(lgDlres, lgDl, log10.(θ0))
lgDls, dlgDl_dθs_ad = DiffResults.value(lgDlres), DiffResults.jacobian(lgDlres)
Dls = 10 .^ lgDls

p = plot(layout=(2,1), size=(800, 1000), left_margin=bottom_margin=30*Plots.px); display(p)
plot!(p[1], ls, Dls / 1e-12; xlabel = "l", ylabel = "Dₗ=l (l+1) Cₗ / 2π / 10⁻¹²", title = "CMB power spectrum")
plot!(p[2], ls, dlgDl_dθs_ad; xlabel = "l", ylabel = "∂ lg(Dₗ) / ∂ lg(θᵢ)", labels = "θᵢ=" .* ["Ωr0" "Ωm0" "Ωb0" "H0" "As" "Yp"])