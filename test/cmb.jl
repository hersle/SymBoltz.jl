using ModelingToolkit
using Plots; Plots.default(label=nothing)
using ForwardDiff, DiffResults, FiniteDiff

model = Symboltz.ΛCDM()
par = Symboltz.CosmologicalParameters()
ls = [2:1:8; 10; 12; 16; 22; 30:15:3000]
θ0 = [par.Ωγ0, par.Ων0, par.Ωc0, par.Ωb0, par.h, par.As, par.Yp]

# differentiated CMB power spectrum
lgDl(lgθ) = log10.(Symboltz.Dl(model, ls, (10 .^ lgθ)...))
lgDlres = DiffResults.JacobianResult(Float64.(ls), θ0)
ForwardDiff.jacobian!(lgDlres, lgDl, log10.(θ0))
lgDls, dlgDl_dθs_ad = DiffResults.value(lgDlres), DiffResults.jacobian(lgDlres)
Dls = 10 .^ lgDls

p = plot(layout=(2,1), size=(800, 1000), left_margin=bottom_margin=30*Plots.px); display(p)
plot!(p[1], ls, Dls / 1e-12; xlabel = "l", ylabel = "Dₗ=l (l+1) Cₗ / 2π / 10⁻¹²", title = "CMB power spectrum")
plot!(p[2], ls, dlgDl_dθs_ad; xlabel = "l", ylabel = "∂ lg(Dₗ) / ∂ lg(θᵢ)", labels = "θᵢ=" .* ["Ωγ0" "Ων0" "Ωm0" "Ωb0" "H0" "As" "Yp"])