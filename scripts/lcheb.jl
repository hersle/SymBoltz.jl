using SymBoltz
using Plots

M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)

mode = :TT # :TT or :EE or :ψψ

ls = unique(round.(range(2, 3000, length = 400))) # ℓ ≥ 2 for T/E, ℓ ≥ 10 for ψ # TODO: piecewise Cheb with Limber?
jl = SphericalBesselCache(ls)
Cls = spectrum_cmb(mode, prob, jl) # no interpolation
plot(log.(8e3 .+ ls), Cls .* ls .^ 5; label = nothing)

jl_cub = SphericalBesselCache(CubicSplineInterpolator(ls[begin], ls[end], 60))
jl_cheb = SphericalBesselCache(ChebyshevInterpolator(ls[begin], ls[end], 60; f = l -> log(8e3 + l)))
jl_piecewise = SphericalBesselCache(PiecewiseChebyshevInterpolator((2, 30, 3000), (5, 45)))
Cls_cub = spectrum_cmb(mode, prob, jl_cub, ls)
Cls_cheb = spectrum_cmb(mode, prob, jl_cheb, ls)
Cls_piecewise = spectrum_cmb(mode, prob, jl_piecewise, ls)

p = begin
    plot(layout = grid(2, 1; heights = [0.75, 0.25]), link = :x)
    plot!(ls, @. Cls * ls .* (ls .+ 1); color = :black, linewidth = 2, label = "direct", ylabel = "Dₗ", subplot = 1, xformatter = _ -> "")
    plot!(ls, @. Cls_cub * ls .* (ls .+ 1); color = 1, label = "cubic spline", subplot = 1)
    plot!(ls, @. Cls_cheb * ls .* (ls .+ 1); color = 2, label = "Chebyshev", subplot = 1)
    plot!(ls, @. Cls_piecewise * ls .* (ls .+ 1); color = 3, label = "Piecewise Chebyshev", subplot = 1)
    plot!(ls, @. abs(Cls_cub / Cls - 1); color = 1, yscale = :log10, ylims = (1e-8, 1e0), label = nothing, ylabel = "rel. err.", subplot = 2, xlabel = "l")
    plot!(ls, @. abs(Cls_cheb / Cls - 1); color = 2, label = nothing, subplot = 2)
    plot!(ls, @. abs(Cls_piecewise / Cls - 1); color = 3, label = nothing, subplot = 2)
end
