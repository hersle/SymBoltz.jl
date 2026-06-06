using SymBoltz, Plots, StaticArrays, LinearAlgebra
Plots.default(dpi = 96, colorbar = nothing, framestyle = :box)
M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
sol = solve(prob)

# Plot source function examples
ks = range(0.1, 3000.0, length = 500)
sol = solve(prob, ks)
camera = (30, 30)
ylims = (0, 3000)
yticks = 0:500:3000
size = (550, 450)
margin = -20*Plots.mm
pT = surface(sol, M.τ, M.k, (M.k*M.χ)^2*M.ST; xlims = (0.035, 0.085), xticks = 0.04:0.01:0.08, ylims, yticks, camera, size, margin)
savefig(pT, "paper_chebyshev/figures/ST.pdf")
pE = surface(sol, M.τ, M.k, (M.k*M.χ)^2*M.SE; xlims = (0.035, 0.085), xticks = 0.04:0.01:0.08, ylims, yticks, camera, size, margin)
savefig(pE, "paper_chebyshev/figures/SE.pdf")
pψ = surface(sol, M.τ, M.k, -(M.k*M.χ)^2*M.Sψ; ylims, yticks, camera, size, margin)
savefig(pψ, "paper_chebyshev/figures/Spsi.pdf")
pm = surface(sol, M.τ, M.k, M.m.Δ; ylims, yticks, camera, size, margin)
savefig(pm, "paper_chebyshev/figures/Sm.pdf")

#S = (M.k*M.χ)^2 * M.ST
#S = (M.k*M.χ)^2 * M.SE
#S = (M.k*M.χ)^2 * M.Sψ
S = M.m.Δ

τs = sol.bg.t[begin:end-1]
ks_raw = range(1e0, 2e3, length = 1000)
Ss_raw = source_grid(prob, SVector(S), τs, ks_raw, sol.bg)
Ss_raw = stack(Ss_raw)[1, :, :]

korder = 32
ks_cheb = SymBoltz.chebpoints(korder, ks_raw[begin], ks_raw[end])
Ss_cheb, _, _ = source_grid_chebyshev(prob, SVector(S), τs, ks_raw, sol.bg; order = korder)
Ss_cheb = stack(Ss_cheb)[1, :, :]
abserr_cheb = abs.(Ss_cheb .- Ss_raw)
relerr_cheb = abserr_cheb ./ abs.(Ss_raw)

# TODO: interpolate with Cubic splines
#klength = 32
#ks_cubic = range(ks_raw[begin], ks_raw[end], klength) # TODO: try using Chebyshev k-nodes
ks_cubic = reverse(ks_cheb)
Ss_cubic = source_grid(prob, SVector(S), τs, ks_cubic, sol.bg)
Ss_cubic = stack(Ss_cubic)[1, :, :]
Ss_cubic = source_grid(Ss_cubic, ks_cubic, ks_raw) # interpolate with Cubic splines
abserr_cubic = abs.(Ss_cubic .- Ss_raw)
relerr_cubic = abserr_cubic ./ abs.(Ss_raw)

xlims = (0.04, 0.07)
p1 = surface(τs, ks_raw, transpose(Ss_raw); xlabel = "τ", ylabel = "k", title = "S(τ, k) = $S, direct integration of $(length(ks_raw))×k", xlims , size = (1000, 600))
p2 = surface(τs, ks_raw, transpose(Ss_cheb); xlabel = "τ", ylabel = "k", title = "S(τ, k) = $S, interpolation from $(korder)-order Chebyshev", xlims, size = (1000, 600))
p3 = surface(τs, ks_raw, transpose(Ss_cubic); xlabel = "τ", ylabel = "k", title = "S(τ, k) = $S, interpolation from Cubic spline with $(length(ks_cubic))×k", xlims, size = (1000, 600))
#plot(p1, p2, p3, layout = (1, 3), size = (2000, 400))



# Plot convergence
klengths = 20:20:500
norms_cheb = Float64[]
norms_cubic = Float64[]
norms_cubic_at_cheb = Float64[]
error(Ss, p = 2) = norm(Ss .- Ss_raw, p) / length(Ss)^(1/p)
for klength in klengths
    # 1) Chebyshev interpolation at chebyshev k-nodes
    korder = klength - 1
    Ss, _, _ = source_grid_chebyshev(prob, SVector(S), τs, ks_raw, sol.bg; order = korder)
    Ss = stack(Ss)[1, :, :]
    push!(norms_cheb, error(Ss))

    # 2) Cubic spline interpolation at uniform k-nodes
    ks = range(ks_raw[begin], ks_raw[end], klength)
    Ss = source_grid(prob, SVector(S), τs, ks, sol.bg)
    Ss = stack(Ss)[1, :, :]
    Ss = source_grid(Ss, ks, ks_raw) # interpolate with Cubic splines
    push!(norms_cubic, error(Ss))

    # 3) Cubic spline interpolation, but at Chebyshev nodes (to rule out effect of k-point selection)
    ks = reverse(SymBoltz.chebpoints(korder, ks_raw[begin], ks_raw[end])) # to ascending order
    Ss = source_grid(prob, SVector(S), τs, ks, sol.bg)
    Ss = stack(Ss)[1, :, :]
    Ss = source_grid(Ss, ks, ks_raw) # interpolate with Cubic splines
    push!(norms_cubic_at_cheb, error(Ss))
end
pconv = plot(klengths, [norms_cheb, norms_cubic, norms_cubic_at_cheb]; yscale = :log10, xticks=20:40:500, xlabel = "number of k-points", ylabel = "‖S(interp) - S(direct)‖", marker = :circle, label = ["Chebyshev interpolation (on Chebyshev grid)" "cubic spline interpolation (on uniform grid)" "cubic spline interpolation (on Chebyshev grid)"], title = "S(τ, k) = $S")
savefig(pconv, "paper_chebyshev/figures/convm.pdf")

# Plot ST at recombination for different interpolation lengths
ks_raw = range(1.0, 2500.0, length = 500)
iτ = argmax(sol[M.b.v])
τ = τs[iτ]
Ss_raw = source_grid(prob, SVector(S), [τ], ks_raw, sol.bg)
Ss_raw = stack(Ss_raw)[1, 1, :]
margin = 3*Plots.mm
title = "S(τ = $(round(τ, digits=3)), k) = $S"
xlims = (0, ks_raw[end])
pinterp_cheb = plot(ks_raw, Ss_raw; title, label = "$(length(ks_raw)) solved k-modes", color = :black, xlims, linewidth = 4, size = (600, 600), layout = grid(2, 1, heights=[0.8, 0.2]), link = :x, margin)
for (color, klength) in enumerate(10:10:70)
    # 1) Chebyshev interpolation at chebyshev k-nodes
    korder = klength - 1
    Ss, _, _ = source_grid_chebyshev(prob, SVector(S), [τ], ks_raw, sol.bg; order = korder)
    Ss = stack(Ss)[1, 1, :]
    plot!(pinterp_cheb, ks_raw, Ss; color, ylabel = S, label = "$klength Chebyshev points", xformatter = _ -> "", bottom_margin = -5*Plots.mm, subplot = 1)
    plot!(pinterp_cheb, ks_raw, log10.(abs.(Ss .- Ss_raw)); color, xlabel = M.k, ylabel = "log10(|error|)", ylims = (0, 10), label = nothing, subplot = 2)
end
pinterp_cheb
savefig(pinterp_cheb, "paper_chebyshev/figures/interpT_cheb.pdf")

pinterp_cubic = plot(ks_raw, Ss_raw; title, label = "$(length(ks_raw)) solved k-modes", color = :black, xlims, linewidth = 4, size = (600, 600), layout = grid(2, 1, heights=[0.8, 0.2]), link = :x, margin)
for (color, klength) in enumerate(10:10:70)
    # 2) Cubic spline interpolation at uniform k-nodes
    ks = range(ks_raw[begin], ks_raw[end], klength)
    #ks = reverse(SymBoltz.chebpoints(klength-1, ks_raw[begin], ks_raw[end])) # at Chebyshev nodes
    Ss = source_grid(prob, SVector(S), [τ], ks, sol.bg)
    Ss = stack(Ss)[1, :, :]
    Ss = source_grid(Ss, ks, ks_raw) # interpolate with Cubic splines
    Ss = Ss[1, :]
    plot!(pinterp_cubic, ks_raw, Ss; color, ylabel = S, label = "$klength cubic spline points", xformatter = _ -> "", bottom_margin = -5*Plots.mm, subplot = 1)
    plot!(pinterp_cubic, ks_raw, log10.(abs.(Ss .- Ss_raw)); xlabel = M.k, color, ylabel = "log10(error)", ylims = (0, 10), label = nothing, subplot = 2)
end
pinterp_cubic
savefig(pinterp_cubic, "paper_chebyshev/figures/interpT_cubic.pdf")

# Plot interpolation of matter power spectrum; using chebyshev; but for linear vs log k-grid
S = -M.m.Δ
#logks_raw = range(-1, 5, length = 200)
logks_raw = range(2, 4, length = 200)
ks_raw = 10 .^ (logks_raw)
τ = τs[end]
Ss_raw = source_grid(prob, SVector(S), [τ], ks_raw, sol.bg)
Ss_raw = stack(Ss_raw)[1, 1, :]
margin = 3*Plots.mm
title = "S(τ = $(round(τ, digits=3)), k) = $S"
xlims = extrema(ks_raw)
pmlin = plot(ks_raw, Ss_raw; xscale = :log10, yscale = :log10, title, label = "$(length(ks_raw)) solved k-modes", color = :black, xlims, linewidth = 4, size = (600, 500), layout = grid(2, 1, heights=[0.8, 0.2]), legend_position = :bottomright, link = :x, margin)
for (color, klength) in enumerate([16,32,64,128])
    # 1) Chebyshev interpolation at chebyshev k-nodes
    korder = klength - 1
    Ss, _, _ = source_grid_chebyshev(prob, SVector(log(S)), [τ], ks_raw, sol.bg; order = korder, f = log, f⁻¹ = exp)
    Ss = exp.(stack(Ss)[1, 1, :])
    plot!(pmlin, ks_raw, Ss; color, ylabel = S, label = "$klength Chebyshev points", xformatter = _ -> "", bottom_margin = -5*Plots.mm, subplot = 1)
    plot!(pmlin, ks_raw, abs.(Ss./Ss_raw.-1); yscale = :log10, color, xlabel = M.k, ylabel = "error", ylims = (1e-6, 1e-1), label = nothing, subplot = 2)
end
pmlin
# TODO: now plot the same, but using a linear grid for k

# Plot Chebyshev points
ppoints = plot(xlabel = "Chebyshev points", ylabel = "number of points", xticks = ([-1, +1], ["kₘᵢₙ", "kₘₐₓ"]), grid = false, xlims = (-1, +1), ylims = (0, 80), size = (600, 300), margin = 2*Plots.mm)
for klength in 2:80
    korder = klength - 1
    ks = SymBoltz.chebpoints(korder, -1, +1)
    scatter!(ppoints, ks, fill(length(ks), klength); label = nothing, color = :black, markersize = 2.5)
end
ppoints
savefig(ppoints, "paper_chebyshev/figures/points.pdf")

# Compare CMB spectra with Chebyshev vs cubic spline interpolation
ls = range(log(2), log(2500), length=200) .|> exp .|> round |> unique .|> Int
jl = SphericalBesselCache(ls)
mode = :TT
kτ0s = 0.1*jl.l[begin]:2π/2:2*jl.l[end] # 2*lmax for TT/EE, 10*lmax for ψψ
normalization = :Dl
kinterpolate = :adaptive
klength0 = 2^10+1 # reference
klengths = [2^n+1 for n in 4:9] # 30:10:80 for TT/EE; 100:100:800 for ψψ
Dls0 = spectrum_cmb(mode, prob, jl; normalization, kτ0s, kinterpolate = :chebyshev, coarse_length = klength0)
Dlss = [spectrum_cmb(mode, prob, jl; normalization, kτ0s, kinterpolate, coarse_length = klength, sourceopts = (refine = false,)) for klength in klengths] # Chebyshev
p = begin
    labelfunc(klength) = "$klength " * (kinterpolate == :chebyshev ? "Chebyshev" : "cubic spline") * " k-points"
    xscale = :log10
    legend_position = :topleft
    plot(ls, Dls0; ylims = (0, 1e-9), title = "mode = $mode", xscale, color = :black, linewidth = 3, label = labelfunc(klength0), layout = grid(2, 1, heights=[0.75, 0.25]), size = (600, 500), legend_position, subplot = 1)
    plot!(ls, Dlss; xscale, color = permutedims(eachindex(klengths)), ylabel = "D(ℓ) = ℓ(ℓ+1) C(ℓ) / 2π", xformatter = _ -> "", label = labelfunc.(permutedims(klengths)), bottom_margin = -3*Plots.mm, subplot = 1)
    plot!(ls, [abs.(Dls./Dls0.-1) for Dls in Dlss]; xscale, xlabel = "ℓ", ylabel = "error", color = permutedims(eachindex(klengths)), yscale = :log10, label = nothing, ylims = (1e-7, 1e-1), subplot = 2)
end
fname = "spectrum" * Dict(:TT => "T", :EE => "E", :ψψ => "psi")[mode] * "_" * (kinterpolate == :chebyshev ? "cheb" : "cubic") * ".pdf"
savefig(p, "paper_chebyshev/figures/$fname")
