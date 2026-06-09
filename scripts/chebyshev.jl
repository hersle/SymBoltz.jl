using SymBoltz, Plots, StaticArrays, LinearAlgebra, Printf
Plots.default(dpi = 96, colorbar = nothing, framestyle = :box, grid = false)
M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)

ks = lingrid(1e-2, 2e3; length = 500)
sol = solve(prob, ks)
τ0 = sol[M.τ0]
τrec = sol[M.τrec]

ST = (M.k*M.χ)*M.ST
SE = (M.k*M.χ)^2*M.SE
SM = M.m.Δ
SL = -(M.g.Φ+M.g.Ψ)
pT = surface(sol, M.τ, M.k, ST; xlims = (0.035, 0.085), ylims = (0, 2000), zticks = nothing, zlabel = nothing, title = ST, cam = (40, 25), seriescolor = :plasma, size = (500, 500), margin = -10*Plots.mm)
pE = surface(sol, M.τ, M.k, SE; xlims = (0.035, 0.085), ylims = (0, 2000), zticks = nothing, zlabel = nothing, title = SE, cam = (40, 25), seriescolor = :plasma, size = (500, 500), margin = -10*Plots.mm)
pM = surface(sol, M.τ, M.k, SM; ylims = (0, 2000), zticks = nothing, zlabel = nothing, title = SM, cam = (40, 25), seriescolor = :plasma, size = (500, 500), margin = -10*Plots.mm)
pL = surface(sol, M.τ, log10(M.k), SL; ylims = (0, log10(2e3)), zticks = nothing, zlabel = nothing, title = SL, cam = (40, 25), seriescolor = :plasma, size = (500, 500), margin = -10*Plots.mm)
p = plot(pT, pE, pM, pL, layout = (1, 4), size = (1000, 300), zticks = nothing, margin = -1*Plots.mm, tickfontsize = 7)
savefig(pT, "paper_chebyshev/figures/ST.pdf")
savefig(pE, "paper_chebyshev/figures/SE.pdf")
savefig(pM, "paper_chebyshev/figures/SM.pdf")
savefig(pL, "paper_chebyshev/figures/SL.pdf")
savefig(p, "paper_chebyshev/figures/S.pdf")

# Plot source function examples
#=
function plot_source(S, τs, ks; fk = identity, kwargs...)
    Ss = source_grid(prob, S, τs, ks, sol.bg)
    p = surface(τs, fk.(ks), transpose(Ss); title = "S = $S", xlabel = "τ", ylabel = "k", kwargs...)
    return p
end
ks = lingrid(0.01, 3e3; length = 500)
τs = lingrid(0.04, 0.08; length = 100)
pT = plot_source((M.k*M.χ)*M.ST, τs, ks, cam = (45, 30))
pE = plot_source((M.k*M.χ)^2*M.SE, τs, ks, cam = (45, 30))
#τs = sol.bg.t[begin:end-1]
#pψ = plot_source((M.g.Φ+M.g.Ψ), sol.bg.t, ks, cam = (45, 30), fk = k -> asinh(k/2000))
pm = plot_source(M.m.Δ, sol.bg.t, ks, cam = (45, 30))


S = (M.k*M.χ)^2 * M.ST
#S = (M.k*M.χ)^2 * M.SE
#S = (M.k*M.χ)^2 * M.Sψ
#S = M.m.Δ

τs = sol.bg.t[begin:end-1]
ks_raw = range(1e0, 2e3, length = 500)
Ss_raw = source_grid(prob, S, τs, ks_raw, sol.bg)

kgrid_cheb = ChebyshevInterpolator(ks_raw[begin], ks_raw[end], 32)
Ss_cheb = source_grid(prob, S, τs, ks_raw, kgrid_cheb, sol.bg)

# TODO: interpolate with Cubic splines
kgrid_cubic = CubicSplineInterpolator(lingrid(ks_raw[begin], ks_raw[end]; length = length(kgrid_cheb.ks))) # TODO: try using Chebyshev k-nodes
Ss_cubic = source_grid(prob, S, τs, ks_raw, kgrid_cubic, sol.bg)

xlims = (0.04, 0.07)
p1 = surface(τs, ks_raw, transpose(Ss_raw); xlabel = "τ", ylabel = "k", title = "S(τ, k) = $S, direct integration of $(length(ks_raw))×k", xlims , size = (1000, 600))
p2 = surface(τs, ks_raw, transpose(Ss_cheb); xlabel = "τ", ylabel = "k", title = "S(τ, k) = $S, interpolation from $(order(kgrid_cheb))-order Chebyshev", xlims, size = (1000, 600))
p3 = surface(τs, ks_raw, transpose(Ss_cubic); xlabel = "τ", ylabel = "k", title = "S(τ, k) = $S, interpolation from Cubic spline with $(length(kgrid_cubic.ks))×k", xlims, size = (1000, 600))
#plot(p1, p2, p3, layout = (1, 3), size = (2000, 400))
=#


# TODO: also surface-plot 
function plot_convergence(S, τs, ks_raw; tol = 1e-7, klengths = 25:25:250, kwargs...)
    ptopts = (abstol = tol, reltol = tol)
    Ss_raw = source_grid(prob, S, τs, ks_raw, sol.bg; ptopts)
    norms_cheb = Float64[]
    norms_cubic = Float64[]
    norms_cubic_at_cheb = Float64[]
    error(Ss, p = 2) = norm(Ss .- Ss_raw, p) / length(Ss)^(1/p)
    for klength in klengths
        # 1) Chebyshev interpolation at chebyshev k-nodes
        korder = klength - 1
        kgrid_cheb = ChebyshevInterpolator(ks_raw[begin], ks_raw[end], korder)
        Ss = source_grid(prob, S, τs, ks_raw, kgrid_cheb, sol.bg; ptopts)
        push!(norms_cheb, error(Ss))

        # 2) Cubic spline interpolation at uniform k-nodes
        kgrid_cubic = CubicSplineInterpolator(range(ks_raw[begin], ks_raw[end], klength))
        Ss = source_grid(prob, S, τs, ks_raw, kgrid_cubic, sol.bg; ptopts)
        push!(norms_cubic, error(Ss))

        # 3) Cubic spline interpolation, but at Chebyshev nodes (to rule out effect of k-point selection)
        kgrid = CubicSplineInterpolator(reverse(kgrid_cheb.xs)) # to ascending order
        Ss = source_grid(prob, S, τs, ks_raw, kgrid, sol.bg; ptopts)
        push!(norms_cubic_at_cheb, error(Ss))
    end
    return plot(klengths, [norms_cheb, norms_cubic, norms_cubic_at_cheb]; xticks = klengths, yscale = :log10, xlabel = "number of k-points", ylabel = "‖S(interp) - S(direct)‖", marker = :circle, label = ["Chebyshev interpolation (on Chebyshev grid)" "cubic spline interpolation (on uniform grid)" "cubic spline interpolation (on Chebyshev grid)"], title = "S(τ, k) = $S", kwargs...)
end
pT = plot_convergence(ST, sol.bg.t[begin:end-1], lingrid(1, 2500; length = 500); yticks = 10.0 .^ (-1:4)) # multiply by χ to smooth out χ=0-behavior
pE = plot_convergence(SE, sol.bg.t[begin:end-1], lingrid(1, 2500; length = 500); yticks = 10.0 .^ (-6:0))
pM = plot_convergence(SM, sol.bg.t[begin:end-1], lingrid(1, 2500; length = 500); yticks = 10.0 .^ (-6:0))
pL = plot_convergence(SL, sol.bg.t[begin:end-1], lingrid(1, 2500; length = 500); yticks = 10.0 .^ (-6:0))
savefig(pT, "paper_chebyshev/figures/convT.pdf")
savefig(pE, "paper_chebyshev/figures/convE.pdf")
savefig(pM, "paper_chebyshev/figures/convM.pdf")
savefig(pL, "paper_chebyshev/figures/convL.pdf")

# Plot heatmap of errors for Chebyshev vs cubic splines
function plot_heatmap_errors(S, τs, ks, kgrid1, kgrid2; xvar = τ, tol = 1e-7, kwargs...)
    ptopts = (abstol = tol, reltol = tol)
    Ss = source_grid(prob, S, τs, ks, sol.bg; ptopts)
    Ss1 = source_grid(prob, S, τs, ks, kgrid1, sol.bg; ptopts)
    Ss2 = source_grid(prob, S, τs, ks, kgrid2, sol.bg; ptopts)
    title = "log10(S(interp)-S(direct)): S=$S"
    p1 = heatmap(sol(xvar, τs), ks, transpose(max.(log10.(abs.(Ss1 .- Ss)), -300.0)); title, xlabel = xvar, kwargs...)
    p2 = heatmap(sol(xvar, τs), ks, transpose(max.(log10.(abs.(Ss2 .- Ss)), -300.0)); title, xlabel = xvar, kwargs...)
    return p1, p2
end
kgrid_cheb = ChebyshevInterpolator(1.0, 2000.0, 99)
kgrid_cubic = CubicSplineInterpolator(1.0, 2000.0, 99)
pTcheb, pTcubic = plot_heatmap_errors(ST, sol.bg.t, lingrid(1, 2000; length = 500), kgrid_cheb, kgrid_cubic; ylabel = "k", clims = (-5,  0), cbar = true, xvar = log10(M.g.a), xlims = (-3.2, 0.0))
pEcheb, pEcubic = plot_heatmap_errors(SE, sol.bg.t, lingrid(1, 2000; length = 500), kgrid_cheb, kgrid_cubic; ylabel = "k", clims = (-9, -4), cbar = true, xvar = log10(M.g.a), xlims = (-3.2, 0.0))
pMcheb, pMcubic = plot_heatmap_errors(SM, sol.bg.t, lingrid(1, 2000; length = 500), kgrid_cheb, kgrid_cubic; ylabel = "k", clims = (-5,  0), cbar = true, xvar = log10(M.g.a), xlims = (-3.2, 0.0))
pLcheb, pLcubic = plot_heatmap_errors(SL, sol.bg.t, lingrid(1, 2000; length = 500), kgrid_cheb, kgrid_cubic; ylabel = "k", clims = (-8, -3), cbar = true, xvar = log10(M.g.a), xlims = (-3.2, 0.0))
pTcheb
pTcubic
pEcheb
pEcubic
pMcheb
pMcubic
pLcheb
pLcubic
savefig(pTcheb, "paper_chebyshev/figures/heaterrTcheb.pdf")
savefig(pTcubic, "paper_chebyshev/figures/heaterrTcubic.pdf")
savefig(pEcheb, "paper_chebyshev/figures/heaterrEcheb.pdf")
savefig(pEcubic, "paper_chebyshev/figures/heaterrEcubic.pdf")
savefig(pMcheb, "paper_chebyshev/figures/heaterrMcheb.pdf")
savefig(pMcubic, "paper_chebyshev/figures/heaterrMcubic.pdf")
savefig(pLcheb, "paper_chebyshev/figures/heaterrLcheb.pdf")
savefig(pLcubic, "paper_chebyshev/figures/heaterrLcubic.pdf")

function plot_time_slice(S, τ, ks_raw, klengths; f = identity, tol = 1e-7, kinterpolate = :chebyshev, yerrlims = (1e-10, 1e0), kwargs...)
    ptopts = (abstol = tol, reltol = tol)
    Ss_raw = source_grid(prob, S, [τ], ks_raw, sol.bg; ptopts)
    Ss_raw = Ss_raw[1, :] # slice in time
    margin = 3*Plots.mm
    title = "S(τ = $(round(τ, digits=3)), k) = $S"
    p = plot(ks_raw, Ss_raw; title, label = "$(length(ks_raw)) solved k-modes", color = :black, linewidth = 4, size = (600, 600), layout = grid(2, 1, heights=[0.8, 0.2]), link = :x, margin)
    for (color, klength) in enumerate(klengths)
        if kinterpolate == :chebyshev
            kgrid = ChebyshevInterpolator(ks_raw[begin], ks_raw[end], klength - 1; f)
            label = "$klength Chebyshev k-points"
        elseif kinterpolate == :cubic
            kgrid = CubicSplineInterpolator(lingrid(ks_raw[begin], ks_raw[end]; length = klength); f)
            label = "$klength cubic spline k-points"
        elseif kinterpolate == :equispaced
            @assert f == identity
            kgrid = EquispacedInterpolator(ks_raw[begin], ks_raw[end], klength - 1)
            label = "$klength equispaced k-points"
        end
        Ss = source_grid(prob, S, [τ], ks_raw, kgrid, sol.bg; ptopts)
        Ss = Ss[1, :] # slice in time
        plot!(p, ks_raw, Ss; color, label, xformatter = _ -> "", bottom_margin = -5*Plots.mm, subplot = 1, kwargs...)
        plot!(p, ks_raw, abs.(Ss .- Ss_raw); color, xlabel = M.k, ylabel = "error", yscale = :log10, label = nothing, ylims = yerrlims, yticks = 10.0 .^ (-8:8), subplot = 2)
    end
    return p
end
pcheb = plot_time_slice(1e-5*ST, τrec, lingrid(1, 2000; length = 500), 20:20:100; kinterpolate = :chebyshev, yerrlims = (1e-9, 1e0), ylims = (-1.3, 1.3), legend_position = :topright) # at recombination (peak of visibility function)
pcubic = plot_time_slice(1e-5*ST, τrec, lingrid(1, 2000; length = 500), [20:20:100; 200:100:400]; kinterpolate = :cubic, yerrlims = (1e-9, 1e0), ylims = (-1.3, 1.3), legend_position = :topright)
pequi = plot_time_slice(1e-5*ST, τrec, lingrid(1, 2000; length = 500), 20:20:60; kinterpolate = :equispaced, yerrlims = (1e-9, 1e0), ylims = (-1.3, 1.3), legend_position = :topright)
#pchebzoom = plot_time_slice(ST, τrec, lingrid(500.0, 750.0; length = 200), 5:5:20; kinterpolate = :chebyshev, yerrlims = (1e-4, 1e5), ylims = (-1.3e5, 1.3e5))
#pcubiczoom = plot_time_slice(ST, τrec, lingrid(500.0, 750.0; length = 200), 5:5:80; kinterpolate = :cubic, yerrlims = (1e-4, 1e5), ylims = (-1.3e5, 1.3e5))
savefig(pcheb, "paper_chebyshev/figures/reccheb.pdf")
savefig(pcubic, "paper_chebyshev/figures/reccubic.pdf")
savefig(pequi, "paper_chebyshev/figures/recequi.pdf")
#savefig(pchebzoom, "paper_chebyshev/figures/recchebzoom.pdf")
#savefig(pcubiczoom, "paper_chebyshev/figures/reccubiczoom.pdf")

#pcheb_zoom = plot_time_slice(M.k*M.χ*M.ST, τrec, lingrid(500.0, 750.0; length = 200), (4:4:20) .+ 1; kinterpolate = :chebyshev, yerrlims = (1e-4, 1e5), ylims = (-1.3e5, 1.3e5))
#pcubic_zoom = plot_time_slice(M.k*M.χ*M.ST, τrec, lingrid(500.0, 750.0; length = 200), (2 .^ (2:8)) .+ 1; kinterpolate = :cubic, yerrlims = (1e-4, 1e5), ylims = (-1.3e5, 1.3e5))
# note similar convergence rate for Chebyshev linear-n vs. cubic spline exponential-n

# Plot decay of Chebyshev coefficients
function plot_chebyshev_coefficients(S, τs, kgrid; tols = 1e-7, kwargs...)
    p = plot(; yscale = :log10, xlabel = "n", ylabel = "|cₙ|", title = "S = $S")
    for tol in tols
        interps = source_grid_interp(prob, S, τs, kgrid, sol.bg; ptopts = (abstol = tol, reltol = tol))
        ns = 0:order(kgrid)
        for (i, interp) in enumerate(interps)
            τ = τs[i]
            cs = abs.(interp.coefs)
            plot!(p, ns, cs; label = @sprintf("τ = %.3f, pert. ODE tol. = %.0e", τ, tol), kwargs...)
        end
    end
    return p
end
kgrid = ChebyshevInterpolator(1.0, 2000.0, 250)
pT = plot_chebyshev_coefficients(ST, [τrec], kgrid; tols = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7], yticks = 10.0 .^ (-6:5), ylims = (1e-6, 1e5))
pM = plot_chebyshev_coefficients(SM, [τ0], kgrid; tols = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7], yticks = 10.0 .^ (-6:5), ylims = (1e-6, 1e5))
savefig(pT, "paper_chebyshev/figures/coeffsT.pdf")
savefig(pM, "paper_chebyshev/figures/coeffsM.pdf")

#p = plot_chebyshev_coefficients(M.k*M.Sψ, [0.9*τ0], kgrid; ylims = (1e0, 1e8), label = nothing)
#savefig(p, "paper_chebyshev/figures/coeffs.pdf")

# Plot interpolation of matter power spectrum; using chebyshev; but for linear vs log k-grid
#=
S = -M.m.Δ
#logks_raw = range(-1, 5, length = 200)
logks_raw = range(2, 4, length = 200)
ks_raw = 10 .^ (logks_raw)
τ = τs[end]
Ss_raw = source_grid(prob, S, [τ], ks_raw, sol.bg; ptopts)
Ss_raw = Ss_raw[1, :] # time slice
margin = 3*Plots.mm
title = "S(τ = $(round(τ, digits=3)), k) = $S"
xlims = extrema(ks_raw)
pmlin = plot(ks_raw, Ss_raw; xscale = :log10, yscale = :log10, title, label = "$(length(ks_raw)) solved k-modes", color = :black, xlims, linewidth = 4, size = (600, 500), layout = grid(2, 1, heights=[0.8, 0.2]), legend_position = :bottomright, link = :x, margin)
for (color, klength) in enumerate([16,32,64,128])
    # 1) Chebyshev interpolation at chebyshev k-nodes
    kgrid = ChebyshevInterpolator(ks_raw[begin], ks_raw[end], klength - 1; f = log, f⁻¹ = exp)
    Ss = source_grid(prob, log(S), [τ], ks_raw, kgrid, sol.bg; ptopts)
    Ss = exp.(stack(Ss)[1, :]) # time slice
    plot!(pmlin, ks_raw, Ss; color, ylabel = S, label = "$klength Chebyshev points", xformatter = _ -> "", bottom_margin = -5*Plots.mm, subplot = 1)
    plot!(pmlin, ks_raw, abs.(Ss./Ss_raw.-1); yscale = :log10, color, xlabel = M.k, ylabel = "error", ylims = (1e-6, 1e-1), label = nothing, subplot = 2)
end
pmlin
# TODO: now plot the same, but using a linear grid for k
=#

# Plot Chebyshev points
function plot_chebyshev_points(; maxorder = 50, kwargs...)
    p = plot(xlabel = "Chebyshev points", ylabel = "polynomial order", grid = false, xlims = (-1, +1), xticks = -1:1, ylims = (0, maxorder), margin = 2*Plots.mm, kwargs...)
    for n in 1:maxorder
        ks = chebgrid(-1, +1; order = n)
        scatter!(p, ks, fill(n, length(ks)); label = nothing, color = :black, markersize = 2.5)
    end
    return p
end
p = plot_chebyshev_points(; maxorder = 50)
savefig(p, "paper_chebyshev/figures/chebpoints.pdf")

n = 0:6
p_poly = plot(x, cos.(n'.*acos.(x)); xticks = -1:1, yticks = -1:1, xformatter = _ -> "", bottom_margin = -3*Plots.mm, label = "T" .* string.(Char.(0x2080 .+ n')) .* "(x)", legend = :outertop, legendcolumns = -1)
p_nodes = plot(; xticks = -1:1, xlabel = "x", ylabel = "order", xlims = (-1, 1), ylims = (0, 51), grid = false)
for order in 1:50
    xnodes = chebgrid(-1, 1; order)
    scatter!(p_nodes, xnodes, fill(order, length(xnodes)); label = nothing, color = :black, markersize = 2.5)
end
p = plot(p_poly, p_nodes; layout = grid(2, 1, heights = [2/3, 1/3]), link = :x, size = (600, 600))

savefig(p, "paper_chebyshev/figures/chebpoly.pdf")

# Plot Chebyshev points

# Compare CMB spectra with Chebyshev vs cubic spline interpolation
ls = range(log(2), log(2500), length=200) .|> exp .|> round |> unique .|> Int
jl = SphericalBesselCache(ls)
mode = :TT
xs = 0.0:0.0002:1.0
Δkτ0 = 2π/16 # fixes low-l jagginess
kmin, kmax = 0.01, 1000.0 # kτ0max = 2*lmax for TT/EE, 10*lmax for ψψ
normalization = :Dl
kgrid_ref = ChebyshevInterpolator(kmin, kmax, 2^10)
klengths = [2^n+1 for n in 4:8]
kgrids_cubic = [CubicSplineInterpolator(lingrid(kmin, kmax; length)) for length in klengths]
kgrids_cheb = [ChebyshevInterpolator(kmin, kmax, length-1) for length in klengths]
Dls_ref = spectrum_cmb(mode, prob, jl; normalization, kgrid = kgrid_ref, Δkτ0, xs)
Dlss_cheb = [spectrum_cmb(mode, prob, jl; normalization, kgrid = kgrid, Δkτ0, xs) for kgrid in kgrids_cheb]
Dlss_cubic = [spectrum_cmb(mode, prob, jl; normalization, kgrid, Δkτ0, xs) for kgrid in kgrids_cubic]
Dlss = Dlss_cheb
p = begin
    labelfunc(klength) = "$klength " * (Dlss === Dlss_cheb ? "Chebyshev" : "cubic spline") * " k-points"
    xscale = :log10
    legend_position = :topleft
    plot(ls, Dls_ref; ylims = (0, 1e-9), title = "mode = $mode", xscale, color = :black, linewidth = 3, label = labelfunc(length(kgrid_ref.ks)), layout = grid(2, 1, heights=[0.75, 0.25]), size = (600, 500), legend_position, subplot = 1)
    plot!(ls, Dlss; xscale, color = permutedims(eachindex(klengths)), ylabel = "D(ℓ) = ℓ(ℓ+1) C(ℓ) / 2π", xformatter = _ -> "", label = labelfunc.(permutedims(klengths)), bottom_margin = -3*Plots.mm, subplot = 1)
    plot!(ls, [abs.(Dls./Dls_ref.-1) for Dls in Dlss]; xscale, xlabel = "ℓ", ylabel = "rel. error", color = permutedims(eachindex(klengths)), yscale = :log10, label = nothing, ylims = (1e-7, 1e-1), subplot = 2)
end
fname = "spectrum" * Dict(:TT => "T", :EE => "E", :ψψ => "psi")[mode] * "_" * (Dlss === Dlss_cheb ? "cheb" : "cubic") * ".pdf"
savefig(p, "paper_chebyshev/figures/$fname")

# TODO: test on arbitrary normal perturbation variables
S = M.g.Φ
τs = sol.bg.t
ks = range(0.1, 3000.0, length = 100)
Ss = source_grid(prob, S, τs, ks, sol.bg)
Ss_cheb = source_grid_chebyshev(prob, S, τs, ks, 20, sol.bg)


# TODO: plot ST at recombination with domain transformation

k0 = 2000.0 # want linear sampling for k < k0, logarithmic for k > k0
f = k -> tanh(k/k0) # tanh(1) ≈ 0.76, tanh(2) ≈ 0.96
f⁻¹ = f -> k0 * atanh(f)
#f = k -> asinh(k/k0)
#γ = 1.0
#f = k -> asinh((k/k0)^γ)
#f⁻¹ = f -> k0 * (sinh(f))^(1/γ)
#f = k -> log(1 + (k/k0)^γ)
#f⁻¹ = f -> k0 * (exp(f) - 1)^(1/γ)
#f = k -> SymBoltz.smoothifelse(k/k0-1, k, k+log(k/k0); k=1.0)
#f = k -> ifelse(k<k0, k, k0+0.01(k-k0)) # not smooth
plot(f, xlims = (0, 1e4))
kinterp = ChebyshevInterpolator(1e-2, 2e4, 250; f)
S = ST
τ = τrec
Sgrid = source_grid(prob, S, [τ], kinterp.xs)[1, :]
ks = lingrid(minimum(kinterp), maximum(kinterp); length = 5000)
Ss = source_grid(prob, S, [τ], ks, kinterp)[1, :]
p = begin
    color = :black
    xformatter = _ -> ""
    yformatter = _ -> ""
    p = plot(layout = (3, 1), size = (600, 1000), margin = 5*Plots.mm)
    plot!(p, kinterp.xs, kinterp.ys; color, marker = :circle, xlabel = "k", ylabel = "f(k)", label = nothing, subplot = 1)
    vline!(p, [k0], subplot = 1, color = :gray, linestyle = :dash, label = "k = k₀")
    hline!(p, [f(k0)], subplot = 1, color = :gray, linestyle = :dash, label = nothing)
    plot!(p, ks, Ss; color, xlabel = "k", ylabel = "S", xticks = kinterp.xs, label = nothing, xformatter, yformatter, subplot = 2)
    scatter!(p, kinterp.xs, Sgrid; color, xformatter, yformatter, label = nothing, subplot = 2)
    vline!(p, [k0], subplot = 2, color = :gray, linestyle = :dash, label = nothing)
    plot!(p, f.(ks), Ss; color, xlabel = "f(k)", ylabel = "S", label = nothing, xticks = kinterp.ys, xformatter, yformatter, subplot = 3)
    scatter!(p, kinterp.ys, Sgrid; color, xformatter, yformatter, label = nothing, subplot = 3)
    vline!(p, [f(k0)], subplot = 3, color = :gray, linestyle = :dash, label = nothing)
end

# Plot convergence (cₙ as a function of n) for fixed N, but different f(k)
S = ST
p = plot(; layout = (3, 1), size = (600, 1000), left_margin = 5*Plots.mm, right_margin = 5*Plots.mm)
for (f, fname) in [
    (identity, "k"),
    (k->log(k/2000), "log(k/2000)"),
    (k->asinh(k/2000), "asinh(k/2000)"),
    (k->tanh(k/2000), "tanh(k/2000)"),
]
    kinterp = ChebyshevInterpolator(1e-2, 1e4, 250; f)
    interp = only(source_grid_interp(prob, S, [τrec], kinterp, sol.bg))
    ns = 0:order(kinterp)
    cs = abs.(interp.coefs)
    plot!(p[1, 1], ns, cs; xlabel = "n", ylabel = "|cₙ|", yscale = :log10, yticks = 10.0 .^ (-5:5), title = "S(τ = $(round(τrec; digits=3)), k) = $S", label = nothing)

    fmin, fmax = extrema(kinterp.ys)
    plot!(p[2, 1], kinterp.xs, (kinterp.ys .- fmin) ./ (fmax - fmin); xlabel = "k", ylabel = "f(k)", label = "f(k) = $fname", xticks = [kinterp.xs[begin], 2000, kinterp.xs[end]], yticks = nothing)

    fs = range(fmin, fmax; length = 1000)
    Ss = interp.(fs)
    plot!(p[3, 1], (fs .- fmin) ./ (fmax - fmin), Ss; xlabel = "f(k)", ylabel = S, label = nothing, xticks = nothing, yticks = nothing)
end
p
