using SymBoltz, Plots, StaticArrays, LinearAlgebra, Printf
Plots.default(dpi = 96, colorbar = nothing, framestyle = :box, grid = false)
M = ΛCDM()
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)

sol = solve(prob; bgopts)
τ0 = sol[M.τ0]
τrec = sol[M.τrec]

# Plot source function examples
function plot_source(S, τs, ks; kwargs...)
    Ss = source_grid(prob, S, τs, ks, sol.bg)
    p = surface(τs, ks, transpose(Ss); title = "S = $S", xlabel = "τ", ylabel = "k", kwargs...)
    return p
end
ks = lingrid(0.01, 3e3; length = 500)
τs = lingrid(0.04, 0.08; length = 100)
pT = plot_source(M.k*M.ST, τs, ks, cam = (45, 30))
pE = plot_source(M.k^2*M.SE, τs, ks, cam = (45, 30))
pψ = plot_source(M.k*M.Sψ, τs[begin:end-1], ks, cam = (45, 30))
pm = plot_source(M.m.Δ, sol.bg.t, ks, cam = (45, 30))

S = (M.k*M.χ)^2 * M.ST
#S = (M.k*M.χ)^2 * M.SE
#S = (M.k*M.χ)^2 * M.Sψ
#S = M.m.Δ

τs = sol.bg.t[begin:end-1]
ks_raw = range(1e0, 2e3, length = 500)
Ss_raw = source_grid(prob, S, τs, ks_raw, sol.bg)

kgrid_cheb = ChebyshevWavenumberGrid(ks_raw[begin], ks_raw[end], 32)
Ss_cheb = source_grid(prob, S, τs, ks_raw, kgrid_cheb, sol.bg)

# TODO: interpolate with Cubic splines
kgrid_cubic = CubicSplineWavenumberGrid(lingrid(ks_raw[begin], ks_raw[end]; length = length(kgrid_cheb.ks))) # TODO: try using Chebyshev k-nodes
Ss_cubic = source_grid(prob, S, τs, ks_raw, kgrid_cubic, sol.bg)

xlims = (0.04, 0.07)
p1 = surface(τs, ks_raw, transpose(Ss_raw); xlabel = "τ", ylabel = "k", title = "S(τ, k) = $S, direct integration of $(length(ks_raw))×k", xlims , size = (1000, 600))
p2 = surface(τs, ks_raw, transpose(Ss_cheb); xlabel = "τ", ylabel = "k", title = "S(τ, k) = $S, interpolation from $(order(kgrid_cheb))-order Chebyshev", xlims, size = (1000, 600))
p3 = surface(τs, ks_raw, transpose(Ss_cubic); xlabel = "τ", ylabel = "k", title = "S(τ, k) = $S, interpolation from Cubic spline with $(length(kgrid_cubic.ks))×k", xlims, size = (1000, 600))
#plot(p1, p2, p3, layout = (1, 3), size = (2000, 400))



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
        kgrid_cheb = ChebyshevWavenumberGrid(ks_raw[begin], ks_raw[end], korder)
        Ss = source_grid(prob, S, τs, ks_raw, kgrid_cheb, sol.bg; ptopts)
        push!(norms_cheb, error(Ss))

        # 2) Cubic spline interpolation at uniform k-nodes
        kgrid_cubic = CubicSplineWavenumberGrid(range(ks_raw[begin], ks_raw[end], klength))
        Ss = source_grid(prob, S, τs, ks_raw, kgrid_cubic, sol.bg; ptopts)
        push!(norms_cubic, error(Ss))

        # 3) Cubic spline interpolation, but at Chebyshev nodes (to rule out effect of k-point selection)
        kgrid = CubicSplineWavenumberGrid(reverse(kgrid_cheb.ks)) # to ascending order
        Ss = source_grid(prob, S, τs, ks_raw, kgrid, sol.bg; ptopts)
        push!(norms_cubic_at_cheb, error(Ss))
    end
    return plot(klengths, [norms_cheb, norms_cubic, norms_cubic_at_cheb]; xticks = klengths, yscale = :log10, xlabel = "number of k-points", ylabel = "‖S(interp) - S(direct)‖", marker = :circle, label = ["Chebyshev interpolation (on Chebyshev grid)" "cubic spline interpolation (on uniform grid)" "cubic spline interpolation (on Chebyshev grid)"], title = "S(τ, k) = $S", kwargs...)
end
pT = plot_convergence(M.k*M.χ*M.ST, sol.bg.t[begin:end-1], lingrid(1, 2500; length = 500); yticks = 10.0 .^ (-1:4)) # multiply by χ to smooth out χ=0-behavior
pE = plot_convergence((M.k*M.χ)^2*M.SE, sol.bg.t[begin:end-1], lingrid(1, 2500; length = 500); yticks = 10.0 .^ (-6:0))

function plot_time_slice(S, τ, ks_raw, klengths; f = identity, tol = 1e-7, kinterpolate = :chebyshev, yerrlims = (1e-10, 1e0), kwargs...)
    ptopts = (abstol = tol, reltol = tol)
    Ss_raw = source_grid(prob, S, [τ], ks_raw, sol.bg; ptopts)
    Ss_raw = Ss_raw[1, :] # slice in time
    margin = 3*Plots.mm
    title = "S(τ = $(round(τ, digits=3)), k) = $S"
    p = plot(ks_raw, Ss_raw; title, label = "$(length(ks_raw)) solved k-modes", color = :black, linewidth = 4, size = (600, 600), layout = grid(2, 1, heights=[0.8, 0.2]), link = :x, margin)
    for (color, klength) in enumerate(klengths)
        if kinterpolate == :chebyshev
            kgrid = ChebyshevWavenumberGrid(ks_raw[begin], ks_raw[end], klength - 1; f)
            Ss = source_grid(prob, S, [τ], ks_raw, kgrid, sol.bg; ptopts)
            label = "$klength Chebyshev k-points"
        elseif kinterpolate == :cubic
            kgrid = CubicSplineWavenumberGrid(lingrid(ks_raw[begin], ks_raw[end]; length = klength); f)
            Ss = source_grid(prob, S, [τ], ks_raw, kgrid, sol.bg; ptopts)
            label = "$klength cubic spline k-points"
        end
        Ss = Ss[1, :] # slice in time
        plot!(p, ks_raw, Ss; color, ylabel = S, label, xformatter = _ -> "", bottom_margin = -5*Plots.mm, subplot = 1, kwargs...)
        plot!(p, ks_raw, abs.(Ss .- Ss_raw); color, xlabel = M.k, ylabel = "error", yscale = :log10, label = nothing, ylims = yerrlims, yticks = 10.0 .^ (-8:8), subplot = 2)
    end
    return p
end
pcheb = plot_time_slice(M.k*M.χ*M.ST, τrec, lingrid(1, 2500; length = 500), 20:20:80; kinterpolate = :chebyshev, yerrlims = (1e-3, 1e5), ylims = (-1.3e5, 1.3e5)) # at recombination (peak of visibility function)
pcubic = plot_time_slice(M.k*M.χ*M.ST, τrec, lingrid(1, 2500; length = 500), 20:20:80; kinterpolate = :cubic, yerrlims = (1e-3, 1e5), ylims = (-1.3e5, 1.3e5))
pcheb_zoom = plot_time_slice(M.k*M.χ*M.ST, τrec, lingrid(500.0, 750.0; length = 200), (4:4:20) .+ 1; kinterpolate = :chebyshev, yerrlims = (1e-4, 1e5))
pcubic_zoom = plot_time_slice(M.k*M.χ*M.ST, τrec, lingrid(500.0, 750.0; length = 200), (2 .^ (2:8)) .+ 1; kinterpolate = :cubic, yerrlims = (1e-4, 1e5))
# note similar convergence rate for Chebyshev linear-n vs. cubic spline exponential-n

# Plot decay of Chebyshev coefficients
function plot_chebyshev_coefficients(S, τs, kgrid; tols = 1e-5, kwargs...)
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
kgrid = ChebyshevWavenumberGrid(1.0, 2500.0, 250)
p = plot_chebyshev_coefficients(M.k*M.χ*M.ST, [τrec], kgrid; tols = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7], yticks = 10.0 .^ (-6:5), ylims = (1e-6, 1e5))
#p = plot_chebyshev_coefficients(M.k*M.Sψ, [0.9*τ0], kgrid; ylims = (1e0, 1e8), label = nothing)
#savefig(p, "paper_chebyshev/figures/coeffs.pdf")

# Plot interpolation of matter power spectrum; using chebyshev; but for linear vs log k-grid
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
    kgrid = ChebyshevWavenumberGrid(ks_raw[begin], ks_raw[end], klength - 1; f = log, f⁻¹ = exp)
    Ss = source_grid(prob, log(S), [τ], ks_raw, kgrid, sol.bg; ptopts)
    Ss = exp.(stack(Ss)[1, :]) # time slice
    plot!(pmlin, ks_raw, Ss; color, ylabel = S, label = "$klength Chebyshev points", xformatter = _ -> "", bottom_margin = -5*Plots.mm, subplot = 1)
    plot!(pmlin, ks_raw, abs.(Ss./Ss_raw.-1); yscale = :log10, color, xlabel = M.k, ylabel = "error", ylims = (1e-6, 1e-1), label = nothing, subplot = 2)
end
pmlin
# TODO: now plot the same, but using a linear grid for k

# Plot Chebyshev points
function plot_chebyshev_points(; kwargs...)
    p = plot(xlabel = "Chebyshev points", ylabel = "number of points", xticks = ([-1, +1], ["kₘᵢₙ", "kₘₐₓ"]), grid = false, xlims = (-1, +1), ylims = (0, 80), size = (600, 300), margin = 2*Plots.mm, kwargs...)
    for klength in 2:80
        korder = klength - 1
        ks = chebgrid(-1, +1; order = korder)
        scatter!(p, ks, fill(length(ks), klength); label = nothing, color = :black, markersize = 2.5)
    end
    return p
end
p = plot_chebyshev_points()

# Compare CMB spectra with Chebyshev vs cubic spline interpolation
ls = range(log(2), log(2500), length=200) .|> exp .|> round |> unique .|> Int
jl = SphericalBesselCache(ls)
mode = :TT
xs = 0.0:0.0002:1.0
Δkτ0 = 2π/16 # fixes low-l jagginess
kmin, kmax = 0.01, 1000.0 # kτ0max = 2*lmax for TT/EE, 10*lmax for ψψ
normalization = :Dl
kgrid_ref = ChebyshevWavenumberGrid(kmin, kmax, 2^10)
klengths = [2^n+1 for n in 4:8]
kgrids_cubic = [CubicSplineWavenumberGrid(lingrid(kmin, kmax; length)) for length in klengths]
kgrids_cheb = [ChebyshevWavenumberGrid(kmin, kmax, length-1) for length in klengths]
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
f = k -> tanh(k/k0+1) # tanh(1) ≈ 0.76, tanh(2) ≈ 0.96
#f = k -> asinh(k/k0)
#γ = 1.0
#f = k -> asinh((k/k0)^γ)
#f⁻¹ = f -> k0 * (sinh(f))^(1/γ)
#f = k -> log(1 + (k/k0)^γ)
#f⁻¹ = f -> k0 * (exp(f) - 1)^(1/γ)
#f = k -> SymBoltz.smoothifelse(k/k0-1, k, k+log(k/k0); k=1.0)
#f = k -> ifelse(k<k0, k, k0+0.01(k-k0)) # not smooth
plot(f, xlims = (0, 1e4))
kgrid = ChebyshevWavenumberGrid(1e-2, 2e4, 250; f)
S = M.ST
τ = τrec
Sgrid = source_grid(prob, S, [τ], kgrid.ks)[1, :]
ks = lingrid(minimum(kgrid), maximum(kgrid); length = 5000)
Ss = source_grid(prob, S, [τ], ks, kgrid)[1, :]
p = begin
    color = :black
    xformatter = _ -> ""
    yformatter = _ -> ""
    p = plot(layout = (3, 1), size = (600, 1000), margin = 5*Plots.mm)
    plot!(p, kgrid.ks, kgrid.fs; color, marker = :circle, xlabel = "k", ylabel = "f(k)", label = nothing, subplot = 1)
    vline!(p, [k0], subplot = 1, color = :gray, linestyle = :dash, label = "k = k₀")
    hline!(p, [f(k0)], subplot = 1, color = :gray, linestyle = :dash, label = nothing)
    plot!(p, ks, Ss; color, xlabel = "k", ylabel = "S", xticks = kgrid.ks, label = nothing, xformatter, yformatter, subplot = 2)
    scatter!(p, kgrid.ks, Sgrid; color, xformatter, yformatter, label = nothing, subplot = 2)
    vline!(p, [k0], subplot = 2, color = :gray, linestyle = :dash, label = nothing)
    plot!(p, f.(ks), Ss; color, xlabel = "f(k)", ylabel = "S", label = nothing, xticks = kgrid.fs, xformatter, yformatter, subplot = 3)
    scatter!(p, kgrid.fs, Sgrid; color, xformatter, yformatter, label = nothing, subplot = 3)
    vline!(p, [f(k0)], subplot = 3, color = :gray, linestyle = :dash, label = nothing)
end

# Plot convergence (cₙ as a function of n) for fixed N, but different f(k)
S = M.k*M.χ*M.ST
p = plot(; layout = (3, 1), size = (600, 1000), left_margin = 5*Plots.mm, right_margin = 5*Plots.mm)
for (f, fname) in [
    (identity, "k"),
    (k->log(k/2000), "log(k/2000)"),
    (k->asinh(k/2000), "asinh(k/2000)"),
    (k->tanh(k/2000), "tanh(k/2000)"),
]
    kgrid = ChebyshevWavenumberGrid(1e-2, 1e4, 250; f)
    interp = only(source_grid_interp(prob, S, [τrec], kgrid, sol.bg))
    ns = 0:order(kgrid)
    cs = abs.(interp.coefs)
    plot!(p[1, 1], ns, cs; xlabel = "n", ylabel = "|cₙ|", yscale = :log10, yticks = 10.0 .^ (-5:5), title = "S(τ = $(round(τrec; digits=3)), k) = $S", label = nothing)

    fmin, fmax = extrema(kgrid.fs)
    plot!(p[2, 1], kgrid.ks, (kgrid.fs .- fmin) ./ (fmax - fmin); xlabel = "k", ylabel = "f(k)", label = "f(k) = $fname", xticks = [kgrid.ks[begin], 2000, kgrid.ks[end]], yticks = nothing)

    fs = range(fmin, fmax; length = 1000)
    Ss = interp.(fs)
    plot!(p[3, 1], (fs .- fmin) ./ (fmax - fmin), Ss; xlabel = "f(k)", ylabel = S, label = nothing, xticks = nothing, yticks = nothing)
end
p
