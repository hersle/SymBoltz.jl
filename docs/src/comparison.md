# Comparison to CLASS

This page compares results from SymBoltz and the fantastic Einstein-Boltzmann solver [CLASS](https://github.com/lesgourg/class_public/) for the ΛCDM model.

```@raw html
<details>
<summary><h2 style="display: inline-block">Setup</h2></summary>
```
```@example class
using SymBoltz
using ModelingToolkit
using DelimitedFiles
using DataInterpolations
using Unitful, UnitfulAstro
using CairoMakie

lmax = 6
M = w0waCDM(; lmax)
pars = merge(parameters_Planck18(M), Dict(
    M.X.w0 => -0.9,
    M.X.wa => 0.1,
    M.X.cₛ² => 0.9
))
h = pars[M.g.h] # needed later
prob = CosmologyProblem(M, pars)

function run_class(in::Dict{String, Any}, exec, inpath, outpath)
    merge!(in, Dict(
        "root" => outpath,
        "overwrite_root" => "yes",
    ))
    in = sort(collect(in), by = keyval -> keyval[1]) # sort by key
    in = prod(["$key = $val\n" for (key, val) in in]) # make string with "key = val" lines
    mkpath(dirname(inpath)) # create directories up to input file, if needed
    write(inpath, in) # write input string to file

    mkpath(dirname(outpath)) # create directories up to output file, if needed
    return run(`$exec $inpath`) # run
end

function solve_class(pars, k = nothing; exec="class", dir = mktempdir(), out = Set(["bg", "th", "pt", "P", "Cl"]))
    isnothing(k) && filter!(s -> s != "pt", out)
    inpath = joinpath(dir, "input.ini")
    outpath = joinpath(dir, "output", "")
    in = Dict(
        "write_background" => "yes",
        "write_thermodynamics" => "yes",
        "background_verbose" => 2,
        "output" => "mPk, tCl, pCl", # need one to evolve perturbations

        "k_output_values" => isnothing(k) ? "" : NoUnits(k / u"1/Mpc"),
        "ic" => "ad",
        "modes" => "s",
        "gauge" => "newtonian",

        # metric
        "h" => pars[M.g.h],

        # photons
        "T_cmb" => pars[M.γ.T₀],
        "l_max_g" => lmax,
        "l_max_pol_g" => lmax,

        # baryons
        "Omega_b" => pars[M.b.Ω₀],
        "YHe" => pars[M.b.rec.Yp],
        "recombination" => "recfast", # TODO: HyREC
        "recfast_Hswitch" => 1,
        "recfast_Heswitch" => 6,
        "reio_parametrization" => "reio_camb",

        # cold dark matter
        "Omega_cdm" => pars[M.c.Ω₀],

        # neutrinos
        "N_ur" => SymBoltz.have(M, :ν) ? pars[M.ν.Neff] : 0.0,
        "N_ncdm" => SymBoltz.have(M, :h) ? 1 : 0,
        "m_ncdm" => SymBoltz.have(M, :h) ? pars[M.h.m_eV] : 0.0, # in eV/c^2
        "T_ncdm" => SymBoltz.have(M, :h) ? (4/11)^(1/3) : 0.0,
        "l_max_ur" => lmax,
        "l_max_ncdm" => lmax,

        # primordial power spectrum
        "A_s" => pars[M.I.As],
        "n_s" => pars[M.I.ns],

        # w0wa dark energy
        "Omega_Lambda" => 0.0, # unspecified
        "w0_fld" => SymBoltz.have(M, :X) ? pars[M.X.w0] : -1.0,
        "wa_fld" => SymBoltz.have(M, :X) ? pars[M.X.wa] : 0.0,
        "cs2_fld" => SymBoltz.have(M, :X) ? pars[M.X.cₛ²] : 1.0,
        "use_ppf" => SymBoltz.have(M, :X) ? "no" : "yes", # use full equations (must be yes with CC to prevent crash)

        # other stuff
        "Omega_k" => SymBoltz.have(M, :K) ? pars[M.K.Ω₀] : 0.0, # curvature
        "Omega_scf" => 0.0,
        "Omega_dcdmdr" => 0.0,

        # approximations (see include/precisions.h)
        "tight_coupling_trigger_tau_c_over_tau_h" => 1e-2, # cannot turn off?
        "tight_coupling_trigger_tau_c_over_tau_k" => 1e-3, # cannot turn off
        "radiation_streaming_approximation" => 3, # turns off RSA
        "ur_fluid_approximation" => 3, # turns off UFA
        "ncdm_fluid_approximation" => 3, # turns off NCDM fluid approximation

        #"l_max_scalars" => 1500,
        #"temperature_contributions" => "pol",
    )

    secs = @elapsed run_class(in, exec, inpath, outpath)
    println("Ran CLASS in $secs seconds in directory $dir")
    output = Dict()
    for (name, filename, skipstart, target_length) in [("bg", "_background.dat", 3, typemax(Int)), ("th", "_thermodynamics.dat", 10, typemax(Int)), ("pt", "_perturbations_k0_s.dat", 1, typemax(Int)), ("P", "_pk.dat", 3, typemax(Int)), ("Cl", "_cl.dat", 6, typemax(Int))]
        name ∉ out && continue
        file = joinpath(outpath, filename)
        data, head = readdlm(file, skipstart=skipstart, header=true)
        head = split(join(head, ""), ":")
        for (n, h) in enumerate(head[begin:end-1])
            head[n] = h[begin:end-length(string(n))]
        end
        head = string.(head[2:end]) # remove #
        data = data[:,begin:length(head)] # remove extra empty data rows because of CLASS' messed up names
        data = Matrix{Float64}(data)
        println("Read ", filename, " with ", join(size(data), "x"), " numbers")
        output[name] = Dict(head[i] => data[:,i] for i in eachindex(head))
        for key in keys(output[name])
            output[name][key] = SymBoltz.reduce_array!(output[name][key], target_length)
        end
    end
    return output
end

k = 1e1 / u"Mpc" # 1/Mpc
sol1 = solve_class(pars, k)
sol2 = solve(prob, k; ptopts = (alg = SymBoltz.Rodas4P(),)) # looks like lower-precision KenCarp4 and Kvaerno5 "emulate" radiation streaming, while higher-precision Rodas5P continues in an exact way

function plot_compare(x1s, x2s, y1s, y2s, xlabel, ylabels; lgx=false, lgy=false, common=false, errtype=:auto, errlim=NaN, tol = nothing, kwargs...)
    if !(ylabels isa AbstractArray)
        y1s = [y1s]
        y2s = [y2s]
        ylabels = [ylabels]
    end

    xplot(x) = lgx ? log10.(abs.(x)) : x
    yplot(y) = lgy ? log10.(abs.(y)) : y
    xlab(x) = lgx ? "lg(|$x|)" : x
    ylab(y) = lgy ? "lg(|$y|)" : y

    fig = Figure(width = 800, height = 600)
    ax1 = Axis(fig[2:3, 1]; xticklabelsvisible = false)
    ax2 = Axis(fig[4, 1]; xlabel = xlab(xlabel), yminorticks = IntervalsBetween(10), yminorgridvisible = true)
    maxerr = 0.0
    xlims = (NaN, NaN)
    linewidth = 2
    lines!(ax1, [NaN]; color = :black, linewidth, linestyle = :solid, label = "CLASS") # dummy for legend entry
    lines!(ax1, [NaN]; color = :black, linewidth, linestyle = :dash, label = "SymBoltz") # dummy for legend entry
    for (i, (y1, y2, ylabel)) in enumerate(zip(y1s, y2s, ylabels))
        x1, x2 = x1s, x2s
        x1min, x1max = extrema(x1)
        x2min, x2max = extrema(x2)

        plotx1 = xplot(x1) # https://discourse.julialang.org/t/makie-axis-limits-with-inf/85784
        plotx2 = xplot(x2)
        ploty1 = replace!(yplot(y1), -Inf => NaN, Inf => NaN)
        ploty2 = replace!(yplot(y2), -Inf => NaN, Inf => NaN)
        color = Makie.wong_colors()[i]
        lines!(ax1, plotx1, ploty1; color, alpha = 0.6, linewidth, linestyle = :solid)
        lines!(ax1, plotx2, ploty2; color, alpha = 0.6, linewidth, linestyle = :dash)
        lines!(ax1, [NaN]; color, linewidth, label = ylab(ylabel)) # dummy for legend entry

        if i == 1
            xlims = (x1min, x1max) # initialize so not NaN
        end
        xlims = if common
            (max(xlims[1], x1min, x2min), min(xlims[2], x1max, x2max))
        else
            (min(xlims[1], x1min, x2min), max(xlims[2], x1max, x2max))
        end

        xmin = max(x1min, x2min) # need common times to compare error
        xmax = min(x1max, x2max) # need common times to compare error
        i1s = xmin .≤ x1 .≤ xmax # select x values in common range and exclude points too close in time (probably due to CLASS switching between approximation schemes)
        i2s = xmin .≤ x2 .≤ xmax # select x values in common range
        x1 = x1[i1s]
        x2 = x2[i2s]
        y1 = y1[i1s]
        y2 = y2[i2s]
        x = x2 # compare ratios at SymBoltz' times
        y1 = LinearInterpolation(y1, x1; extrapolation = ExtrapolationType.Linear).(x) # TODO: use built-in CosmoloySolution interpolation
        y2 = LinearInterpolation(y2, x2; extrapolation = ExtrapolationType.Linear).(x)

        !isnothing(tol) && @assert all(isapprox.(y1, y2; atol = tol)) "$ylabel does not match within absolute tolerance $tol. Maximum difference was $(maximum(abs.(y1.-y2)))."

        # Compare absolute error if quantity crosses zero, otherwise relative error (unless overridden)
        abserr = (errtype == :abs) || (errtype == :auto && (any(y1 .<= 0) || any(y2 .<= 0)))
        if abserr
            err = y2 .- y1
            ylabel = "SymBoltz - CLASS"
        else
            err = y2 ./ y1 .- 1
            ylabel = "SymBoltz / CLASS - 1"
        end
        ax2.ylabel = ylabel
        maxerr = max(maxerr, maximum(abs.(err)))
        lines!(ax2, xplot(x), err; color, alpha = 1.0, linewidth)
    end

    # write maximum error like a * 10^b (for integer a and b)
    if isnan(errlim)
        b = floor(log10(maxerr))
        a = ceil(maxerr / 10^b)
        errlim = a * 10^b
    end
    Makie.ylims!(ax2, (-errlim, +errlim))

    xlims = xplot.(xlims) # transform to log?
    Makie.xlims!(ax1, xlims)
    Makie.xlims!(ax2, xlims)

    #axislegend(ax1)
    fig[1, 1] = Legend(fig, ax1; padding = (0, 0, 0, 0), margin = (0, 0, 0, 0), framevisible = false, orientation = :horizontal)

    return fig
end
nothing # hide
```
```@raw html
</details>
```

## Background

### Conformal time
```@example class
a1 = (1 ./ (sol1["bg"]["z"] .+ 1))
a2 = sol2[M.g.a]
χ1 = sol1["bg"]["conf.time[Mpc]"][end].-sol1["bg"]["conf.time[Mpc]"]
χ2 = sol2[M.χ] / (h*SymBoltz.k0)
plot_compare(a1, a2, χ1, χ2, "a", "χ"; tol = 1e-2)
```
### Hubble function
```@example class
E1 = sol1["bg"]["H[1/Mpc]"]./sol1["bg"]["H[1/Mpc]"][end]
E2 = sol2[M.g.E]
plot_compare(a1, a2, E1, E2, "a", "E"; lgx=true, lgy=true, tol = 1e9)
```
### Energy densities
```@example class
ρ1 = map(s -> sol1["bg"]["(.)rho_$s"], ["g", "ur", "cdm", "b", "fld", "ncdm[0]"])
ρ2 = map(s -> sol2[s.ρ] * 8π/3*(h*SymBoltz.k0)^2, [M.γ, M.ν, M.c, M.b, M.X, M.h])
plot_compare(a1, a2, ρ1, ρ2, "a", ["ργ", "ρb", "ρc", "ρX", "ρν", "ρh"]; lgx=true, lgy=true, tol = 1e16)
```
### Equations of state
```@example class
wh1 = sol1["bg"]["(.)p_ncdm[0]"] ./ sol1["bg"]["(.)rho_ncdm[0]"]
wh2 = sol2[M.h.w]
wX1 = sol1["bg"]["(.)w_fld"]
wX2 = sol2[M.X.w]
plot_compare(a1, a2, [wh1, wX1], [wh2, wX2], "a", ["wh", "wX"]; lgx=true, tol = 1e-3)
```
### Photon-baryon sound horizon
```@example class
rs1 = sol1["bg"]["comov.snd.hrz."]
rs2 = sol2[M.rₛ] ./ (h*SymBoltz.k0)
plot_compare(a1, a2, rs1, rs2, "a", "rₛ"; lgx = true, tol = 1e-2)
```
### Luminosity distance
```@example class
dL1 = sol1["bg"]["lum.dist."]
dL2 = SymBoltz.distance_luminosity(sol2) / SymBoltz.Mpc
plot_compare(a1, a2, dL1, dL2, "a", "dL"; lgx=true, lgy=true)
```

## Thermodynamics

### Optical depth derivative
```@example class
a1 = reverse(sol1["th"]["scalefactora"])
a2 = sol2[M.g.a]
dκ1 = reverse(sol1["th"]["kappa'[Mpc^-1]"])
dκ2 = -sol2[M.b.rec.κ̇] * (h*SymBoltz.k0)
plot_compare(a1, a2, dκ1, dκ2, "a", "κ̇"; lgx=true, lgy=true, tol = 1e4)
```
### Optical depth exponential
```@example class
expmκ1 = reverse(sol1["th"]["exp(-kappa)"])
expmκ2 = sol2[M.b.rec.I]
plot_compare(a1, a2, expmκ1, expmκ2, "a", "exp(-κ)"; lgx=true, tol = 1e-3)
```
### Visibility function
```@example class
v1 = reverse(sol1["th"]["g[Mpc^-1]"])
v2 = sol2[M.b.rec.v] * (h*SymBoltz.k0)
plot_compare(a1, a2, v1, v2, "a", "v"; lgx=true, lgy=false, tol = 1e-4)
```
### Free electron fraction
```@example class
Xe1 = reverse(sol1["th"]["x_e"])
Xe2 = sol2[M.b.rec.Xe]
plot_compare(a1, a2, Xe1, Xe2, "a", "Xe"; lgx=true, lgy=false, tol = 1e-3)
```
### Baryon temperature
```@example class
Tb1 = reverse(sol1["th"]["Tb[K]"])
Tb2 = sol2[M.b.rec.Tb]
dTb1 = reverse(sol1["th"]["dTb[K]"])
dTb2 = sol2[M.b.rec.DTb] ./ -sol2[M.g.E] # convert my dT/dt̂ to CLASS' dT/dz = -1/H * dT/dt
plot_compare(a1, a2, [Tb1, dTb1], [Tb2, dTb2], "a", ["Tb", "dTb"]; lgx=true, lgy=true, tol = 1e1)
```
### Baryon equation of state
```@example class
# baryon equation of state parameter (e.g. https://arxiv.org/pdf/1906.06831 eq. (B10))
wb1 = reverse(sol1["th"]["w_b"])
wb2 = sol2[SymBoltz.kB*M.b.rec.Tb/(M.b.rec.μ*SymBoltz.c^2)]
plot_compare(a1, a2, wb1, wb2, "a", "wb"; lgx=true, lgy=true, tol = 1e-9)
```
### Baryon sound speed
```@example class
csb²1 = reverse(sol1["th"]["c_b^2"])
csb²2 = sol2[M.b.rec.cₛ²]
plot_compare(a1, a2, csb²1, csb²2, "a", "csb²"; lgx=true, lgy=true, tol = 1e-8)
```

## Perturbations

### Metric potentials
```@example class
a1 = sol1["pt"]["a"]
a2 = sol2[1, M.g.a]
Φ1, Ψ1 = sol1["pt"]["phi"], sol1["pt"]["psi"]
Φ2, Ψ2 = sol2[1, M.g.Φ], sol2[1, M.g.Ψ]
plot_compare(a1, a2, [Φ1, Ψ1], [Φ2, Ψ2], "a", ["Ψ", "Φ"]; lgx=true, tol = 1e-1)
```
### Energy overdensities
```@example class
δ1 = map(s -> sol1["pt"]["delta_$s"], ["b", "cdm", "g", "ur", "ncdm[0]"])
δ2 = map(s -> sol2[1, s.δ], [M.b, M.c, M.γ, M.ν, M.h])
plot_compare(a1, a2, δ1, δ2, "a", ["δb", "δc", "δγ", "δν", "δh"]; lgx=true, lgy=true)
```
### Momenta
```@example class
θ1 = map(s -> sol1["pt"]["theta_$s"], ["b", "cdm", "g", "ur", "ncdm[0]"])
θ2 = map(s -> sol2[1, s.θ] * (h*SymBoltz.k0), [M.b, M.c, M.γ, M.ν, M.h])
plot_compare(a1, a2, θ1, θ2, "a", ["θb", "θc", "θγ", "θν", "θh"]; lgx=true, lgy=true)
```
### Dark energy overdensity
```@example class
δρX1 = sol1["pt"]["delta_rho_fld"]
δρX2 = sol2[1, M.X.δ*M.X.ρ] * 8π/3*(h*SymBoltz.k0)^2
plot_compare(a1, a2, δρX1, δρX2, "a", "δρX"; lgx=true, lgy=true, tol = 1e-4)
```
### Dark energy momentum
```@example class
pX1 = sol1["pt"]["rho_plus_p_theta_fld"]
pX2 = sol2[1, (M.X.ρ+M.X.P)*M.X.θ * 8π/3*(h*SymBoltz.k0)^3]
plot_compare(a1, a2, pX1, pX2, "a", "pX"; lgx=true, lgy=true, tol = 1e-4)
```
### Shear stresses
```@example class
σ1 = [sol1["pt"]["shear_g"], sol1["pt"]["shear_ur"]]
σ2 = [sol2[1, M.γ.σ], sol2[1, M.ν.F[2]/2]]
plot_compare(a1, a2, σ1, σ2, "a", ["σγ", "σν"]; lgx=true, tol = 1e-0)
```
### Polarization
```@example class
P1 = map(n -> sol1["pt"]["pol$(n)_g"], 0:2)
P2 = [sol2[1, var] for var in [M.γ.G0, M.γ.G[1], M.γ.G[2]]]
plot_compare(a1, a2, P1, P2, "a", ["P0", "P1", "P2"]; lgx=true, tol = 1e-1)
```

## Matter power spectrum
```@example class
function P_class(pars)
    sol = solve_class(pars; out = ["P"])
    h = pars[M.g.h]
    k = sol["P"]["k(h/Mpc)"] * h
    P = sol["P"]["P(Mpc/h)^3"] / h^3
    return k, P
end
function P_class(k, pars)
    k′, P′ = P_class(pars)
    P = exp.(LinearInterpolation(log.(P′), log.(k′)).(log.(k)))
    return P
end
function P_symboltz(k, pars)
    prob′ = SymBoltz.parameter_updater(prob, collect(keys(pars)))(collect(values(pars))) # TODO: move outside; common for Pk and Cl
    P = spectrum_matter(prob′, k / u"Mpc") / u"Mpc^3"
    return P
end
k, P1 = P_class(pars)
P2 = P_symboltz(k, pars)
plot_compare(k, k, P1, P2, "k/Mpc⁻¹", "P/Mpc³"; lgx = true, lgy = true, tol = 1e2)
```
```@example class
using ForwardDiff, FiniteDiff

k = 10 .^ range(-3, 0, length=100) # 1/Mpc
function plot_compare_P_diff(par, val; relstep = 1e-2, kwargs...)
    out = zeros(length(k))
    f = function(out, logval)
        val = exp(logval)
        out .= log.(P_class(k, merge(pars, Dict(par => val))))
        return out
    end
    Δt1 = @elapsed ∂logP1_∂logθ = FiniteDiff.finite_difference_gradient!(out, f, log(val), Val{:central}; relstep)
    println("Computed CLASS derivatives in $Δt1 seconds")
    Δt2 = @elapsed ∂logP2_∂logθ = ForwardDiff.derivative(logval -> log.(P_symboltz(k, Dict(par => exp(logval)))), log(val))
    println("Computed SymBoltz derivatives in $Δt2 seconds")
    return plot_compare(k, k, ∂logP1_∂logθ, ∂logP2_∂logθ, "k", "∂(log(P))/∂(log($(replace(string(par), "₊" => "."))))"; lgx = true, kwargs...)
end

#∂logP1_∂θ = FiniteDiff.finite_difference_jacobian(θ -> log.(P_class(merge(pars, Dict(diffpars .=> θ)))[2]), θ, Val{:central}; relstep = 1e-4) # hide
#∂logP2_∂θ = ForwardDiff.jacobian(θ -> log.(P_symboltz(k, Dict(diffpars .=> θ))[2]), θ) # hide
#plot_compare(k, k, eachcol(∂logP1_∂θ), eachcol(∂logP2_∂θ), "k", ["∂(lg(P))/∂($par)" for par in ["Ωc0", "Ωb0", "h"]]; lgx = true) # hide

plot_compare_P_diff(M.c.Ω₀, pars[M.c.Ω₀]; relstep = 1e-3, tol = 1e-1) # smaller relstep is noisier
```
```@example class
plot_compare_P_diff(M.b.Ω₀, pars[M.b.Ω₀]; relstep = 1e-3, tol = 1e-1) # smaller relstep is noisier
```
```@example class
plot_compare_P_diff(M.g.h, pars[M.g.h]; relstep = 1e-3, tol = 1e-1) # smaller relstep is noisier
```

## CMB power spectrum
```@example class
function Dl_class(modes, l, pars)
    sol = solve_class(pars; out = ["Cl"])
    lout = sol["Cl"]["l"]
    Dl = [sol["Cl"][string(mode)] for mode in modes]
    i = findall(l′ -> l′ in l, lout)
    Dl = [Dl[j][i] for j in eachindex(modes)]
    return stack(Dl)
end
function Dl_symboltz(modes, l, pars; kwargs...)
    prob′ = SymBoltz.parameter_updater(prob, collect(keys(pars)))(collect(values(pars)))
    return spectrum_cmb(modes, prob′, l; normalization = :Dl, kwargs...)
end

l = 20:20:2000 # CLASS default is lmax = 2500
Dl1 = Dl_class([:TT, :TE, :EE], l, pars)
Dl2 = Dl_symboltz([:TT, :TE, :EE], l, pars)
plot_compare(l, l, Dl1[:, 1], Dl2[:, 1], "l", "Dₗ(TT)"; tol = 7e-13)
```
```@example class
plot_compare(l, l, Dl1[:, 2], Dl2[:, 2], "l", "Dₗ(TE)"; tol = 8e-14)
```
```@example class
plot_compare(l, l, Dl1[:, 3], Dl2[:, 3], "l", "Dₗ(EE)"; tol = 2e-14)
```
```@example class
diffpars = [M.c.Ω₀, M.b.Ω₀, M.g.h]
θ = [pars[par] for par in diffpars]
modes = [:TT, :TE, :EE]

Δt1 = @elapsed ∂Dl1_∂θ = FiniteDiff.finite_difference_jacobian(θ -> Dl_class(modes, l, merge(pars, Dict(diffpars .=> θ))), θ, Val{:central}; relstep = 1e-3)
println("Computed CLASS derivatives in $Δt1 seconds")
Δt2 = @elapsed ∂Dl2_∂θ = ForwardDiff.jacobian(θ -> Dl_symboltz(modes, l, Dict(diffpars .=> θ)), θ)
println("Computed SymBoltz derivatives in $Δt2 seconds")

# returned matrices have size (length(l)*length(modes), length(diffpars)); reshape to (length(l), length(modes), length(diffpars))
∂Dl1_∂θ_3d = reshape(∂Dl1_∂θ, (length(l), length(modes), length(diffpars)))
∂Dl2_∂θ_3d = reshape(∂Dl2_∂θ, (length(l), length(modes), length(diffpars)))

plot_compare(l, l, eachcol(∂Dl1_∂θ_3d[:,1,:]), eachcol(∂Dl2_∂θ_3d[:,1,:]), "l", ["∂(Dₗ)/∂($(replace(string(par), "₊" => "."))) (TT)" for par in diffpars]; tol = 1e-10)
```
```@example class
plot_compare(l, l, eachcol(∂Dl1_∂θ_3d[:,2,:]), eachcol(∂Dl2_∂θ_3d[:,2,:]), "l", ["∂(Dₗ)/∂($(replace(string(par), "₊" => "."))) (TE)" for par in diffpars]; tol = 1e-11)
```
```@example class
plot_compare(l, l, eachcol(∂Dl1_∂θ_3d[:,3,:]), eachcol(∂Dl2_∂θ_3d[:,3,:]), "l", ["∂(Dₗ)/∂($(replace(string(par), "₊" => "."))) (EE)" for par in diffpars]; tol = 1e-12)
```

