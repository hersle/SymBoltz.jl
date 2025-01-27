# Comparison to CLASS

[This comparison script](https://github.com/hersle/SymBoltz.jl/blob/main/docs/src/comparison.md) automatically builds the latest version of the fantastic Einstein-Boltzmann solver [CLASS](https://github.com/lesgourg/class_public/) and compares its results to the latest version of SymBoltz:

```@example class
using SymBoltz
lmax = 6
M = SymBoltz.ΛCDM(; lmax)
pars = SymBoltz.parameters_Planck18(M)
prob = CosmologyProblem(M, pars; tspan = (1e-6, 1e3))
```
```@setup class
using ModelingToolkit
using DelimitedFiles
using DataInterpolations
using Unitful, UnitfulAstro
using Plots

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

function solve_class(pars, k; exec="class", inpath="/tmp/symboltz_class/input.ini", outpath="/tmp/symboltz_class/output/")
    k = NoUnits(k / u"1/Mpc")
    in = Dict(
        "write_background" => "yes",
        "write_thermodynamics" => "yes",
        "background_verbose" => 2,
        "output" => "mPk, tCl", # need one to evolve perturbations

        "k_output_values" => k,
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
        "N_ur" => pars[M.ν.Neff],
        "N_ncdm" => SymBoltz.have(M, :h) ? 1 : 0,
        "m_ncdm" => SymBoltz.have(M, :h) ? pars[M.h.m] / (SymBoltz.eV/SymBoltz.c^2) : 0.0, # in eV/c^2
        "T_ncdm" => SymBoltz.have(M, :h) ? (4/11)^(1/3) : 0.0, # TODO: CLASS uses something slightly different by default
        "l_max_ur" => lmax,
        "l_max_ncdm" => lmax,

        # primordial power spectrum
        "A_s" => pars[M.I.As],
        "n_s" => pars[M.I.ns],

        # other stuff
        "Omega_k" => 0.0,
        "Omega_fld" => 0.0,
        "Omega_scf" => 0.0,
        "Omega_dcdmdr" => 0.0,

        # approximations (see include/precisions.h)
        "tight_coupling_trigger_tau_c_over_tau_h" => 1e-2, # cannot turn off?
        "tight_coupling_trigger_tau_c_over_tau_k" => 1e-3, # cannot turn off
        #"radiation_streaming_approximation" => 3, # turns off RSA; commented because the number of perturbation points explodes without RSA
        "ur_fluid_approximation" => 3, # turns off UFA

        #"l_max_scalars" => 1500,
        #"temperature_contributions" => "pol", # TODO: make agree with polarization
    )

    run_class(in, exec, inpath, outpath)
    output = Dict()
    for (name, filename, skipstart, target_length) in [("bg", "_background.dat", 3, typemax(Int)), ("th", "_thermodynamics.dat", 10, typemax(Int)), ("pt", "_perturbations_k0_s.dat", 1, 25000), ("P", "_pk.dat", 3, typemax(Int)), ("Cl", "_cl.dat", 6, typemax(Int))]
        file = outpath * filename
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
sol2 = solve(prob, k; ptopts = (alg = SymBoltz.Rodas5P(),)) # looks like lower-precision KenCarp4 and Kvaerno5 "emulate" radiation streaming, while higher-precision Rodas5P continues in an exact way

# map results from both codes to common convention
h = pars[M.g.h]
sols = Dict(
    # background
    "a_bg" => (1 ./ (sol1["bg"]["z"] .+ 1), sol2[M.g.a]),
    "t" => (sol1["bg"]["conf.time[Mpc]"], sol2[M.t] / (h * SymBoltz.k0)),
    "E" => (sol1["bg"]["H[1/Mpc]"] ./ sol1["bg"]["H[1/Mpc]"][end], sol2[M.g.E]),
    "ργ" => (sol1["bg"]["(.)rho_g"] / sol1["bg"]["(.)rho_crit"][end], sol2[M.γ.ρ] / (3/8π)),
    "ρν" => (sol1["bg"]["(.)rho_ur"] / sol1["bg"]["(.)rho_crit"][end], sol2[M.ν.ρ] / (3/8π)),
    "ρc" => (sol1["bg"]["(.)rho_cdm"] / sol1["bg"]["(.)rho_crit"][end], sol2[M.c.ρ] / (3/8π)),
    "ρb" => (sol1["bg"]["(.)rho_b"] / sol1["bg"]["(.)rho_crit"][end], sol2[M.b.ρ] / (3/8π)),
    "ρΛ" => (sol1["bg"]["(.)rho_lambda"] / sol1["bg"]["(.)rho_crit"][end], sol2[M.Λ.ρ] / (3/8π)),
    "ρh" => (sol1["bg"]["(.)rho_ncdm[0]"] / sol1["bg"]["(.)rho_crit"][end], sol2[M.h.ρ] / (3/8π)),
    "wh" => (sol1["bg"]["(.)p_ncdm[0]"] ./ sol1["bg"]["(.)rho_ncdm[0]"], sol2[M.h.w]),

    # thermodynamics
    "a_th" => (reverse(sol1["th"]["scalefactora"]), sol2[M.g.a]),
    "τ̇" => (reverse(sol1["th"]["kappa'[Mpc^-1]"]), -sol2[M.b.rec.τ̇] * (SymBoltz.k0*h)),
    "csb²" => (reverse(sol1["th"]["c_b^2"]), sol2[M.b.rec.cₛ²]),
    "Xe" => (reverse(sol1["th"]["x_e"]), sol2[M.b.rec.Xe]),
    "Tb" => (reverse(sol1["th"]["Tb[K]"]), sol2[M.b.rec.Tb]),
    "exp(-τ)" => (reverse(sol1["th"]["exp(-kappa)"]), sol2[exp(-M.b.rec.τ)]),
    "v" => (reverse(sol1["th"]["g[Mpc^-1]"]), sol2[M.b.rec.v] * (SymBoltz.k0*h)),
    #"Tb′" => (reverse(sol1["th"]["dTb[K]"]), sol2[M.b.rec.DTb] ./ -sol2[M.g.E]), # convert my dT/dt̂ to CLASS' dT/dz = -1/H * dT/dt 

    # perturbations
    "a_pt" => (sol1["pt"]["a"], sol2[1, M.g.a]),
    "Φ" => (sol1["pt"]["phi"], sol2[1, M.g.Φ]),
    "Ψ" => (sol1["pt"]["psi"], sol2[1, M.g.Ψ]),
    "δb" => (sol1["pt"]["delta_b"], sol2[1, M.b.δ]),
    "δc" => (sol1["pt"]["delta_cdm"], sol2[1, M.c.δ]),
    "δγ" => (sol1["pt"]["delta_g"], sol2[1, M.γ.δ]),
    "δν" => (sol1["pt"]["delta_ur"], sol2[1, M.ν.δ]),
    "δh" => (sol1["pt"]["delta_ncdm[0]"], sol2[1, M.h.δ]),
    "θb" => (sol1["pt"]["theta_b"], sol2[1, M.b.θ] * (h*SymBoltz.k0)),
    "θc" => (sol1["pt"]["theta_cdm"], sol2[1, M.c.θ] * (h*SymBoltz.k0)),
    "θγ" => (sol1["pt"]["theta_g"], sol2[1, M.γ.θ] * (h*SymBoltz.k0)),
    "θν" => (sol1["pt"]["theta_ur"], sol2[1, M.ν.θ] * (h*SymBoltz.k0)),
    "θh" => (sol1["pt"]["theta_ncdm[0]"], sol2[1, M.h.θ] * (h*SymBoltz.k0)),
    "σγ" => (sol1["pt"]["shear_g"], sol2[1, M.γ.F[2] / 2]),
    "σν" => (sol1["pt"]["shear_ur"], sol2[1, M.ν.F[2] / 2]),
    "P0" => (sol1["pt"]["pol0_g"], sol2[1, M.γ.G[0]]),
    "P1" => (sol1["pt"]["pol1_g"], sol2[1, M.γ.G[1]]),
    "P2" => (sol1["pt"]["pol2_g"], sol2[1, M.γ.G[2]]),
)

# matter power spectrum
ks = sol1["P"]["k(h/Mpc)"] * h # 1/Mpc
Ps_class = sol1["P"]["P(Mpc/h)^3"] / h^3
Ps = spectrum_matter(prob, ks / u"Mpc"; verbose=true) / u"Mpc^3"

# CMB power spectrum
ls = sol1["Cl"]["l"]
ls = Int.(ls[begin:10:end])
Dls_class = sol1["Cl"]["TT"]
Dls_class = Dls_class[begin:10:end]
Cls = spectrum_cmb(:TT, prob, ls; verbose=true)
Dls = SymBoltz.Dl(Cls, ls)

sols = merge(sols, Dict(
    "k" => (ks, ks),
    "P" => (Ps_class, Ps),
    "l" => (ls, ls),
    "Dl" => (Dls_class, Dls)
))

function plot_compare(xlabel, ylabels; lgx=false, lgy=false, errtype=:auto, errlim=NaN, alpha=1.0, kwargs...)
    if !(ylabels isa AbstractArray)
        ylabels = [ylabels]
    end

    xplot(x) = lgx ? log10.(abs.(x)) : x
    yplot(y) = lgy ? log10.(abs.(y)) : y
    xlab(x) = lgx ? "lg(|$x|)" : x
    ylab(y) = lgy ? "lg(|$y|)" : y

    p = plot(; layout=grid(2, 1, heights=(3/4, 1/4)), size = (800, 600), kwargs...)
    plot!(p[1]; titlefontsize = 8, ylabel = join(ylab.(ylabels), ", "))
    maxerr = 0.0
    for (i, ylabel) in enumerate(ylabels)
        x1, x2 = sols[xlabel]
        x1min, x1max = extrema(x1)
        x2min, x2max = extrema(x2)
        xmin = max(x1min, x2min)
        xmax = min(x1max, x2max)
        i1s = xmin .≤ x1 .≤ xmax # select x values in common range and exclude points too close in time (probably due to CLASS switching between approximation schemes)
        i2s = xmin .≤ x2 .≤ xmax # select x values in common range
        x1 = x1[i1s]
        x2 = x2[i2s]
        x = x2 # compare ratios at SymBoltz' times

        y1, y2 = sols[ylabel]
        y1 = y1[i1s]
        y2 = y2[i2s]
        plot!(p[1], xplot(x1), yplot(y1); color = :grey, linewidth = 2, alpha, linestyle = :solid, label = i == 1 ? "CLASS" : nothing, xformatter = _ -> "")
        plot!(p[1], xplot(x2), yplot(y2); color = :black, linewidth = 2, alpha, linestyle = :dash,  label = i == 1 ? "SymBoltz" : nothing)

        # TODO: use built-in CosmoloySolution interpolation
        y1 = LinearInterpolation(y1, x1; extrapolation = ExtrapolationType.Linear).(x)
        y2 = LinearInterpolation(y2, x2; extrapolation = ExtrapolationType.Linear).(x)

        # Compare absolute error if quantity crosses zero, otherwise relative error (unless overridden)
        abserr = (errtype == :abs) || (errtype == :auto && (any(y1 .<= 0) || any(y2 .<= 0)))
        if abserr
            err = y2 .- y1
            ylabel = "SymBoltz - CLASS"
        else
            err = y2 ./ y1 .- 1
            ylabel = "SymBoltz / CLASS - 1"
        end
        maxerr = max(maxerr, maximum(abs.(err)))
        plot!(p[end], xplot(x), err; linewidth = 2, yminorticks = 10, yminorgrid = true, color = :black, xlabel = xlab(xlabel), ylabel, top_margin = -5*Plots.mm, label = nothing)
    end

    # write maximum error like a * 10^b (for integer a and b)
    if isnan(errlim)
        b = floor(log10(maxerr))
        a = ceil(maxerr / 10^b)
        errlim = a * 10^b
    end
    plot!(p[end], ylims = (-errlim, +errlim))

    return p
end
nothing # hide
```

### Background

```@example class
plot_compare("a_bg", "t"; lgx=true, lgy=true) # hide
```
```@example class
plot_compare("a_bg", "E"; lgx=true, lgy=true) # hide
```
```@example class
plot_compare("a_bg", ["ργ", "ρb", "ρc", "ρΛ", "ρν", "ρh"]; lgx=true, lgy=true) # hide
```
```@example class
plot_compare("a_bg", ["wh"]; lgx=true) # hide
```

### Thermodynamics
```@example class
plot_compare("a_th", "τ̇"; lgx=true, lgy=true) # hide
```
```@example class
plot_compare("a_th", "exp(-τ)"; lgx=true) # hide
```
```@example class
plot_compare("a_th", "v"; lgx=true, lgy=false) # hide
```
```@example class
plot_compare("a_th", "Xe"; lgx=true, lgy=false) # hide
```
```@example class
plot_compare("a_th", "Tb"; lgx=true, lgy=true) # hide
```
```@example class
plot_compare("a_th", "csb²"; lgx=true, lgy=true) # hide
```

### Perturbations

```@example class
plot_compare("a_pt", ["Ψ", "Φ"]; lgx=true) # hide
```
```@example class
plot_compare("a_pt", ["δb", "δc", "δγ", "δν", "δh"]; lgx=true, lgy=true) # hide
```
```@example class
plot_compare("a_pt", ["θb", "θc", "θγ", "θν", "θh"]; lgx=true, lgy=true) # hide
```
```@example class
plot_compare("a_pt", ["σγ", "σν"]; lgx=true) # hide
```
```@example class
plot_compare("a_pt", ["P0", "P1", "P2"]; lgx=true) # hide
```

### Power spectrum

```@example class
plot_compare("k", "P"; lgx=true, lgy=true) # hide
```
```@example class
plot_compare("l", "Dl"; errlim=5e-2) # TODO: fix # hide
```
