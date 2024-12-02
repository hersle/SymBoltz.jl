# Comparison to other codes

TODO: rethink table

| Feature                  | [SymBoltz.jl](https://github.com/hersle/SymBoltz.jl) | [CAMB](https://camb.info/)    | [CLASS](https://lesgourg.github.io/class_public/class.html) |
| :----------------------: | :--------------------------------------------------: | :---------------------------: | :---------------------------------------------------------: |
| **Language**             | Julia                                                | Fortran + Python              | C + Python                                                  |
| **Modularity**           | Every component is fully modular                     | Class-based for simple models | Requires raw code modification in many places               |
| **Approximations**       | None                                                 | Mandatory (?)                 | Mandatory (?)                                               |
| **Speed**                | Fast                                                 | Faster                        | Fastest                                                     |
| **Sensitivity analysis** | Automatic differentiation                            | Finite differences            | Finite differences                                          |

## Numerical comparison to CLASS

[This comparison script](https://github.com/hersle/SymBoltz.jl/blob/main/docs/src/comparison.md) automatically builds the latest version of the fantastic Einstein-Boltzmann solver [CLASS](https://github.com/lesgourg/class_public/) and compares its results to the latest version of SymBoltz:

```@setup class
using SymBoltz
using ModelingToolkit
using DelimitedFiles
using DataInterpolations
using Unitful, UnitfulAstro
using Plots; Plots.default(label=nothing)
using Printf

lmax = 6
M = SymBoltz.ΛCDM(; lmax, h = nothing, Λanalytical = true)
pars = SymBoltz.parameters_Planck18(M)

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
        "reio_parametrization" => "reio_none", # TODO: enable

        # cold dark matter
        "Omega_cdm" => pars[M.c.Ω₀],

        # neutrinos # TODO: set neutrino stuff to 0 unless otherwise specified
        "N_ur" => pars[M.ν.Neff],
        "N_ncdm" => 0.0, # TODO
        "m_ncdm" => 0.00, # TODO
        "T_ncdm" => 0.0, # TODO (pars.Neff/3)^(1/4) * (4/11)^(1/3), # TODO: set properly
        #"Omega_ur" => 7/8 * (pars[M.ν.Neff]/3) * (4/11)^(4/3) * pars.Ωγ0, # TODO: set properly! # massless neutrinos # TODO: proper Neff
        "l_max_ur" => lmax,
        "l_max_ncdm" => lmax,

        # primordial power spectrum
        "A_s" => 2e-9, # TODO
        "n_s" => 1.0, # TODO

        # other stuff
        "Omega_k" => 0.0,
        "Omega_fld" => 0.0,
        "Omega_scf" => 0.0,
        "Omega_dcdmdr" => 0.0,

        # approximations (see include/precisions.h)
        "tight_coupling_trigger_tau_c_over_tau_h" => 1e-2, # cannot turn off?
        "tight_coupling_trigger_tau_c_over_tau_k" => 1e-3, # cannot turn off
        "radiation_streaming_approximation" => 3, # turn off
        "ur_fluid_approximation" => 3, # turn off
    )

    run_class(in, exec, inpath, outpath)
    output = Dict()
    for (name, filename, skipstart) in [("bg", "_background.dat", 3), ("th", "_thermodynamics.dat", 10), ("pt", "_perturbations_k0_s.dat", 1), ("P", "_pk.dat", 3), ("Cl", "_cl.dat", 6)]
        file = outpath * filename
        data, head = readdlm(file, skipstart=skipstart, header=true)
        head = split(join(head, ""), ":")
        for (n, h) in enumerate(head[begin:end-1])
            head[n] = h[begin:end-length(string(n))]
        end
        head = string.(head[2:end]) # remove #
        data = data[:,begin:length(head)] # remove extra empty data rows because of CLASS' messed up names
        data = Matrix{Float64}(data)
        output[name] = Dict(head[i] => data[:,i] for i in eachindex(head))
    end
    return output
end

k = 1e0 / u"Mpc" # 1/Mpc # disagreement on smaller scales
sol1 = solve_class(pars, k)
sol2 = solve(M, pars, k) # looks like lower-precision KenCarp4 and Kvaerno5 "emulate" radiation streaming, while higher-precision Rodas5P continues in an exact way

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
    #"ρmν" => (sol1["bg"]["(.)rho_ncdm[0]"] / sol1["bg"]["(.)rho_crit"][end], sol2[M.h.ρ] / (3/8π)),

    # thermodynamics
    "a_th" => (reverse(sol1["th"]["scalefactora"]), sol2[M.g.a]),
    "τ̇" => (reverse(sol1["th"]["kappa'[Mpc^-1]"]), .- sol2[M.b.rec.τ̇] * (SymBoltz.k0 * h)),
    "csb²" => (reverse(sol1["th"]["c_b^2"]), sol2[M.b.rec.cₛ²]),
    "Xe" => (reverse(sol1["th"]["x_e"]), sol2[M.b.rec.Xe]),
    "Tb" => (reverse(sol1["th"]["Tb[K]"]), sol2[M.b.rec.Tb]),
    #"Tb′" => (reverse(sol1["th"]["dTb[K]"]), sol2[M.b.rec.DTb] ./ -sol2[M.g.E]), # convert my dT/dt̂ to CLASS' dT/dz = -1/H * dT/dt 

    # perturbations
    "a_pt" => (sol1["pt"]["a"], sol2[1, M.g.a]),
    "Φ" => (sol1["pt"]["phi"], sol2[1, M.g.Φ]),
    "Ψ" => (sol1["pt"]["psi"], sol2[1, M.g.Ψ]),
    "δb" => (sol1["pt"]["delta_b"], sol2[1, M.b.δ]),
    "δc" => (sol1["pt"]["delta_cdm"], sol2[1, M.c.δ]),
    "δγ" => (sol1["pt"]["delta_g"], sol2[1, M.γ.δ]),
    "δν" => (sol1["pt"]["delta_ur"], sol2[1, M.ν.δ]),
    #"δmν" => (sol1["pt"]["delta_ncdm[0]"], sol2[1, M.h.δ]),
    "θb" => (sol1["pt"]["theta_b"], sol2[1, M.b.θ] * (h*SymBoltz.k0)),
    "θc" => (sol1["pt"]["theta_cdm"], sol2[1, M.c.θ] * (h*SymBoltz.k0)),
    "θγ" => (sol1["pt"]["theta_g"], sol2[1, M.γ.θ] * (h*SymBoltz.k0)),
    "θν" => (sol1["pt"]["theta_ur"], sol2[1, M.ν.θ] * (h*SymBoltz.k0)),
    #"θmν" => (sol1["pt"]["theta_ncdm[0]"], sol2[1, M.h.θ] * (h*SymBoltz.k0)), # TODO: correct???
    "σγ" => (sol1["pt"]["shear_g"], sol2[1, M.γ.F[2] / 2]), # TODO: factor 2 difference
    "σν" => (sol1["pt"]["shear_ur"], sol2[1, M.ν.F[2] / 2]), # TODO: factor 2 difference
    "P0" => (sol1["pt"]["pol0_g"], sol2[1, M.γ.G[0]]),
    "P1" => (sol1["pt"]["pol1_g"], sol2[1, M.γ.G[1]]),
    "P2" => (sol1["pt"]["pol2_g"], sol2[1, M.γ.G[2]]),
)

# matter power spectrum
ks = sol1["P"]["k(h/Mpc)"] * h # 1/Mpc
Ps_class = sol1["P"]["P(Mpc/h)^3"] / h^3
Ps = power_spectrum(M, pars, ks / u"Mpc"; verbose=true) / u"Mpc^3"

# CMB power spectrum
ls = sol1["Cl"]["l"]
ls = Int.(ls[begin:10:end])
Dls_class = sol1["Cl"]["TT"]
Dls_class = Dls_class[begin:10:end]
Cls = Cl(M, pars, ls; verbose=true)
Dls = SymBoltz.Dl(Cls, ls)

sols = merge(sols, Dict(
    "k" => (ks, ks),
    "P" => (Ps_class, Ps),
    "l" => (ls, ls),
    "Dl" => (Dls_class, Dls)
))

function plot_compare(xlabel, ylabels; lgx=false, lgy=false, alpha=1.0)
    if !(ylabels isa AbstractArray)
        ylabels = [ylabels]
    end

    xplot(x) = lgx ? log10.(abs.(x)) : x
    yplot(y) = lgy ? log10.(abs.(y)) : y
    xlab(x) = lgx ? "lg(|$x|)" : x
    ylab(y) = lgy ? "lg(|$y|)" : y

    # TODO: relative or absolute comparison (of quantities close to 0)
    p = plot(; layout=grid(2, 1, heights=(2/3, 1/3)), size = (800, 700))
    title = join([(@sprintf "%s = %.3f" s Dict(pars)[s]) for s in keys(Dict(pars))], ", ") * (@sprintf ", k = %.2e / Mpc" k / u"1/Mpc")
    plot!(p[1]; title, titlefontsize = 8, ylabel = join(ylab.(ylabels), ", "))
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
        x = x2[begin+1:end-1] # compare ratios at SymBoltz' times # TODO: why must I exclude first and last points to avoid using extrapolate=true in splines?

        y1, y2 = sols[ylabel]
        y1 = y1[i1s]
        y2 = y2[i2s]
        plot!(p[1], xplot(x1), yplot(y1); color = :grey, linewidth = 2, alpha, linestyle = :solid, label = i == 1 ? "CLASS" : nothing, xlabel = xlab(xlabel), legend_position = :topright)
        plot!(p[1], xplot(x2), yplot(y2); color = :black, linewidth = 2, alpha, linestyle = :dash,  label = i == 1 ? "SymBoltz" : nothing)

        # TODO: use built-in CosmoloySolution interpolation
        y1 = LinearInterpolation(y1, x1; extrapolate=true).(x)
        y2 = LinearInterpolation(y2, x2; extrapolate=true).(x)
        r = @. abs(y2-y1) / max(abs(y1), abs(y2))
        plot!(p[end], xplot(x), r * 100; linewidth = 2, yminorticks = 10, yminorgrid = true, color = :black, ylims=(0, 10))
    end
    hline!(p[end], [0.0]; color = :black, linewidth = 2, linestyle = :dash, ylabel = "|y₂-y₁| / max(|y₁|, |y₂|) / %", z_order = 1)
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
plot_compare("a_bg", ["ργ", "ρν", "ρb", "ρc", "ρΛ"]; lgx=true, lgy=true) # hide
```

### Thermodynamics
```@example class
plot_compare("a_th", "τ̇"; lgx=true, lgy=true) # hide
```
```@example class
plot_compare("a_th", "Xe"; lgx=true, lgy=false) # hide
```
```@example class
plot_compare("a_th", "Tb"; lgx=true, lgy=true) # hide
```
```@example class
plot_compare("a_th", "csb²"; lgx=true, lgy=true) # hide
# TODO: Ṫ # hide
```

### Perturbations

```@example class
plot_compare("a_pt", ["Ψ", "Φ"]; lgx=true) # hide
```
```@example class
plot_compare("a_pt", ["δb", "δc", "δγ", "δν"]; lgx=true, lgy=true) # hide
```
```@example class
plot_compare("a_pt", ["θb", "θc", "θγ", "θν"]; lgx=true, lgy=true) # hide
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
plot_compare("l", "Dl") # TODO: fix # hide
```
