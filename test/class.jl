using ModelingToolkit
using DelimitedFiles
using DataInterpolations
using Plots; Plots.default(label=nothing)

@kwdef struct Parameters
    Ωγ0 = 5.5e-5
    Ων0 = 3.046 * 7/8 * (4/11)^(4/3) * 5.5e-5 # 0.0
    Ωc0 = 0.267
    Ωb0 = 0.03
    h = 0.67
    As = 2.1e-9
    Yp = 0.245
end
lmax = 6

par = Parameters()

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

function output_class(par::Parameters, k::Real; exec="class", inpath="/tmp/symboltz_class/input.ini", outpath="/tmp/symboltz_class/output/")
    in = Dict(
        "write_background" => "yes",
        "write_thermodynamics" => "yes",
        "k_output_values" => k,
        "output" => "mPk", # need one to evolve perturbations
        "ic" => "ad",
        "modes" => "s",
        "gauge" => "newtonian",
        "h" => par.h,
        "Omega_g" => par.Ωγ0,
        "Omega_b" => par.Ωb0,
        "Omega_cdm" => par.Ωc0,
        "Omega_ur" => par.Ων0, # massless neutrinos
        "Omega_dcdmdr" => 0.0,
        "Omega_k" => 0.0,
        "Omega_fld" => 0.0,
        "Omega_scf" => 0.0,
        "N_ncdm" => 0.0,
        "YHe" => par.Yp, # TODO: disable recombination and reionization?
        "recombination" => "recfast", # or HyREC
        "recfast_Hswitch" => 1,
        "recfast_Heswitch" => 6,
        "reio_parametrization" => "reio_none",
        "l_max_g" => lmax,
        "l_max_pol_g" => lmax,
        "l_max_ur" => lmax,
        "l_max_ncdm" => lmax,
    )

    run_class(in, exec, inpath, outpath)
    output = Dict()
    for (name, filename, skipstart) in [("bg", "_background.dat", 3), ("th", "_thermodynamics.dat", 10), ("pt", "_perturbations_k0_s.dat", 1)]
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

kMpc = 1e0 # 1/Mpc # disagreement on smaller scales
sol1 = output_class(par, kMpc)

k = kMpc ./ (par.h * Symboltz.k0) # h/Mpc -> code units
@named bg = Symboltz.background_ΛCDM()
@named th = Symboltz.thermodynamics_ΛCDM(bg)
@named pt = Symboltz.perturbations_ΛCDM(th, lmax)
sol2 = Symboltz.solve(pt, [k], par.Ωγ0, par.Ων0, par.Ωc0, par.Ωb0, par.h, par.Yp)[1]

# map results from both codes to common convention
results = Dict(
    # background
    "lg(a_bg)" => (log10.(1 ./ (sol1["bg"]["z"] .+ 1)), log10.(sol2[bg.sys.g.a])),
    "E" => (sol1["bg"]["H[1/Mpc]"] ./ sol1["bg"]["H[1/Mpc]"][end], sol2[bg.sys.g.E]),

    # thermodynamics
    "lg(a_th)" => (log10.(reverse(sol1["th"]["scalefactora"])), log10.(sol2[bg.sys.g.a])),
    "lg(τ′)" => (log10.(reverse(sol1["th"]["kappa'[Mpc^-1]"])), log10.(.- sol2[th.sys.dτ] * (Symboltz.k0 * par.h))),
    "lg(csb²)" => (log10.(reverse(sol1["th"]["c_b^2"])), log10.(sol2[th.sys.cs²])),
    "csb²" => (reverse(sol1["th"]["c_b^2"]), sol2[th.sys.cs²]),
    "Xe" => (reverse(sol1["th"]["x_e"]), sol2[th.sys.Xe]),
    "Tb" => (reverse(sol1["th"]["Tb[K]"]), sol2[th.sys.Tb]),

    # perturbations
    "lg(a_pt)" => (log10.(sol1["pt"]["a"]), log10.(sol2[bg.sys.g.a])),
    "a_pt" => (sol1["pt"]["a"], sol2[bg.sys.g.a]),
    "Φ" => (sol1["pt"]["phi"], sol2[pt.ssys.g1.Φ]), # TODO: same?
    "Ψ" => (sol1["pt"]["psi"], -sol2[pt.ssys.g1.Ψ]), # TODO: same?
    "δb" => (sol1["pt"]["delta_b"], -sol2[pt.ssys.bar.δ]), # TODO: sign?
    "δc" => (sol1["pt"]["delta_cdm"], -sol2[pt.ssys.cdm.δ]), # TODO: sign?
    "δγ" => (sol1["pt"]["delta_g"], -sol2[pt.ssys.ph.δ]),
    "θb" => (sol1["pt"]["theta_b"], -sol2[pt.ssys.bar.u] * kMpc),
    "θc" => (sol1["pt"]["theta_cdm"], -sol2[pt.ssys.cdm.u] * kMpc),
    "θγ" => (sol1["pt"]["theta_g"], -sol2[pt.ssys.ph.Θ[1]] * 3 * kMpc),
    "Π" => (sol1["pt"]["shear_g"], sol2[pt.ssys.ph.Θ[2]] * -2),
    "P0" => (sol1["pt"]["pol0_g"], sol2[pt.ssys.ph.ΘP0] * -4), # TODO: is -4 correct ???
    "P1" => (sol1["pt"]["pol1_g"], sol2[pt.ssys.ph.ΘP[1]] * -4), # TODO: is -4 correct ???
    "P2" => (sol1["pt"]["pol2_g"], sol2[pt.ssys.ph.ΘP[2]] * -4), # TODO: is -4 correct ???
)

# TODO: relative or absolute comparison (of quantities close to 0)
xlabels, ylabels = ["lg(a_th)", "lg(a_pt)"], ["Xe", "Φ"]
p = plot(; layout = (length(ylabels)+1, 1), size = (700, 800))
title = join(["$s = $(getfield(par, s))" for s in fieldnames(Parameters)], ", ") * ", k = $(kMpc) / Mpc"
plot!(p[1]; title, titlefontsize = 9)
for (i, (xlabel, ylabel)) in enumerate(zip(xlabels, ylabels))
    x1, x2 = results[xlabel]
    x1min, x1max = extrema(x1)
    x2min, x2max = extrema(x2)
    xmin = max(x1min, x2min)
    xmax = min(x1max, x2max)
    i1s = xmin .≤ x1 .≤ xmax .&& [BitVector([true]); x1[begin+1:end] .- x1[begin:end-1] .> 1e-5] # select x values in common range and exclude points too close in time (probably due to CLASS switching between approximation schemes)
    i2s = xmin .≤ x2 .≤ xmax # select x values in common range
    x1 = x1[i1s]
    x2 = x2[i2s]
    x = x2[begin+1:end-1] # compare ratios at Symboltz' times # TODO: why must I exclude first and last points to avoid using extrapolate=true in splines?

    y1, y2 = results[ylabel]
    y1 = y1[i1s]
    y2 = y2[i2s]
    color = i
    plot!(p[i], x1, y1; color, linestyle = :dash, label = "CLASS (y₁)", xlabel, ylabel)
    plot!(p[i], x2, y2; color, linestyle = :solid, label = "Symboltz (y₂)")
    y1 = CubicSpline(y1, x1; extrapolate=true).(x)
    y2 = CubicSpline(y2, x2; extrapolate=true).(x)
    r = @. abs(y2-y1) / max(abs(y1), abs(y2))
    plot!(p[end], x, r * 100; yminorticks = 10, yminorgrid = true, color)
end
hline!(p[end], [0.0]; color = :black, linestyle = :dash, ylabel = "|y₂-y₁| / max(|y₁|, |y₂|) / %", z_order = 1)