include("../Symboltz.jl")
import .Symboltz
using ModelingToolkit
using DelimitedFiles
using DataInterpolations
using Plots; Plots.default(label=nothing)

@kwdef struct Parameters
    Ωr0 = 5.5e-5
    Ωc0 = 0.267
    Ωb0 = 0.03
    h = 0.67
    As = 2.1e-9
    Yp = 0.01 # 245 # TODO: more
end
lmax = 6

par = Parameters()

function run_class(in::Dict{String, Any}; exec="class/class_public-3.2.3/class", inpath="class/input.ini", outpath="class/output/")
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

function output_class(par::Parameters, k::Real; kwargs...)
    in = Dict(
        "write_background" => "yes",
        "write_thermodynamics" => "yes",
        "k_output_values" => k,
        "output" => "mPk", # need one to evolve perturbations
        "ic" => "ad",
        "modes" => "s",
        "gauge" => "newtonian",
        "h" => par.h,
        "Omega_g" => par.Ωr0,
        "Omega_b" => par.Ωb0,
        "Omega_cdm" => par.Ωm0 - par.Ωb0,
        "Omega_dcdmdr" => 0.0,
        "Omega_k" => 0.0,
        "Omega_fld" => 0.0,
        "Omega_scf" => 0.0,
        "N_ur" => 0.0,
        "N_ncdm" => 0.0,
        "YHe" => par.Yp, # TODO: disable recombination and reionization?
        "recombination" => "hyrec", # or HyREC
        "l_max_g" => lmax,
        "l_max_pol_g" => lmax,
        "l_max_ur" => lmax,
        "l_max_ncdm" => lmax,
    )

    run_class(in; kwargs...)
    output = Dict()
    for (name, file, skipstart) in [("bg", "class/output/_background.dat", 3), ("th", "class/output/_thermodynamics.dat", 10), ("pt", "class/output/_perturbations_k0_s.dat", 1)]
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

kMpc = 1e-1 # 1/Mpc # disagreement on smaller scales
sol1 = output_class(par, kMpc)

k = kMpc ./ (par.h * Symboltz.k0) # h/Mpc -> code units
@named bg = Symboltz.background_ΛCDM()
@named th = Symboltz.thermodynamics_ΛCDM(bg)
@named pt = Symboltz.perturbations_ΛCDM(th, lmax)
sol2_th = Symboltz.solve(th, par.Ωr0, par.Ωm0, par.Ωb0, par.h, par.Yp; reltol = 1e-10)
sol2_pt = Symboltz.solve(pt, [k], par.Ωr0, par.Ωm0, par.Ωb0, par.h, par.Yp; reltol = 1e-10)[1]

# map results from both codes to common convention
results = Dict(
    # background
    "lg(a_bg)" => (log10.(1 ./ (sol1["bg"]["z"] .+ 1)), log10.(sol2_pt[bg.sys.g.a])),
    "E" => (sol1["bg"]["H[1/Mpc]"] ./ sol1["bg"]["H[1/Mpc]"][end], sol2_pt[bg.sys.g.E]),

    # thermodynamics
    "lg(a_th)" => (log10.(reverse(sol1["th"]["scalefactora"])), log10.(sol2_th[bg.sys.g.a])),
    "lg(τ′)" => (log10.(reverse(sol1["th"]["kappa'[Mpc^-1]"])), log10.(.- sol2_th[th.ssys.dτ] * (Symboltz.k0 * par.h))),
    "lg(|csb²|)" => (log10.(abs.(reverse(sol1["th"]["c_b^2"]))), log10.(abs.(sol2_th[th.sys.cs²]))),
    "Xe" => (reverse(sol1["th"]["x_e"]), sol2_th[th.sys.Xe]),
    "Tb" => (reverse(sol1["th"]["Tb[K]"]), sol2_th[th.sys.Tb]),

    # perturbations
    "lg(a_pt)" => (log10.(sol1["pt"]["a"]), log10.(sol2_pt[bg.sys.g.a])),
    "a_pt" => (sol1["pt"]["a"], sol2_pt[bg.sys.g.a]),
    "Φ" => (sol1["pt"]["phi"], sol2_pt[pt.ssys.g1.Φ]), # TODO: same?
    "Ψ" => (sol1["pt"]["psi"], -sol2_pt[pt.ssys.g1.Ψ]), # TODO: same?
    "δb" => (sol1["pt"]["delta_b"], -sol2_pt[pt.ssys.bar.δ]), # TODO: sign?
    "δc" => (sol1["pt"]["delta_cdm"], -sol2_pt[pt.ssys.cdm.δ]), # TODO: sign?
    "δγ" => (sol1["pt"]["delta_g"], -sol2_pt[pt.ssys.ph.δ]),
    "θb" => (sol1["pt"]["theta_b"], -sol2_pt[pt.ssys.bar.u] * kMpc),
    "θc" => (sol1["pt"]["theta_cdm"], -sol2_pt[pt.ssys.cdm.u] * kMpc),
    "θγ" => (sol1["pt"]["theta_g"], -sol2_pt[pt.ssys.ph.Θ[1]] * 3 * kMpc),
    "Π" => (sol1["pt"]["shear_g"], sol2_pt[pt.ssys.ph.Θ[2]] * -2),
    "P0" => (sol1["pt"]["pol0_g"], sol2_pt[pt.ssys.ph.ΘP0] * -4), # TODO: is -4 correct ???
    "P1" => (sol1["pt"]["pol1_g"], sol2_pt[pt.ssys.ph.ΘP[1]] * -4), # TODO: is -4 correct ???
    "P2" => (sol1["pt"]["pol2_g"], sol2_pt[pt.ssys.ph.ΘP[2]] * -4), # TODO: is -4 correct ???
)

# TODO: relative or absolute comparison (of quantities close to 0)
xlabel, ylabels = "lg(a_pt)", ["Φ", "Ψ"]
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
p = plot(; layout = (length(ylabels)+1, 1), size = (700, 800))
title = join(["$s = $(getfield(par, s))" for s in fieldnames(Parameters)], ", ") * ", k = $(kMpc) / Mpc"
plot!(p[1]; title, titlefontsize = 9)
for (i, ylabel) in enumerate(ylabels)
    y1, y2 = results[ylabel]
    y1 = y1[i1s]
    y2 = y2[i2s]
    color = i
    plot!(p[i], x1, y1; color, linestyle = :dash, label = "CLASS", ylabel)
    plot!(p[i], x2, y2; color, linestyle = :solid, label = "Symboltz")
    y1 = CubicSpline(y1, x1).(x)
    y2 = CubicSpline(y2, x2).(x)
    r = y2 ./ y1
    plot!(p[end], x, r; yminorticks = 10, yminorgrid = true, color)
end
hline!(p[end], [1.0]; color = :black, linestyle = :dash, xlabel, ylabel = "Symboltz / CLASS", z_order = 1)