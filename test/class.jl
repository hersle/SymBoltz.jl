import .Symboltz
using ModelingToolkit
using DelimitedFiles
using DataInterpolations
using Plots; Plots.default(label=nothing)

@kwdef struct Parameters
    Ωr0 = 5.5e-5
    Ωm0 = 0.317
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
sol2 = Symboltz.solve(pt, [k], par.Ωr0, par.Ωm0, par.Ωb0, par.h, par.Yp; reltol = 1e-10)[1]

# map results from both codes to common convention
results = Dict(
    # background
    "lg(a_bg)" => (log10.(1 ./ (sol1["bg"]["z"] .+ 1)), log10.(sol2[bg.sys.g.a])),
    "E" => (sol1["bg"]["H[1/Mpc]"] ./ sol1["bg"]["H[1/Mpc]"][end], sol2[bg.sys.g.E]),

    # thermodynamics
    "lg(a_th)" => (log10.(reverse(sol1["th"]["scalefactora"])), log10.(sol2[bg.sys.g.a])),
    "lg(τ′)" => (log10.(reverse(sol1["th"]["kappa'[Mpc^-1]"])), log10.(.- sol2[pt.ssys.dτ] * (Symboltz.k0 * par.h))),

    # perturbations
    "lg(a_pt)" => (log10.(sol1["pt"]["a"]), log10.(sol2[bg.sys.g.a])),
    "a_pt" => (sol1["pt"]["a"], sol2[bg.sys.g.a]),
    "Φ" => (sol1["pt"]["phi"], sol2[pt.ssys.gpt.Φ]), # TODO: same?
    "Ψ" => (sol1["pt"]["psi"], -sol2[pt.ssys.gpt.Ψ]), # TODO: same?
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

xlabel = "lg(a_pt)"
ylabel = "Φ"
x1, x2 = results[xlabel]
y1, y2 = results[ylabel]
x1min, x1max = extrema(x1)
x2min, x2max = extrema(x2)
x = range(max(x1min, x2min), min(x1max, x2max), length = 1000)
y1min, y1max = extrema(y1)
y2min, y2max = extrema(y2)
ymin = min(y1min, y2min)
ymax = max(y1max, y2max)
ylims = (ymin, ymax)

r = CubicSpline(y1, x1; extrapolate=true).(x) ./ CubicSpline(y2, x2; extrapolate=true).(x) # ratio at common x # TODO: verify extrapolate is ok
xlims = extrema(x)

p = plot(layout = (2, 1), size = (700, 800))
plot!(p[1], x1, y1; label = "CLASS", ylabel, title = "k = $(kMpc) / Mpc")
plot!(p[1], x2, y2; label = "Symboltz", xlims, ylims)
plot!(p[2], x, r; ylims = (0.8, 1.2), yticks = [0.8, 0.9, 1.0, 1.1, 1.2], yminorticks = 0.8:0.01:1.2, yminorgrid = true, xlims, xlabel, ylabel = "Symboltz / CLASS")