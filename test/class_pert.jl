import .Symboltz
using ModelingToolkit
using DelimitedFiles
using DataInterpolations
using Plots; Plots.default(label=nothing)

@kwdef struct Parameters
    Ωr0 = 5.5e-5
    Ωm0 = 0.317
    Ωb0 = 0.001 # 5 # TODO: more
    h = 0.67
    As = 2.1e-9
    Yp = 0.01 # 245 # TODO: more
end

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

function perts_class(par::Parameters, k::Real; kwargs...)
    in = Dict(
        "output" => "mPk", # need one to evolve perturbations
        "k_output_values" => k,
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
    )

    run_class(in; kwargs...)
    data, head = readdlm("class/output/_perturbations_k0_s.dat", skipstart=1, header=true)
    head = split(join(head, ""), ":")
    for (n, h) in enumerate(head[begin:end-1])
        head[n] = h[begin:end-length(string(n))]
    end
    head = string.(head[2:end]) # remove #
    data = data[:,begin:length(head)] # remove extra empty data rows because of CLASS' messed up names
    data = Matrix{Float64}(data)
    return Dict(head[i] => data[:,i] for i in 1:length(head))
end

k = 1e-4 # 1/Mpc # agree on large scales, disagree on small scales
sol1 = perts_class(par, k)

k = k ./ par.h / Symboltz.k0 # h/Mpc -> code units
@named bg = Symboltz.background_ΛCDM()
@named th = Symboltz.thermodynamics_ΛCDM(bg)
@named pt = Symboltz.perturbations_ΛCDM(th, 6)
sol2 = Symboltz.solve(pt, [k], par.Ωr0, par.Ωm0, par.Ωb0, par.h, par.Yp)[1]

results = Dict(
    "lg(a)" => (log10.(sol1["a"]), log10.(sol2[bg.sys.g.a])),
    "a" => (sol1["a"], sol2[bg.sys.g.a]),
    "Φ" => (sol1["phi"], sol2[pt.ssys.gpt.Φ]), # TODO: same?
    "Ψ" => (sol1["psi"], -sol2[pt.ssys.gpt.Ψ]), # TODO: same?
    "δb" => (sol1["delta_b"], -sol2[pt.ssys.bar.δ]), # TODO: sign?
    "δc" => (sol1["delta_cdm"], -sol2[pt.ssys.cdm.δ]), # TODO: sign?
    "δγ" => (sol1["delta_g"], -sol2[pt.ssys.ph.δ])
)

xlabel = "lg(a)"
ylabel = "Ψ"
x1, x2 = results[xlabel]
y1, y2 = results[ylabel]
x = x1
r = CubicSpline(y1, x1).(x) ./ CubicSpline(y2, x2).(x) # ratio at common x
xlims = extrema(x)

p = plot(layout = (2, 1), size = (600, 800))
plot!(p[1], x1, y1; label = "CLASS", ylabel, title = "k = $(k*Symboltz.k0) h/Mpc")
plot!(p[1], x2, y2; label = "Symboltz", xlims)
plot!(p[2], x, r; ylims=(0.9, 1.1), xlims, xlabel, ylabel = "Symboltz / CLASS")