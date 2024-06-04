import .Symboltz
using ModelingToolkit
using DelimitedFiles
using DataInterpolations
using DifferentialEquations
using Plots; Plots.default(label=nothing)

@kwdef struct Parameters
    Ωr0 = 5.5e-5
    Ωm0 = 0.317
    Ωb0 = 0.03
    h = 0.67
    As = 2.1e-9
    Yp = 0.01 # 245 # TODO: more
    lmax = 6
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
        "l_max_g" => par.lmax,
        "l_max_pol_g" => par.lmax,
        "l_max_ur" => par.lmax,
        "l_max_ncdm" => par.lmax,
        "radiation_streaming_approximation" => 0,
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

kMpc = 1e-1 # 1/Mpc # disagreement on smaller scales
sol1 = perts_class(par, kMpc)

k = kMpc ./ (par.h * Symboltz.k0) # h/Mpc -> code units
@named bg = Symboltz.background_ΛCDM()
@named th = Symboltz.thermodynamics_ΛCDM(bg)
@named pt = Symboltz.perturbations_ΛCDM(th, par.lmax)
sol2 = Symboltz.solve(pt, [k], par.Ωr0, par.Ωm0, par.Ωb0, par.h, par.Yp; solver = KenCarp47(), reltol = 1e-10)[1]

# map results from both codes to common convention
results = Dict(
    "η" => (sol1["tau[Mpc]"], sol2[Symboltz.η] / (par.h*Symboltz.k0)),
    "a" => (sol1["a"], sol2[bg.sys.g.a]),
    "Φ" => (sol1["phi"], sol2[pt.ssys.gpt.Φ]), # TODO: same?
    "Ψ" => (sol1["psi"], -sol2[pt.ssys.gpt.Ψ]), # TODO: same?
    "δb" => (sol1["delta_b"], -sol2[pt.ssys.bar.δ]), # TODO: sign?
    "δc" => (sol1["delta_cdm"], -sol2[pt.ssys.cdm.δ]), # TODO: sign?
    "δγ" => (sol1["delta_g"], -sol2[pt.ssys.ph.δ]),
    "θb" => (sol1["theta_b"], -sol2[pt.ssys.bar.u] * kMpc),
    "θc" => (sol1["theta_cdm"], -sol2[pt.ssys.cdm.u] * kMpc),
    "θγ" => (sol1["theta_g"], -sol2[pt.ssys.ph.Θ[1]] * 3 * kMpc),
    "Π" => (sol1["shear_g"], sol2[pt.ssys.ph.Θ[2]] * -2),
    "P0" => (sol1["pol0_g"], sol2[pt.ssys.ph.ΘP0] * -4), # TODO: is -4 correct ???
    "P1" => (sol1["pol1_g"], sol2[pt.ssys.ph.ΘP[1]] * -4), # TODO: is -4 correct ???
    "P2" => (sol1["pol2_g"], sol2[pt.ssys.ph.ΘP[2]] * -4), # TODO: is -4 correct ???

    "lg(η)" => (log.(sol1["tau[Mpc]"]), log.(sol2[Symboltz.η] / (par.h*Symboltz.k0))),
    "lg(a)" => (log10.(sol1["a"]), log10.(sol2[bg.sys.g.a])),
)

xlabel = "lg(η)"
ylabel = "a"
x1, x2 = results[xlabel]
y1, y2 = results[ylabel]
x = x1
r = CubicSpline(y1, x1; extrapolate=true).(x) ./ CubicSpline(y2, x2; extrapolate=true).(x) # ratio at common x # TODO: verify extrapolate is ok
xlims = extrema(x)

p = plot(layout = (2, 1), size = (600, 800))
plot!(p[1], x1, y1; label = "CLASS", ylabel, title = "k = $(k*Symboltz.k0) h/Mpc")
plot!(p[1], x2, y2; label = "Symboltz", xlims)
plot!(p[2], x, r; ylims=(0.9, 1.1), xlims, xlabel, ylabel = "Symboltz / CLASS")