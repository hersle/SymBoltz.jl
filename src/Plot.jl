using RecipesBase

@recipe function plot(sol::CosmologySolution, x, y; N = 500)
    ts = exp.(range(log.(extrema(sol[t]))..., length=N))
    xs = sol(ts, x)
    ys = sol(ts, y)
    xlabel --> (x isa AbstractArray ? "" : x)
    ylabel --> (y isa AbstractArray ? "" : y)
    label --> y'
    return xs, ys
end

@recipe function plot(sol::CosmologySolution, k, x, y; N = 500)
    ts = exp.(range(log.(extrema(sol[t]))..., length=N))
    xs = sol(k, ts, x)
    ys = sol(k, ts, y)

    if k isa AbstractArray
        # stack multiple k into several timeseries y1(k1) y2(k1) y3(k1) ... y1(k2) y2(k2) y3(k2) ...
        xs = hcat((xs[i, :, :] for i in axes(xs, 1))...)
        ys = hcat((ys[i, :, :] for i in axes(ys, 1))...)
        color --> vec([j for i in 1:length(y), j in 1:length(k)])' # per-k color
        linestyle --> [:solid :dash :dot :dashdot :dashdotdot][vec([i for i in 1:length(y), j in 1:length(k)])'] # per-variable linestyle
        ks = k[vec([j for i in 1:length(y), j in 1:length(k)])]
        label --> permutedims(repeat(string.(y), length(k)) .* (", k = $k H₀/c" for k in ks))
    else
        label --> y'
    end

    xlabel --> (x isa AbstractArray ? "" : x)
    ylabel --> (y isa AbstractArray ? "" : y)
    label --> y'
    return xs, ys
end