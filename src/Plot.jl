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

    xlabel --> (x isa AbstractArray ? "" : x)
    ylabel --> (y isa AbstractArray ? "" : y)
    label --> y'
    return xs, ys
end