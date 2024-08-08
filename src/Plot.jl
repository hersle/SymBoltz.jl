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

@recipe function plot(sol::CosmologySolution, k::AbstractArray, x::AbstractArray, y::AbstractArray; N = 500)
    ts = exp.(range(log.(extrema(sol[t]))..., length=N))
    xs = sol(k, ts, x)
    ys = sol(k, ts, y)

    for iv in eachindex(y)
        linestyle = [:solid :dash :dot :dashdot :dashdotdot][iv]
        for ik in eachindex(k)
            color = ik

            @series begin
                color --> color
                linestyle --> linestyle
                label --> ""
                xs[ik, :, iv], ys[ik, :, iv]
            end

            # dummy plot for wavenumber label
            if iv == 1
                @series begin
                    color := color
                    linestyle := :solid
                    label := "k = $(k[ik]) H₀/c"
                    [NaN], [NaN]
                end
            end
        end

        # dummy plot for variable label
        @series begin
            color := :black
            linestyle := linestyle
            label := y[iv]
            [NaN], [NaN]
        end
    end
end

@recipe function plot(sol::CosmologySolution, k, x::Number, y; N = 500)
    sol, k, fill(x, length(y)), y
end