using RecipesBase

function displayname(var)
    name = string(var)
    name = replace(name, Symbolics.NAMESPACE_SEPARATOR => '.')
    return name
end

@recipe function plot(sol::CosmologySolution, x, y; Nextra = 0)
    if !(y isa AbstractArray)
        y = [y]
    end
    τs = timeseries(sol; Nextra)
    xs = sol(x, τs)
    ys = sol(y, τs)
    xlabel --> (x isa AbstractArray ? "" : displayname(x))
    ylabel --> (y isa AbstractArray ? "" : displayname(y))
    line_z = get(plotattributes, :line_z, nothing)
    if line_z isa Num
        line_z := sol(line_z, τs)
        colorbar_title --> displayname(line_z)
    end
    label --> displayname.(y')
    return xs, permutedims(ys)
end

@recipe function plot(sol::CosmologySolution, x, y, k; Nextra = 0, klabel = true)
    if !(y isa AbstractArray)
        y = [y]
    end

    xlabel --> (x isa AbstractArray ? "" : displayname(x))
    ylabel --> (y isa AbstractArray ? "" : displayname(y))

    for iv in eachindex(y)
        linestyle = [:solid :dash :dot :dashdot :dashdotdot][iv]
        for ik in eachindex(k)
            τs = timeseries(sol, k[ik]; Nextra)
            color = ik
            xs = sol(x, τs, k[ik])
            ys = sol(y[iv], τs, k[ik])
            @series begin
                linestyle --> linestyle
                color --> color
                label := ""
                xs, ys
            end

            # label wavenumber with dummy plot
            if iv == 1 && klabel
                @series begin
                    linestyle := :solid
                    color := color
                    label := "k = $(k[ik])"
                    [NaN], [NaN]
                end
            end
        end

        # label variable with dummy plot
        @series begin
            linestyle := linestyle
            color := :black
            label := displayname(y[iv])
            [NaN], [NaN]
        end
    end
end

# plot Systems as a hierarchical tree
using GraphRecipes
@recipe function plot(::Type{T}, sys::T) where {T <: System}
    nodeshape --> :rect
    nodesize --> 0.12
    fontsize --> 15
    TreePlot(sys)
end
