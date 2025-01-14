using RecipesBase
using Roots

@recipe function plot(sol::CosmologySolution, x, y; Nextra = 0)
    ts = timeseries(sol; Nextra)
    xs = sol(ts, x)
    ys = sol(ts, y)
    xlabel --> (x isa AbstractArray ? "" : x)
    ylabel --> (y isa AbstractArray ? "" : y)
    line_z = get(plotattributes, :line_z, nothing)
    if line_z isa Num # TODO: add to perturbations?
        line_z := sol(ts, line_z)
        colorbar_title --> line_z
    end
    label --> y'
    return xs, ys
end

@recipe function plot(sol::CosmologySolution, k, x, y; Nextra = 0, klabel = true)
    xlabel --> (x isa AbstractArray ? "" : x)
    ylabel --> (y isa AbstractArray ? "" : y)

    for iv in eachindex(y)
        linestyle = [:solid :dash :dot :dashdot :dashdotdot][iv]
        for ik in eachindex(k)
            ts = timeseries(sol, k[ik]; Nextra)
            color = ik
            xs = sol(k[ik], ts, x)
            ys = sol(k[ik], ts, y[iv])
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
            label := y[iv]
            [NaN], [NaN]
        end
    end
end

# plot ODESystems as a hierarchical tree
using GraphRecipes
@recipe function plot(::Type{T}, sys::T) where {T <: ODESystem}
    nodeshape --> :rect
    nodesize --> 0.12
    fontsize --> 15
    TreePlot(sys)
end
