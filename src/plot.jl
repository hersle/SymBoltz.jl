import RecipesBase

function displayname(var)
    name = string(var)
    name = replace(name, Symbolics.NAMESPACE_SEPARATOR => '.')
    return name
end

RecipesBase.@recipe function plot(sol::CosmologySolution, x, y; Nextra = 0)
    xlabel --> (x isa AbstractArray ? "" : displayname(x))
    ylabel --> (y isa AbstractArray ? "" : displayname(y))
    if !(y isa AbstractArray)
        y = [y]
    end
    τs = timeseries(sol; Nextra)
    xs = sol(x, τs)
    ys = sol(y, τs)
    line_z = get(plotattributes, :line_z, nothing)
    if line_z isa Num
        line_z := sol(line_z, τs)
        colorbar_title --> displayname(line_z)
    end
    label --> displayname.(y')
    return xs, permutedims(stack(ys))
end

RecipesBase.@recipe function plot(sol::CosmologySolution, x, y, k; Nextra = 0, transpose_cycle = false)
    xlabel --> (x isa AbstractArray ? "" : displayname(x))
    ylabel --> (y isa AbstractArray ? "" : displayname(y))
    if !(y isa AbstractArray)
        y = [y]
    end
    if !(k isa AbstractArray)
        k = [k]
    end
    linestyles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    for iv in eachindex(y)
        for ik in eachindex(k)
            τs = timeseries(sol, k[ik]; Nextra)
            xs = sol(x, τs, k[ik])
            ys = sol(y[iv], τs, k[ik])
            RecipesBase.@series begin
                linestyle --> (linestyles[mod1(transpose_cycle ? iv : ik, 5)])
                color --> (transpose_cycle ? ik : iv)
                label --> "$(displayname(y[iv])), k = $(k[ik])"
                xs, ys
            end
        end
    end
end

# TODO: Makie recipes (currently don't work well because of compatibility constraints; e.g. PairPlots requires Makie )
#=
import MakieCore
MakieCore.@recipe CosmologyPlot begin
    downcolor = :red
    upcolor = :green
end
function MakieCore.plot!(cp::CosmologyPlot{<:Tuple{AbstractVector{<:Real}, AbstractVector{<:Real}}})
    lines!(cp, [0.0, 1.0], [0.0, 1.0])
    return cp
end
=#

# plot Systems as a hierarchical tree
import GraphRecipes
GraphRecipes.@recipe function plot(::Type{T}, sys::T) where {T <: System}
    nodeshape --> :rect
    nodesize --> 0.12
    fontsize --> 15
    TreePlot(sys)
end
