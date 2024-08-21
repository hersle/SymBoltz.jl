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

    for iv in eachindex(y)
        linestyle = [:solid :dash :dot :dashdot :dashdotdot][iv]
        for ik in eachindex(k)
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
            if iv == 1
                @series begin
                    linestyle := :solid
                    color := color
                    label := "k = $(round(k[ik]; digits=3)) H₀/c"
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

# plot ODESystems as a hierarchical tree # TODO: contribute back to ModelingToolkit.jl?
using AbstractTrees
using SymbolicUtils
function AbstractTrees.children(sys::ODESystem; variables = false)
    syss = collect(sys.systems)
    if variables && isempty(syss)
        return union(unknowns(sys), observed(sys), parameters(sys))
    else
        return syss
    end
end
function AbstractTrees.children(var::SymbolicUtils.BasicSymbolic)
    return []
end
AbstractTrees.printnode(io::IO, sys::ODESystem) = print(io, sys.name)
AbstractTrees.printnode(io::IO, var::SymbolicUtils.BasicSymbolic) = print(io, var)

using GraphRecipes
@recipe function plot(::Type{T}, sys::T) where {T <: ODESystem}
    nodeshape --> :rect
    nodesize --> 0.12
    fontsize --> 15
    TreePlot(sys)
end

@recipe function plot(::Type{T}, M::T) where {T <: CosmologyModel}
    return (M.sys,) # call ODESystem recipe
end
