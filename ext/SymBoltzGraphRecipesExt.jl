module SymBoltzGraphRecipesExt

import GraphRecipes
using SymBoltz: System

# plot Systems as a hierarchical tree
GraphRecipes.@recipe function plot(::Type{T}, sys::T) where {T <: System}
    nodeshape --> :rect
    nodesize --> 0.12
    fontsize --> 15
    GraphRecipes.TreePlot(sys)
end

end
