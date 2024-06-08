include("../Symboltz.jl")
import .Symboltz
using ModelingToolkit
using Plots; Plots.default(label=nothing)
using GraphRecipes

@named bg = Symboltz.background_ΛCDM()
@named th = Symboltz.thermodynamics_ΛCDM(bg)
@named pt = Symboltz.perturbations_ΛCDM(th, 6)

tstate = ModelingToolkit.get_tearing_state(pt.ssys)
imat = ModelingToolkit.BipartiteGraphs.incidence_matrix(tstate.structure.graph)
function plot_dependency_graph(sys)
    tstate = ModelingToolkit.get_tearing_state(sys) # works on simplified systems # TODO: what to do on unsimplified systems?

    names = string.(tstate.fullvars)
    names = replace.(names, "Differential(η)" => "Dη", "(η)" => "") # TODO: consider ()′ with something like https://stackoverflow.com/a/35271017

    # color each variable based on its subsystem
    names_cleaned = replace.(names, (["Dη", "(", ")"] .=> "")...) # remove any derivative signs
    subsystems = ['₊' in name ? name[begin:findlast('₊', name)-1] : "" for name in names_cleaned]
    usubsystems = unique(subsystems)
    groups = [findfirst(usub -> usub == sub, usubsystems) for sub in subsystems]

    plot(
        tstate.structure.graph, names = names, method = :stress,
        fontsize = 5, nodeshape = :hexagon, nodecolor = groups,
        self_edge_size = 0 # hide self-to-self nodes
    ) # see https://docs.juliaplots.org/stable/generated/graph_attributes/
end

plot_dependency_graph(pt.ssys)