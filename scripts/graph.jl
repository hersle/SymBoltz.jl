using SymBoltz, ModelingToolkit, Plots, GraphRecipes
Plots.default(label=nothing)

M = ΛCDM(h = nothing, ν = nothing, K = nothing, lmax = 3)
Ms = structural_simplify(M)
tstate = ModelingToolkit.get_tearing_state(Ms)
imat = ModelingToolkit.BipartiteGraphs.incidence_matrix(tstate.structure.graph)
function plot_dependency_graph(sys)
    tstate = ModelingToolkit.get_tearing_state(sys) # works on simplified systems # TODO: what to do on unsimplified systems?

    names = string.(tstate.fullvars)
    names = replace.(names, "Differential(t)" => "D", "(t)" => "") # TODO: consider ()′ with something like https://stackoverflow.com/a/35271017

    # color each variable based on its subsystem
    names_cleaned = replace.(names, (["D", "(", ")"] .=> "")...) # remove any derivative signs
    subsystems = ['₊' in name ? name[begin:prevind(name, findlast('₊', name))] : "" for name in names_cleaned]
    usubsystems = unique(subsystems)
    groups = [findfirst(usub -> usub == sub, usubsystems) for sub in subsystems]

    plot(
        tstate.structure.graph, names = names, method = :stress,
        fontsize = 3, nodesize = 0.3, nodeshape = :circle, nodecolor = groups,
        self_edge_size = 0 # hide self-to-self nodes
    ) # see https://docs.juliaplots.org/stable/generated/graph_attributes/
end

plot_dependency_graph(Ms)