module SymBoltzMakieExt

using Makie
import SymBoltz: plot_interactive # import (not using) so the method below extends this, rather than shadowing it locally
using SymBoltz: CosmologyProblem, displayname, parameter_updater

"""
    plot_interactive(f::Function, prob::CosmologyProblem, pars::AbstractVector{<:Pair}; xlabel = "", ylabel = "", kwargs...)

Interactively plot data computed by `f(prob′)`, with one slider per parameter in `pars` (pairs of a parameter and its slider range).
The plot updates live as sliders are dragged: `prob` is remade into `prob′` with the current slider values and passed to `f`, which should return points to plot, e.g. `zip(xs, ys)`.
"""
function plot_interactive(f::Function, prob::CosmologyProblem, pars::AbstractVector{<:Pair}; xlabel = "", ylabel = "", kwargs...)
    fig = Figure(size = (800, 800), fontsize = 18)
    ax = Axis(fig[1, 1]; xlabel, ylabel)
    sg = SliderGrid(fig[2, 1], [(label = displayname(par), range = range, startvalue = prob.bg.ps[par]) for (par, range) in pars]...)

    probgen = parameter_updater(prob, first.(pars))
    θ = Observable([prob.bg.ps[par] for (par, _) in pars])
    for (i, slider) in enumerate(sg.sliders)
        on(slider.value) do val
            θ[][i] = val
            notify(θ)
        end
    end

    lines!(ax, @lift(f(probgen($θ))); kwargs...)
    return fig
end

end
