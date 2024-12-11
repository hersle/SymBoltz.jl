using RecipesBase
using Roots

function get_ts(sol::CosmologySolution, Nextra)t
    ts = sol[t]
    return extend_ts(ts, Nextra)
end
function get_ts(sol::CosmologySolution, k, Nextra)
    k = k_dimensionless.(k, sol.bg.ps[:h])
    i1, i2 = get_neighboring_wavenumber_indices(sol, k)
    i1 = max(i1, 1)
    t1, t2 = sol[i1, t], sol[i2, t]
    ts = sort!(unique!([t1; t2]))
    return extend_ts(ts, Nextra)
end
function extend_ts(ts, Nextra)
    return exp.(extend_array(log.(ts), Nextra))
end

# Get the time when some variable equals some value
# TODO: use for something
function get_ts_when(sol::CosmologySolution, var, val)
    allequal(sign.(diff(sol[var]))) || error("$var is not monotonically increasing/decreasing")
    f(t) = sol(t, var) - val # var(t) == val when f(t) == 0
    t0 = sum(sol[t][[begin,end]]) / 2 # initial guess at middle time
    return find_zero(f, t0)
end

@recipe function plot(sol::CosmologySolution, x, y; Nextra = 0)
    ts = get_ts(sol, Nextra)
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
            ts = get_ts(sol, k[ik], Nextra)
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

@recipe function plot(::Type{T}, M::T) where {T <: CosmologyModel}
    return (M.sys,) # call ODESystem recipe
end
