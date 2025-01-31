import SpecialFunctions: zeta as ζ
import Symbolics: taylor
import Base: identity
using QuadGK

∫(f, a, b) = quadgk(f, a, b)[1] # TODO: use integrator that supports dual numbers
∫(f, w) = sum(w .*  f) # ≈ ∫f(x)dx over weights found from QuadGK.gauss()

# callback for terminating an integrator when var == val0
# TODO: replace with symbolic continuous_events
function callback_terminator(sys, var, val0)
    varindex = ModelingToolkit.variable_index(sys, var)
    return ContinuousCallback((u, _, _) -> (val = u[varindex]; val - val0), terminate!)
end

function transform(f::Function, sys::ODESystem; fullname=string(sys.name))
    subs = [transform(f, sub; fullname = (ModelingToolkit.iscomplete(sys) ? "" : fullname * "₊") * string(sub.name)) for sub in sys.systems]
    sys = f(sys, fullname)
    return compose(sys, subs)
end

#=
function basename(sys::ODESystem)
    name = ModelingToolkit.get_name(sys) |> string
    prefix = ModelingToolkit.get_namespace(sys) * "₊"
    return chopprefix(name, prefix) |> Symbol
end
=#

function replace(sys::ODESystem, old_new_subsys::Pair{ODESystem, ODESystem})
    old_subsys, new_subsys = old_new_subsys # unpack
    fullname_target = ModelingToolkit.get_name(old_subsys) |> string
    return transform((sys, fullname) -> (fullname == fullname_target ? new_subsys : identity(sys)), sys)
end

# for testing: transform(identity, sys) should do no harm to a system
function identity(sys::ODESystem)
    eqs = ModelingToolkit.get_eqs(sys)
    ieqs = ModelingToolkit.get_initialization_eqs(sys)
    vars = ModelingToolkit.get_unknowns(sys)
    pars = ModelingToolkit.get_ps(sys)
    defs = ModelingToolkit.get_defaults(sys)
    guesses = ModelingToolkit.get_guesses(sys)
    pars = filter(p -> Symbol(p) != Symbol("DEF"), pars) # TODO: remove once fixed: https://github.com/SciML/ModelingToolkit.jl/issues/3322
    return ODESystem(eqs, t, vars, pars; initialization_eqs=ieqs, defaults=defs, guesses=guesses, name=sys.name, description=sys.description)
end

function debugize(sys::ODESystem)
    return transform((s, _) -> length(s.systems) == 0 ? debug_system(s) : identity(s), sys)
end

O(x, ϵⁿ) = x * ϵⁿ
O(eq::Equation, ϵⁿ) = O(eq.lhs, ϵⁿ) ~ O(eq.rhs, ϵⁿ)
O(ϵⁿ) = x -> O(x, ϵⁿ)

function taylor(sys::ODESystem, ϵ, orders)
    eqs = ModelingToolkit.get_eqs(sys)
    ieqs = ModelingToolkit.get_initialization_eqs(sys)
    vars = ModelingToolkit.get_unknowns(sys)
    pars = ModelingToolkit.get_ps(sys)
    defs = ModelingToolkit.get_defaults(sys)
    guesses = ModelingToolkit.get_guesses(sys)

    # extract requested orders
    eqs = taylor(eqs, ϵ, orders)
    ieqs = taylor(ieqs, ϵ, orders)

    # remove resulting trivial equations
    trivial_eqs = [0 ~ 0, 0 ~ -0.0]
    eqs = filter(eq -> !(eq in trivial_eqs), eqs)
    ieqs = filter(eq -> !(eq in trivial_eqs), ieqs)

    pars = filter(p -> Symbol(p) != Symbol("DEF"), pars) # TODO: remove once fixed: https://github.com/SciML/ModelingToolkit.jl/issues/3322

    return ODESystem(eqs, t, vars, pars; initialization_eqs=ieqs, defaults=defs, guesses=guesses, name=sys.name, description=sys.description)
end

have(sys, s::Symbol) = s in Base.propertynames(sys)
have(s) = !isnothing(s) # shorthand for checking if we have a given species

function spline(y, x)
    dx = diff(x)
    dx[end] == 0 && return spline(y[begin:end-1], x[begin:end-1]) # endpoints are duplicated when ODE solver ends with callback; in that case remove it
    all(diff(x) .< 0) && return spline(reverse(y), reverse(x)) # reverse if monotonically decreasing
    all(diff(x) .> 0) || error("x is not monotonically increasing")
    return CubicSpline(y, x; extrapolation = ExtrapolationType.Linear)
end

# compute dy/dx by splining y(x)
function D_spline(y, x; order = 1)
    y_spline = spline(y, x)
    return map(x -> DataInterpolations.derivative(y_spline, x, order), x)
end

value(s, t) = s(t)
derivative(s, t, order=1) = DataInterpolations.derivative(s, t, order)
@register_symbolic value(s::CubicSpline, t)
@register_symbolic derivative(s::CubicSpline, t, order)

# create a range, optionally skipping the first point
range_until(start, stop, step; skip_start=false) = range(skip_start ? start+step : start, step=step, length=Int(ceil((stop-start)/step+1)))

"""
    extend_array(a::AbstractArray, N)

Add `Nextra` linearly separated points between each point in `a`.
"""
function extend_array(a::AbstractArray, Nextra)
    length(a) >= 1 || error("Cannot extend empty array")
    b = similar(a, length(a) + (length(a) - 1) * Nextra)
    for i in 1:length(a)-1
        i1 = 1+(Nextra+1)*(i-1)
        i2 = i1 + Nextra
        b[i1:i2] .= range(a[i], a[i+1], length = Nextra+2)[begin:end-1]
    end
    b[end] = a[end]
    return b
end

"""
    reduce_array(a::AbstractArray, target_length::Integer)

Reduce the length of the array `a` by keeping only every N-th element such that its final length is roughly `target_length`.
This modified `a`.
"""
function reduce_array!(a::AbstractArray, target_length::Integer)
    if target_length > length(a)
        return a
    end
    N = Int(ceil(length(a) / target_length))
    is = range(1, length(a), step = N)
    a[1:length(is)] .= a[is]
    if is[end] != length(a)
        a[length(is)+1] = a[end]
        return a[1:length(is)+1]
    else
        return a[1:length(is)]
    end
end
