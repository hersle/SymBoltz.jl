import DataInterpolations: AbstractInterpolation, CubicSpline, CubicHermiteSpline
import SpecialFunctions: zeta as ζ
import Symbolics: taylor, operation, sorted_arguments, unwrap
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
    eqs = taylor(eqs, ϵ, orders) # TODO: want taylor_coeff?
    ieqs = taylor(ieqs, ϵ, orders) # TODO: want taylor_coeff?

    # remove resulting trivial equations
    trivial_eqs = [0 ~ 0, 0 ~ -0.0]
    eqs = filter(eq -> !(eq in trivial_eqs), eqs)
    ieqs = filter(eq -> !(eq in trivial_eqs), ieqs)

    pars = filter(p -> Symbol(p) != Symbol("DEF"), pars) # TODO: remove once fixed: https://github.com/SciML/ModelingToolkit.jl/issues/3322

    return ODESystem(eqs, t, vars, pars; initialization_eqs=ieqs, defaults=defs, guesses=guesses, name=sys.name, description=sys.description)
end

have(sys, s::Symbol) = s in nameof.(ModelingToolkit.get_systems(sys))
have(s) = !isnothing(s) # shorthand for checking if we have a given species

function spline(y, x)
    dx = diff(x)
    dx[end] == 0 && return spline(y[begin:end-1], x[begin:end-1]) # endpoints are duplicated when ODE solver ends with callback; in that case remove it
    all(diff(x) .< 0) && return spline(reverse(y), reverse(x)) # reverse if monotonically decreasing
    all(diff(x) .> 0) || error("x is not monotonically increasing")
    return CubicSpline(y, x; extrapolation = ExtrapolationType.Linear)
end

function spline(sol::ODESolution, var, dvar = nothing)
    t = sol.t
    y = sol(t, Val{0}; idxs=var).u
    if isnothing(dvar)
        ẏ = sol(t, Val{1}; idxs=var).u # evaluates derivative from ODE f (?)
    else
        ẏ = sol(t, Val{0}; idxs=dvar).u # evaluates derivative from custom expression
    end
    return CubicHermiteSpline(ẏ, y, t; extrapolation = ExtrapolationType.Linear)
end

value(s::AbstractInterpolation, t) = s(t)
derivative(s::AbstractInterpolation, t, order=1) = DataInterpolations.derivative(s, t, order)
@register_symbolic value(s::AbstractInterpolation, t)
@register_symbolic derivative(s::AbstractInterpolation, t, order)

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

# TODO: handle higher-order (3+) derivatives
function structural_simplify_spline(sys::ODESystem, vars; verbose = false)
    isempty(ModelingToolkit.get_systems(sys)) || error("System must be flattened")

    eqs = ModelingToolkit.get_eqs(sys)
    ieqs = ModelingToolkit.get_initialization_eqs(sys)
    pars = ModelingToolkit.get_ps(sys)
    defs = ModelingToolkit.get_defaults(sys)
    guesses = ModelingToolkit.get_guesses(sys)

    # 1) Recreate sys with spline parameters
    dummyspline = CubicHermiteSpline([NaN, NaN, NaN], [NaN, NaN, NaN], 0.0:1.0:2.0)
    extraeqs = Equation[]
    var2spl = Dict() # list of additional spline parameters
    varnames = nameof.(operation.(unwrap.(vars))) # names of variables to spline
    for (i, eq) in enumerate(eqs) # TODO: iterate over vars instead? needs a little different buildup
        var = eq.lhs
        op = operation(var)
        if op == D # is time derivative
            var = only(sorted_arguments(var)) # e.g. a(t)
            varname = nameof(operation(var)) # e.g. :a
            if varname in varnames
                splname = Symbol(varname, :_spline) # e.g. :a_spline
                spl = only(@parameters $splname::CubicHermiteSpline)
                verbose && println("Splining $var -> $spl")
                eqs[i] = var ~ value(spl, t) # replace equation with splined version
                delete!(defs, var) # remove any existing initial conditions (would cause overdetermined initialization system)
                var2spl[var] = spl # add spline parameter
                defs = merge(defs, Dict(spl => dummyspline))
            end
        end
    end
    @set! sys.eqs = [eqs; extraeqs]

    # 2) Let structural_simplify do its magic, including inserting dummy derivative expressions
    sys = structural_simplify(sys)

    # 3) Finally substitute derivatives for the spline for the dummy derivatives
    eqs = ModelingToolkit.get_eqs(sys)
    obs = ModelingToolkit.get_observed(sys)
    for spl in values(var2spl)
        for order in reverse(1:2) # start with highest-order derivative (so D(D(y)) is handled as a whole, instead of the inner D(y) first)
            for list in [eqs, obs]
                for i in eachindex(list)
                    list[i] = substitute(list[i], (D^order)(value(spl, t)) => derivative(spl, t, order)) # overwrite equations
                end
            end
        end
    end

    # Return splined system and variable-to-parameter map
    vars = ModelingToolkit.get_unknowns(sys)
    pars = [pars; collect(values(var2spl))] # add splines as system parameters
    pars = filter(p -> Symbol(p) != Symbol("DEF"), pars) # TODO: remove once fixed: https://github.com/SciML/ModelingToolkit.jl/issues/3322
    sys = ODESystem([eqs; obs], t, vars, pars; initialization_eqs=ieqs, defaults=defs, guesses=guesses, name=sys.name, description=sys.description) # copy structurally simplified system into a non-simplified one (2x simplify is forbidden)
    return sys, var2spl
end
