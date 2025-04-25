import DataInterpolations: AbstractInterpolation, CubicSpline, CubicHermiteSpline
import Symbolics: taylor, operation, sorted_arguments, unwrap
import Base: identity, replace
using QuadGK

∫(f, a, b) = quadgk(f, a, b)[1]
∫(f, w) = sum(w .*  f) # ≈ ∫f(x)dx over weights found from QuadGK.gauss()

# callback for terminating an integrator when var >= val0 or var == val0
# TODO: replace with symbolic events?
function callback_today(sys, var, val0; continuous = true, terminate = true, save_positions = (true, false), kwargs...)
    varindex = ModelingToolkit.variable_index(sys, var)
    t0idx = sys.t0
    if continuous
        T = ContinuousCallback
        f = (u, _, _) -> (val = u[varindex]; val - val0)
    else
        T = DiscreteCallback
        f = (u, _, _) -> (val = u[varindex]; val >= val0)
    end
    function affect!(integrator)
        integrator.ps[t0idx] = integrator.t # set time today to time when a == 1
        terminate && terminate!(integrator) # stop integration if desired
    end
    return T(f, affect!; save_positions, kwargs...)
end

function transform(f::Function, sys::ODESystem; fullname=string(sys.name))
    subs = [transform(f, sub; fullname = (ModelingToolkit.iscomplete(sys) ? "" : fullname * "₊") * string(sub.name)) for sub in sys.systems]
    sys = f(sys, fullname)
    return compose(sys, subs)
end

function replace(sys::ODESystem, old_new_subsys::Pair{ODESystem, ODESystem})
    old_subsys, new_subsys = old_new_subsys # unpack
    fullname_target = ModelingToolkit.get_name(old_subsys) |> string
    return transform((sys, fullname) -> (fullname == fullname_target ? new_subsys : identity(sys)), sys)
end

# for testing: transform(identity, sys) should do no harm to a system
function identity(sys::ODESystem)
    iv = ModelingToolkit.get_iv(sys)
    eqs = ModelingToolkit.get_eqs(sys)
    ieqs = ModelingToolkit.get_initialization_eqs(sys)
    vars = ModelingToolkit.get_unknowns(sys)
    pars = ModelingToolkit.get_ps(sys)
    defs = ModelingToolkit.get_defaults(sys)
    guesses = ModelingToolkit.get_guesses(sys)
    pars = filter(p -> Symbol(p) != Symbol("DEF"), pars) # TODO: remove once fixed: https://github.com/SciML/ModelingToolkit.jl/issues/3322
    return ODESystem(eqs, iv, vars, pars; initialization_eqs=ieqs, defaults=defs, guesses=guesses, name=sys.name, description=sys.description)
end

function debugize(sys::ODESystem)
    return transform((s, _) -> length(s.systems) == 0 ? debug_system(s) : identity(s), sys)
end

O(x, ϵⁿ) = x * ϵⁿ
O(eq::Equation, ϵⁿ) = O(eq.lhs, ϵⁿ) ~ O(eq.rhs, ϵⁿ)
O(ϵⁿ) = x -> O(x, ϵⁿ)

function taylor(sys::ODESystem, ϵ, orders; kwargs...)
    iv = ModelingToolkit.get_iv(sys)
    eqs = ModelingToolkit.get_eqs(sys)
    ieqs = ModelingToolkit.get_initialization_eqs(sys)
    vars = ModelingToolkit.get_unknowns(sys)
    pars = ModelingToolkit.get_ps(sys)
    defs = ModelingToolkit.get_defaults(sys)
    guesses = ModelingToolkit.get_guesses(sys)

    # extract requested orders
    eqs = taylor(eqs, ϵ, orders; kwargs...)
    ieqs = taylor(ieqs, ϵ, orders; kwargs...)

    # remove resulting trivial equations
    trivial_eqs = [0 ~ 0, 0 ~ -0.0]
    eqs = filter(eq -> !(eq in trivial_eqs), eqs)
    ieqs = filter(eq -> !(eq in trivial_eqs), ieqs)

    pars = filter(p -> Symbol(p) != Symbol("DEF"), pars) # TODO: remove once fixed: https://github.com/SciML/ModelingToolkit.jl/issues/3322

    return ODESystem(eqs, iv, vars, pars; initialization_eqs=ieqs, defaults=defs, guesses=guesses, name=sys.name, description=sys.description)
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

function spline(y, ẏ, x)
    dx = diff(x)
    all(diff(x) .< 0) && return spline(reverse(y), reverse(ẏ), reverse(x)) # reverse if monotonically decreasing
    all(diff(x) .> 0) || error("x is not monotonically increasing")
    return CubicHermiteSpline(ẏ, y, x)
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

function structural_simplify_spline(sys::ODESystem, vars; maxorder = 2)
    vars = ModelingToolkit.unwrap.(vars)

    # Build mapping from variables to spline parameters
    function spline(var)
        splname = Symbol(nameof(operation(var)), :ˍspline) # e.g. :aˍspline
        spl, = @parameters $splname::CubicHermiteSpline
        return spl
    end
    var2spl = Dict(var => spline(var) for var in vars)

    # Replace variable equations in system by spline evaluations
    function spline(sys::ODESystem)
        # TODO: initialization_eqs? guesses?
        iv = ModelingToolkit.get_iv(sys)
        eqs = ModelingToolkit.get_eqs(sys) # will be modified in-place
        obs = ModelingToolkit.get_observed(sys) # will be modified in-place

        diffvar(expr) = operation(expr) == D ? diffvar(only(sorted_arguments(expr))) : expr # return e.g. a(t) for D(D(a(t)))

        for var in vars
            spl = var2spl[var]

            # Find and replace the unknown equation by observed spline evaluation
            i = findfirst(eq -> isequal(diffvar(eq.lhs), var), eqs)
            isnothing(i) && error("$var is not an unknown in the system $(nameof(sys))")
            deleteat!(eqs, i) # delete unknown equation
            insert!(obs, 1, var ~ value(spl, iv)) # add observed equation (at the top, for safety, since observed equations are already topsorted in structural_simplify; alternatively consider calling ModelingToolkit.topsort_equations)

            # Find and replace any observed derivative expressions with spline derivative evaluations
            for order in 1:maxorder
                dvar = ModelingToolkit.diff2term((D^order)(var))
                i = findfirst(eq -> isequal(eq.lhs, dvar), obs)
                if !isnothing(i)
                    obs[i] = dvar ~ derivative(spl, iv, order)
                end
            end
        end

        # Remove splined variables from unknowns (they no longer need to be solved for)
        unknowns = ModelingToolkit.get_unknowns(sys)
        unknowns = filter(var -> !Symbolics.contains_var(var, vars), unknowns)
        @set! sys.unknowns = unknowns

        # Add splines as system parameters
        ps = ModelingToolkit.get_ps(sys)
        ps = [ps; collect(values(var2spl))]
        @set! sys.ps = ps

        # Do not solve for splined variables during initialization, and add dummy defaults for all splines
        defs = ModelingToolkit.get_defaults(sys)
        defs = remove_initial_conditions!(defs, keys(var2spl))
        defs = merge(defs, Dict(spl => CubicHermiteSpline([NaN, NaN], [NaN, NaN], [0.0, 1.0]) for spl in values(var2spl)))
        @set! sys.defaults = defs

        return sys
    end

    sys = structural_simplify(sys; additional_passes = [spline])

    return sys, var2spl
end

function remove_initial_conditions!(ics::Dict, vars; maxorder=2)
    for var in vars
        for order in 0:maxorder
            delete!(ics, (D^order)(var)) # remove any existing initial conditions (would cause overdetermined initialization system)
        end
    end
    return ics
end
