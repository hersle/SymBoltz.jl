import SpecialFunctions: zeta as ζ
using SymbolicIndexingInterface: variable_index
using QuadGK

∫(f, a, b) = quadgk(f, a, b)[1] # TODO: use integrator that supports dual numbers
∫(f, w) = sum(w .*  f) # ≈ ∫f(x)dx over weights found from QuadGK.gauss()

# callback for terminating an integrator when var == val0
function callback_terminator(sys, var, val0)
    varindex = variable_index(sys, var)
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
function identity(sys)
    eqs = ModelingToolkit.get_eqs(sys)
    ieqs = ModelingToolkit.get_initialization_eqs(sys)
    vars = ModelingToolkit.get_unknowns(sys)
    pars = ModelingToolkit.get_ps(sys)
    defs = ModelingToolkit.get_defaults(sys)
    guesses = ModelingToolkit.get_guesses(sys)
    return ODESystem(eqs, t, vars, pars; initialization_eqs=ieqs, defaults=defs, guesses=guesses, name=sys.name)
end

function debugize(sys)
    return transform((s, _) -> length(s.systems) == 0 ? debug_system(s) : identity(s), sys)
end

O(x, ϵⁿ) = x * ϵⁿ
O(eq::Equation, ϵⁿ) = O(eq.lhs, ϵⁿ) ~ O(eq.rhs, ϵⁿ)
O(ϵⁿ) = x -> O(x, ϵⁿ)

function extract_order(expr, order)
    if order == 0
        return substitute(expr, ϵ => 0)
    else
        expr = Differential(ϵ)(expr) |> expand_derivatives # differentiate away one power of ϵ^order -> order*ϵ^(order-1)
        expr = expr / order # remove prefactor from differentiation
        return extract_order(expr, order - 1)
    end
end

function extract_order(eq::Equation, order)
    return extract_order(eq.lhs, order) ~ extract_order(eq.rhs, order)
end

function extract_order(sys::ODESystem, orders)
    eqs = ModelingToolkit.get_eqs(sys)
    ieqs = ModelingToolkit.get_initialization_eqs(sys)
    vars = ModelingToolkit.get_unknowns(sys)
    pars = ModelingToolkit.get_ps(sys)
    defs = ModelingToolkit.get_defaults(sys)
    guesses = ModelingToolkit.get_guesses(sys)

    # extract requested orders
    eqs = vcat((extract_order.(eqs, order) for order in orders)...)
    ieqs = vcat((extract_order.(ieqs, order) for order in orders)...)

    # remove resulting trivial equations
    trivial_eqs = [0 ~ 0, 0 ~ -0.0]
    eqs = filter(eq -> !(eq in trivial_eqs), eqs)
    ieqs = filter(eq -> !(eq in trivial_eqs), ieqs)

    sys0 = ODESystem(eqs, t, vars, pars; initialization_eqs=ieqs, defaults=defs, guesses=guesses, name=sys.name)
    return sys0
end

function spline(y, x)
    # remove duplicate x values
    i = unique(i -> (x[i], y[i]), eachindex(x)) # indices of unique values
    x, y = x[i], y[i] # pick them out

    # sort x # TODO: add to DataInterpolations.jl?
    i = sortperm(x) # indices that sorts x
    x, y = x[i], y[i]

    return CubicSpline(y, x; extrapolate=true)
end

# compute dy/dx by splining y(x)
function D_spline(y, x; order = 1)
    y_spline = spline(y, x)
    y′ = DataInterpolations.derivative.(Ref(y_spline), x, order)
    return y′
end

# create a range, optionally skipping the first point
range_until(start, stop, step; skip_start=false) = range(skip_start ? start+step : start, step=step, length=Int(ceil((stop-start)/step+1)))