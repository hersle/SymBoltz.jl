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

# proxy function for evaluating a spline
@register_symbolic spleval(x, spline::CubicSpline)
spleval(x, spline) = spline(x)

function spline(y, x)
    i = unique(i -> x[i], eachindex(x)) # get indices of unique x values
    x, y = x[i], y[i] # remove duplicate x values
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