using SymbolicIndexingInterface: variable_index

# callback for terminating an integrator when var == val0
function callback_terminator(sys, var, val0)
    varindex = variable_index(sys, var)
    return ContinuousCallback((u, _, _) -> (val = u[varindex]; val - val0), terminate!)
end

# proxy function for evaluating a spline
@register_symbolic spleval(x, spline::CubicSpline)
spleval(x, spline) = spline(x)

# compute dy/dx by splining y(x)
function D_spline(y, x; Spline = CubicSpline, order = 1)
    y_spline = Spline(y, x; extrapolate=true)
    y′ = DataInterpolations.derivative.(Ref(y_spline), x, order)
    return y′
end

# create a range, optionally skipping the first point
range_until(start, stop, step; skip_start=false) = range(skip_start ? start+step : start, step=step, length=Int(ceil((stop-start)/step+1)))