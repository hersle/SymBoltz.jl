import DataInterpolations: AbstractInterpolation, CubicSpline, CubicHermiteSpline
import Symbolics: taylor, operation, sorted_arguments, unwrap
import Base: identity, replace
using QuadGK
using ModelingToolkit: get_description, get_systems

∫(f, a, b) = quadgk(f, a, b)[1]
∫(f, w) = sum(w .*  f) # ≈ ∫f(x)dx over weights found from QuadGK.gauss()

smoothifelse(x, v1, v2; k=1) = 1/2 * ((v1+v2) + (v2-v1)*tanh(k*x)) # smooth transition/step function from v1 at x<0 to v2 at x>0; smoothifelse(x, v1, v2, k→∞) → ifelse(x < 0, v1, v2)

function transform(f::Function, sys::System; fullname=nameof(sys))
    subs = [transform(f, sub; fullname = ModelingToolkit.iscomplete(sys) ? Symbol() : Symbol(fullname, :₊, nameof(sub))) for sub in get_systems(sys)]
    sys = f(sys, fullname)
    return compose(sys, subs)
end

function replace(sys::System, old_new_subsys::Pair{System, System})
    old_subsys, new_subsys = old_new_subsys # unpack
    fullname_target = ModelingToolkit.get_name(old_subsys) |> string
    return transform((sys, fullname) -> (fullname == fullname_target ? new_subsys : identity(sys)), sys)
end

# for testing: transform(identity, sys) should do no harm to a system
function identity(sys::System)
    iv = ModelingToolkit.get_iv(sys)
    eqs = ModelingToolkit.get_eqs(sys)
    ieqs = ModelingToolkit.get_initialization_eqs(sys)
    vars = ModelingToolkit.get_unknowns(sys)
    pars = ModelingToolkit.get_ps(sys)
    defs = ModelingToolkit.get_defaults(sys)
    guesses = ModelingToolkit.get_guesses(sys)
    return System(eqs, iv, vars, pars; initialization_eqs=ieqs, defaults=defs, guesses=guesses, name=nameof(sys), description=get_description(sys))
end

function debugize(sys::System)
    return transform((s, _) -> length(get_systems(s)) == 0 ? debug_system(s) : identity(s), sys)
end

function isbackground(expr)
    vars = ModelingToolkit.vars(expr; op = Nothing) # don't collect Differential(τ)(f(τ, k))
    for var in vars
        if iscall(var)
            if operation(var) === getindex
                var = arguments(var)[1] # e.g. F(τ, k)[1] to F(τ, k)
            end
            for arg in arguments(var)
                isequal(arg, k) && return false # function of k, e.g. f(τ, k)?
            end
        end
    end
    return true
end
function isperturbation(expr)
    return true # function of k⁰ or k¹? always yes
end

function filter_system(f::Function, sys::System)
    iv = ModelingToolkit.get_iv(sys)
    eqs = ModelingToolkit.get_eqs(sys)
    ieqs = ModelingToolkit.get_initialization_eqs(sys)
    vars = ModelingToolkit.get_unknowns(sys)
    pars = ModelingToolkit.get_ps(sys)
    defs = ModelingToolkit.get_defaults(sys)
    guesses = ModelingToolkit.get_guesses(sys)

    # extract requested orders
    eqs = filter(f, eqs)
    ieqs = filter(f, ieqs)
    vars = filter(f, vars)
    pars = filter(!isinitial, pars) # remove Initial(...) # TODO: shouldn't have to touch this
    defs = filter(pair -> #=!isinitial(pair.first) &&=# f(pair.first), defs)

    return System(eqs, iv, vars, pars; initialization_eqs=ieqs, defaults=defs, guesses=guesses, name=nameof(sys), description=get_description(sys))

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
    return CubicHermiteSpline(ẏ, y, x) # TODO: use PCHIP instead? https://docs.sciml.ai/DataInterpolations/stable/methods/#PCHIP-Interpolation
end

function spline(sol::ODESolution)
    ts = sol.t
    Nu, _ = size(sol)
    us = map(u -> SVector{Nu}(u), sol(ts, Val{0}))
    dus = map(u -> SVector{Nu}(u), sol(ts, Val{1}))
    return CubicHermiteSpline(dus, us, ts; extrapolation = ExtrapolationType.Extension) # TODO: use PCHIP instead? https://docs.sciml.ai/DataInterpolations/stable/methods/#PCHIP-Interpolation
end

# TODO: takes up a lot of time in solvept; refactor so all splines are computed simultaneously for the same ODE time t
value(s::AbstractInterpolation, t, i) = s(t)[i]
@register_symbolic value(s::AbstractInterpolation, t, i)

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

# TODO: generate_jacobian fails on systems returned from this function
# TODO: Use MTKStdLib Interpolation blocks? https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/input_component/#Interpolation-Block
function mtkcompile_spline(sys::System, vars)
    vars = ModelingToolkit.unwrap.(vars)

    # Build mapping from variables to spline parameters
    splname = :bgspline
    spl, = @parameters $splname::CubicHermiteSpline

    # Replace variable equations in system by spline evaluations
    function spline(sys::System)
        # TODO: initialization_eqs? guesses?
        iv = ModelingToolkit.get_iv(sys)
        eqs = ModelingToolkit.get_eqs(sys) # will be modified in-place
        obs = ModelingToolkit.get_observed(sys) # will be modified in-place

        diffvar(expr) = operation(expr) == D ? diffvar(only(sorted_arguments(expr))) : expr # return e.g. a(τ) for D(D(a(τ)))

        for (vari, var) in enumerate(vars)
            # Find and replace the unknown equation by observed spline evaluation
            i = findfirst(eq -> isequal(diffvar(eq.lhs), var), eqs)
            isnothing(i) && error("$var is not an unknown in the system $(nameof(sys))")
            deleteat!(eqs, i) # delete unknown equation
            insert!(obs, 1, var ~ value(spl, iv, vari)) # add observed equation (at the top, for safety, since observed equations are already topsorted in mtkcompile; alternatively consider calling ModelingToolkit.topsort_equations)
        end

        # Remove splined variables from unknowns (they no longer need to be solved for)
        unknowns = ModelingToolkit.get_unknowns(sys)
        unknowns = filter(var -> !Symbolics.contains_var(var, vars), unknowns)
        @set! sys.unknowns = unknowns

        # Add splines as system parameters
        ps = ModelingToolkit.get_ps(sys)
        ps = [ps; spl]
        @set! sys.ps = ps

        # Do not solve for splined variables during initialization, and add dummy defaults for all splines
        defs = ModelingToolkit.get_defaults(sys)
        defs = remove_initial_conditions!(defs, vars)
        nans = SVector{length(vars)}(zeros(length(vars)))
        spldummy = CubicHermiteSpline([nans, nans], [nans, nans], [0.0, 1.0])
        defs = merge(defs, Dict(spl => spldummy))
        @set! sys.defaults = defs

        return sys
    end

    sys = mtkcompile(sys; additional_passes = [spline])

    return sys, spl
end

function remove_initial_conditions!(ics::Dict, vars; maxorder=2)
    for var in vars
        for order in 0:maxorder
            delete!(ics, (D^order)(var)) # remove any existing initial conditions (would cause overdetermined initialization system)
        end
    end
    return ics
end
