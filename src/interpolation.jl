abstract type AbstractInterpolator{T} end

struct CubicSplineInterpolator{T, F <: Function} <: AbstractInterpolator{T}
    xs::Vector{T} # points in input domain: x = f⁻¹(y) (e.g. wavenumbers k)
    ys::Vector{T} # points in interpolation domain: y = f(x)
    f::F
end

struct EquispacedInterpolator{T <: Real, F} <: AbstractInterpolator{T}
    xs::Vector{T} # points in input domain: x = f⁻¹(y) (e.g. wavenumbers k)
    ys::Vector{T} # points in interpolation domain: y = f(x)
    ws::Vector{T} # Barycentric interpolation weights
    f::F
end

struct ChebyshevInterpolator{T <: Real, F} <: AbstractInterpolator{T}
    xs::Vector{T} # points in input domain: x = f⁻¹(y) (e.g. wavenumbers k)
    ys::Vector{T} # points in interpolation domain: y = f(x)
    ws::Vector{T} # Barycentric interpolation weights
    f::F
end

struct PiecewiseChebyshevInterpolator{T <: Real, G <: Tuple} <: AbstractInterpolator{T}
    subgrids::G # NTuple of ChebyshevInterpolator in ascending x-order
    xs::Vector{T} # all unique coarse x-values, in descending x-order
    ys::Vector{T} # x = y # TODO: generalize with f-transform?
    iranges::Vector{UnitRange{Int}} # index range into xs for each subgrid
end

struct ChebyshevIntegerInterpolator{Tx <: Integer, Tw <: Real, F} <: AbstractInterpolator{Tx}
    xs::Vector{Tx}
    ys::Vector{Tx}
    ws::Vector{Tw}
    f::F
end

function CubicSplineInterpolator(xs; f = identity)
    issorted(xs) || throw(ArgumentError("Input points must be sorted in ascending order"))
    xs = collect(xs) # to array
    ys = f.(xs)
    return CubicSplineInterpolator(xs, ys, f)
end
function CubicSplineInterpolator(xmin, xmax, n)
    xs = range(xmin, xmax, length = n+1)
    return CubicSplineInterpolator(xs)
end

function EquispacedInterpolator(xmin, xmax, order)
    xmax > xmin || throw(ArgumentError("Interval $((xmin, xmax)) is not sorted"))
    xs = lingrid(xmin, xmax; length = order + 1)
    ys = xs

    # Precompute Barycentric interpolation weights
    T = eltype(xs)
    n = length(xs) - 1
    ws = T[binomial(n, j) for j in 0:n]
    ws .*= [iseven(j) ? T(+1) : T(-1) for j in 0:n]

    return EquispacedInterpolator(xs, ys, ws, identity)
end

function ChebyshevInterpolator(xmin, xmax, order; f = identity, f⁻¹ = nothing)
    xmax > xmin || throw(ArgumentError("Interval $((xmin, xmax)) is not sorted"))
    ymin, ymax = f(xmin), f(xmax)
    ys = chebpoints(order, ymin, ymax)
    issorted(ys; rev = true) || throw(ArgumentError("Domain transformation f(x) is not monotonically increasing"))
    if f == identity && isnothing(f⁻¹)
        f⁻¹ = identity
    end
    if isnothing(f⁻¹)
        # invert numerically
        xs = map(ys) do y
            prob = IntervalNonlinearProblem((x, _) -> f(x) - y, (xmin, xmax))
            sol = solve(prob)
            return sol.u
        end
    else
        # invert analytically
        xs = f⁻¹.(ys)
    end
    xs[end] ≈ xmin && xs[begin] ≈ xmax || throw(ArgumentError("f(x) and f⁻¹(x) are not inverses"))
    xs[end], xs[begin] = xmin, xmax  # prevent floating point bounds errors from f⁻¹(f(k))

    # Precompute Barycentric interpolation weights
    T = eltype(xs)
    n = length(xs) - 1
    ws = [iseven(j) ? T(+1) : T(-1) for j in 0:n]
    ws[begin] /= 2
    ws[end] /= 2

    # Change to ascending order
    reverse!(xs)
    reverse!(ys)
    reverse!(ws)

    return ChebyshevInterpolator(xs, ys, ws, f)
end

function PiecewiseChebyshevInterpolator(xbreaks, orders; f = identity, f⁻¹ = identity)
    N = length(orders) # number of piecewise subgrids
    length(xbreaks) == N + 1 || throw(ArgumentError("Need $(N+1) x-breaks for $N intervals, got $(length(xbreaks))"))
    if !(f isa Tuple)
        f = ntuple(_ -> f, N)
    end
    if !(f⁻¹ isa Tuple)
        f⁻¹ = ntuple(_ -> f⁻¹, N)
    end
    length(f)  == N || throw(ArgumentError("Need $N f, got $(length(f))"))
    length(f⁻¹) == N || throw(ArgumentError("Need $N f⁻¹, got $(length(f⁻¹))"))

    subgrids = ntuple(j -> ChebyshevInterpolator(xbreaks[j], xbreaks[j+1], orders[j]; f = f[j], f⁻¹ = f⁻¹[j]), N)
    xs = reduce(vcat, subgrids[j].xs[2:end] for j in 2:N; init = subgrids[begin].xs) # combine unique x-points in descending order (boundaries share x points)
    iranges = Vector{UnitRange{Int}}(undef, N)
    i = 1
    for j in 1:N
        n = length(subgrids[j].xs)
        iranges[j] = i : i + n - 1 # index range into xs corresponding to subgrid j
        i += n - 1
    end
    return PiecewiseChebyshevInterpolator{eltype(xs), typeof(subgrids)}(subgrids, xs, xs, iranges)
end

function ChebyshevIntegerInterpolator(xmin, xmax, order::Integer)
    xmax > xmin || throw(ArgumentError("Interval $((xmin, xmax)) is not sorted"))
    order ≥ 1 || throw(ArgumentError("Order must be ≥ 1, got $order"))

    xs = reverse!(chebpoints(order, xmin, xmax)) # start with Chebyshev points
    xs = round.(Int, xs) # round each point to its nearest integer
    allunique(xs) || throw(ArgumentError(
        "Integer-rounded Chebyshev nodes on ($xmin, $xmax) of order $order collide. Reduce the order or widen the interval."
    ))

    ys = xs
    ws = baryweights(ys)

    return ChebyshevIntegerInterpolator(xs, ys, ws, identity)
end

Base.eltype(::Type{<:AbstractInterpolator{T}}) where {T} = T # type of x-points
Base.extrema(interp::AbstractInterpolator) = (minimum(interp), maximum(interp))
Base.firstindex(interp::AbstractInterpolator) = firstindex(interp.xs)
Base.lastindex(interp::AbstractInterpolator) = lastindex(interp.xs)
Base.minimum(interp::AbstractInterpolator) = interp[begin]
Base.maximum(interp::AbstractInterpolator) = interp[end]
Base.getindex(interp::AbstractInterpolator, i::Int) = interp.xs[i]
Base.iterate(interp::AbstractInterpolator, args...; kwargs...) = iterate(interp.xs, args...; kwargs...)

# Compute Barycentric interpolation weights wᵢ = 1 / ∏_{j≠i}(xᵢ - xⱼ) for arbitrary points
# See https://people.maths.ox.ac.uk/trefethen/barycentric.pdf (section 7)
# TODO: consider exp(sum(log(...))) trick in https://github.com/chebfun/chebfun/blob/master/baryWeights.m
function baryweights(x::AbstractVector)
    C = 4 / (maximum(x) - minimum(x)) # capacity: used to keep product close to 1
    w = [1 / prod(C*(x[i]-x[j]) for j in eachindex(x) if j != i) for i in eachindex(x)]
    w ./= maximum(abs, w) # normalize so largest weight is 1 for stability
    return w
end

Base.length(interp::AbstractInterpolator) = length(interp.xs)
order(interp::AbstractInterpolator) = length(interp) - 1

# Barycentric interpolation formula https://epubs.siam.org/doi/10.1137/S0036144502417715
function (interp::AbstractInterpolator)(f::AbstractVector, y::Number)
    # promote e.g. integer input to float weights
    num = zero(eltype(f))
    den = zero(typeof(y))
    @fastmath @inbounds for j in eachindex(interp.ys)
        dy = y - interp.ys[j]
        iszero(dy) && return f[j]
        t = interp.ws[j] / dy
        num += t * f[j]
        den += t
    end
    return num / den
end

function (interp::CubicSplineInterpolator)(f::AbstractVector, y::AbstractArray)
    return CubicSpline(f, interp.ys)(y)
end

function (interp::PiecewiseChebyshevInterpolator)(f::AbstractVector, ys_fine::AbstractVector)
    out = similar(f, length(ys_fine))
    for j in eachindex(interp.subgrids)
        subgrid = interp.subgrids[j]
        ymin, ymax = extrema(subgrid)
        in_range = findall(y -> ymin ≤ y ≤ ymax, ys_fine)
        irange = interp.iranges[j]
        for i in in_range
            out[i] = subgrid(f[irange], ys_fine[i])
        end
    end
    return out
end

function (interp::AbstractInterpolator)(f::AbstractVector, x::AbstractArray)
    return interp.(Ref(f), x)
end

function (interp::AbstractInterpolator)(out::AbstractVector, f::AbstractVector, x::AbstractArray)
    length(out) == length(x) || throw(ArgumentError("out and x have different lengths $(length(out)) and $(length(x))"))
    for i in eachindex(x)
        out[i] = interp(f, x[i])
    end
    return out
end

interpolate(x::AbstractInterpolator, y, x′) = x(y, x.f.(x′))
interpolate(x::PiecewiseChebyshevInterpolator, y, x′) = x(y, x′) # no f field
interpolate(x::AbstractVector, y, x′) = interpolate(CubicSplineInterpolator(x), y, x′)

Base.show(io::IO, interp::CubicSplineInterpolator) = print(io, "Cubic spline interpolator: domain = $(extrema(interp)), order = $(order(interp))")
Base.show(io::IO, interp::EquispacedInterpolator) = print(io, "Equispaced polynomial interpolator: domain = $(extrema(interp)), order = $(order(interp))")
Base.show(io::IO, interp::ChebyshevInterpolator) = print(io, "Chebyshev polynomial interpolator: domain = $(extrema(interp)), order = $(order(interp))")
Base.show(io::IO, interp::PiecewiseChebyshevInterpolator) = print(io, "Piecewise Chebyshev polynomial interpolator: domain = $(join(extrema.(interp.subgrids), " + ")), order = $(join(order.(interp.subgrids), " + "))")
Base.show(io::IO, interp::ChebyshevIntegerInterpolator) = print(io, "Integer-rounded Chebyshev interpolator: domain = $(extrema(interp)), order = $(order(interp))")
