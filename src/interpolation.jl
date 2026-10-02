# Interpolators sample a function at nodes x, and interpolates in a transformed coordinate y = f(x)
abstract type AbstractInterpolator{T} end

struct CubicSplineInterpolator{T, Y, F} <: AbstractInterpolator{T}
    xs::Vector{T} # points in input domain: x = f⁻¹(y) (e.g. wavenumbers k)
    ys::Vector{Y} # points in interpolation domain: y = f(x)
    f::F
end

struct BarycentricInterpolator{T, W, F} <: AbstractInterpolator{T}
    xs::Vector{T} # points in input domain: x = f⁻¹(y) (e.g. wavenumbers k)
    ys::Vector{T} # points in interpolation domain: y = f(x)
    ws::Vector{W} # Barycentric interpolation weights
    f::F
end

struct PiecewiseInterpolator{T, P <: Tuple} <: AbstractInterpolator{T}
    pieces::P # interpolators in ascending x-order
    xs::Vector{T} # all unique nodes, in ascending x-order
    iranges::Vector{UnitRange{Int}} # index range into xs for each piece
end

function CubicSplineInterpolator(xs; f = identity)
    issorted(xs) || throw(ArgumentError("Input points must be sorted in ascending order"))
    xs = collect(xs) # to array
    ys = f.(xs)
    return CubicSplineInterpolator(xs, ys, f)
end

# Map nodes ys in the interpolation domain [f(xmin), f(xmax)] back to the input domain [xmin, xmax]
function inverse_nodes(ys, xmin, xmax, f, f⁻¹)
    issorted(ys) || throw(ArgumentError("Domain transformation f(x) is not monotonically increasing"))
    if isnothing(f⁻¹)
        f⁻¹ = f == identity ? identity : y -> solve(IntervalNonlinearProblem((x, _) -> f(x) - y, (xmin, xmax))).u # invert numerically
    end
    xs = f⁻¹.(ys)
    xs[begin] ≈ xmin && xs[end] ≈ xmax || throw(ArgumentError("f(x) and f⁻¹(x) are not inverses"))
    xs[begin], xs[end] = xmin, xmax # prevent floating point bounds errors from f⁻¹(f(k))
    return xs
end

function EquispacedInterpolator(xmin, xmax, order; f = identity, f⁻¹ = nothing)
    xmax > xmin || throw(ArgumentError("Interval $((xmin, xmax)) is not sorted"))
    ys = collect(lingrid(f(xmin), f(xmax); length = order + 1))
    xs = inverse_nodes(ys, xmin, xmax, f, f⁻¹)
    ws = eltype(ys)[(-1)^j * binomial(order, j) for j in 0:order]
    return BarycentricInterpolator(xs, ys, ws, f)
end

function ChebyshevInterpolator(xmin, xmax, order; f = identity, f⁻¹ = nothing)
    xmax > xmin || throw(ArgumentError("Interval $((xmin, xmax)) is not sorted"))
    ys = chebgrid(f(xmin), f(xmax); order)
    xs = inverse_nodes(ys, xmin, xmax, f, f⁻¹)
    ws = eltype(ys)[(-1)^j for j in 0:order]
    ws[begin] /= 2
    ws[end] /= 2
    return BarycentricInterpolator(xs, ys, ws, f)
end

function ChebyshevIntegerInterpolator(xmin, xmax, order::Integer)
    xmax > xmin || throw(ArgumentError("Interval $((xmin, xmax)) is not sorted"))
    order ≥ 1 || throw(ArgumentError("Order must be ≥ 1, got $order"))
    xs = round.(Int, chebgrid(xmin, xmax; order)) # round each Chebyshev point to its nearest integer
    allunique(xs) || throw(ArgumentError(
        "Integer-rounded Chebyshev nodes on ($xmin, $xmax) of order $order collide. Reduce the order or widen the interval."
    ))
    return BarycentricInterpolator(xs, xs, baryweights(xs), identity)
end

function PiecewiseChebyshevInterpolator(xbreaks, orders; f = identity, f⁻¹ = nothing)
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

    pieces = ntuple(j -> ChebyshevInterpolator(xbreaks[j], xbreaks[j+1], orders[j]; f = f[j], f⁻¹ = f⁻¹[j]), N)
    return PiecewiseInterpolator(pieces...)
end

# Combine interpolators on adjacent intervals; a node shared by neighboring pieces is stored (and sampled) only once
function PiecewiseInterpolator(pieces::AbstractInterpolator...)
    all(pieces[j][end] ≤ pieces[j+1][begin] for j in 1:length(pieces)-1) || throw(ArgumentError("Pieces must be sorted and non-overlapping"))
    xs = unique(vcat((piece.xs for piece in pieces)...)) # promotes types; drops duplicate shared endpoints
    iranges = [searchsortedfirst(xs, piece[begin]):searchsortedlast(xs, piece[end]) for piece in pieces]
    return PiecewiseInterpolator(pieces, xs, iranges)
end

Base.eltype(::Type{<:AbstractInterpolator{T}}) where {T} = T # type of x-points
Base.extrema(interp::AbstractInterpolator) = (minimum(interp), maximum(interp))
Base.firstindex(interp::AbstractInterpolator) = firstindex(interp.xs)
Base.lastindex(interp::AbstractInterpolator) = lastindex(interp.xs)
Base.minimum(interp::AbstractInterpolator) = interp[begin]
Base.maximum(interp::AbstractInterpolator) = interp[end]
Base.getindex(interp::AbstractInterpolator, i::Int) = interp.xs[i]
Base.iterate(interp::AbstractInterpolator, args...; kwargs...) = iterate(interp.xs, args...; kwargs...)
Base.length(interp::AbstractInterpolator) = length(interp.xs)
order(interp::AbstractInterpolator) = length(interp) - 1

# Compute Barycentric interpolation weights wᵢ = 1 / ∏_{j≠i}(xᵢ - xⱼ) for arbitrary points
# See https://people.maths.ox.ac.uk/trefethen/barycentric.pdf (section 7)
# TODO: consider exp(sum(log(...))) trick in https://github.com/chebfun/chebfun/blob/master/baryWeights.m
function baryweights(x::AbstractVector)
    C = 4 / (maximum(x) - minimum(x)) # capacity: used to keep product close to 1
    w = [1 / prod(C*(x[i]-x[j]) for j in eachindex(x) if j != i) for i in eachindex(x)]
    w ./= maximum(abs, w) # normalize so largest weight is 1 for stability
    return w
end

# Barycentric interpolation formula https://epubs.siam.org/doi/10.1137/S0036144502417715
function barycentric(ys, ws, vals, y)
    num = zero(eltype(vals)) # promote e.g. integer input to float weights
    den = zero(typeof(y))
    @fastmath @inbounds for j in eachindex(ys)
        dy = y - ys[j]
        iszero(dy) && return vals[j]
        t = ws[j] / dy
        num += t * vals[j]
        den += t
    end
    return num / den
end

# Interpolating function of y = f(x) through values vals at the nodes
interpolant(interp::CubicSplineInterpolator, vals) = CubicSpline(vals, interp.ys)
interpolant(interp::BarycentricInterpolator, vals) = y -> barycentric(interp.ys, interp.ws, vals, y)

(interp::AbstractInterpolator)(vals::AbstractVector, x) = interpolant(interp, vals).(interp.f.(x))

# Index of the piece containing x (the left one at shared endpoints)
pieceindex(interp::PiecewiseInterpolator, x) = something(findfirst(piece -> x ≤ maximum(piece), interp.pieces), lastindex(interp.pieces))

function (interp::PiecewiseInterpolator)(vals::AbstractVector, x::Number)
    j = pieceindex(interp, x)
    return interp.pieces[j](@view(vals[interp.iranges[j]]), x)
end
function (interp::PiecewiseInterpolator)(vals::AbstractVector, xs::AbstractArray)
    out = similar(vals, size(xs))
    js = pieceindex.(Ref(interp), xs)
    for j in eachindex(interp.pieces)
        in_piece = js .== j
        out[in_piece] .= interp.pieces[j](@view(vals[interp.iranges[j]]), xs[in_piece]) # evaluate each piece once on all its points
    end
    return out
end

# Interpolation is linear in vals, so interp(vals, x) == W * vals; column j of W interpolates the j-th unit vector
function interpolation_matrix(interp::AbstractInterpolator, x::AbstractVector; thread = true)
    T = float(promote_type(eltype(interp), eltype(x)))
    W = Matrix{T}(undef, length(x), length(interp))
    @tasks for j in axes(W, 2)
        @set scheduler = thread ? :dynamic : :serial
        @local e = zeros(T, length(interp)) # one unit vector per task
        e[j] = 1
        W[:, j] .= interp(e, x)
        e[j] = 0
    end
    return W
end

# Whether type E consists only of scalars of type R (e.g. SVector{N, R} or ForwardDiff.Dual{Tag, R})
isflat(::Type{E}, ::Type{R}) where {E, R} = E === R || (isbitstype(E) && isstructtype(E) && all(T -> isflat(T, R), fieldtypes(E)))

# Interpolate each row of V (values at the nodes along columns) with one matrix multiplication
function (interp::AbstractInterpolator)(V::AbstractMatrix, x::AbstractVector; thread = true)
    W = interpolation_matrix(interp, x; thread)
    R = eltype(W)
    V isa Array && isflat(eltype(V), R) || return V * transpose(W)
    # multiply elements as consecutive scalars in one real matrix (much faster than with e.g. SVector or Dual elements)
    out = similar(V, size(V, 1), length(x))
    flat(A) = reshape(reinterpret(R, A), :, size(A, 2))
    mul!(flat(out), flat(V), transpose(W))
    return out
end

interpolate(interp::AbstractInterpolator, vals, x) = interp(vals, x)
interpolate(xs::AbstractVector, vals, x) = interpolate(CubicSplineInterpolator(xs), vals, x) # use cubic splines when only passing an array

Base.show(io::IO, interp::CubicSplineInterpolator) = print(io, "Cubic spline interpolator: type = $(eltype(interp)), domain = $(extrema(interp)), order = $(order(interp))")
Base.show(io::IO, interp::BarycentricInterpolator) = print(io, "Barycentric polynomial interpolator: type = $(eltype(interp)), domain = $(extrema(interp)), order = $(order(interp))")
Base.show(io::IO, interp::PiecewiseInterpolator) = print(io, "Piecewise interpolator: type = $(eltype(interp)), domain = $(join(extrema.(interp.pieces), " + ")), order = $(join(order.(interp.pieces), " + "))")
