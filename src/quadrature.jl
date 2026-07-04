struct Quadrature{T}
    x::Vector{T} # integration points on [-1, +1]
    w::Vector{T} # integration weights on [-1, +1]
    name::Symbol

    function Quadrature(x::AbstractVector{T}, w::AbstractVector{T}; name = Symbol()) where {T}
        all(-1 .≤ x .≤ 1) || throw(ArgumentError("Quadrature points must be on the canonical interval [-1, 1]"))
        issorted(x) || throw(ArgumentError("Quadrature points must be sorted in ascending order"))
        sum(w) ≈ 2 || throw(ArgumentError("Quadrature weights must sum to 2, but sums to $(sum(w))"))
        length(x) == length(w) || throw(ArgumentError("Quadrature nodes and weights must have same length"))
        new{T}(x, w, name)
    end
end

function TrapezoidalQuadrature(x::AbstractArray)
    N = length(x)
    N ≥ 2 || throw(ArgumentError("Trapezoidal quadrature needs at least 2 points"))
    xmin, xmax = extrema(x)
    x = collect(x)
    x .= -1 .+ 2 .* (x .- xmin) / (xmax - xmin) # normalize to canonical [-1, 1]
    w = zeros(N)
    # each interval contributes (y1+y2)*(x2-x1)/2
    for i in 1:N-1
        dx = x[i+1] - x[i]
        w[i] += dx / 2
        w[i+1] += dx / 2
    end
    return Quadrature(x, w; name = Symbol("Trapezoidal"))
end

TrapezoidalQuadrature(N::Integer) = TrapezoidalQuadrature(range(-1, 1, length = N)) # reduces to the standard uniform-grid rule (1, 2, 2, ..., 2, 1)

function SimpsonQuadrature(x::AbstractArray)
    N = length(x)
    isodd(N) && N ≥ 3 || throw(ArgumentError("Simpson's rule needs an odd number of points ≥ 3"))
    xmin, xmax = extrema(x)
    x = collect(x)
    x .= -1 .+ 2 .* (x .- xmin) / (xmax - xmin) # normalize to canonical [-1, 1]
    w = zeros(N)
    @inbounds for i in 1:2:N-2
        # fit a quadratic through each triple of points (possibly unequally spaced) and integrate it exactly
        x0, x1, x2 = x[i], x[i+1], x[i+2]
        h1 = x1 - x0
        h2 = x2 - x1
        H = h1 + h2
        w[i] += H/6 * (2h1 - h2) / h1
        w[i+1] += H/6 * H^2 / (h1*h2)
        w[i+2] += H/6 * (2h2 - h1) / h2
    end
    return Quadrature(x, w; name = Symbol("Simpson"))
end

SimpsonQuadrature(N::Integer) = SimpsonQuadrature(range(-1, 1, length = N)) # reduces to the standard uniform-grid rule (1, 4, 2, 4, 2, ..., 4, 1)

Base.nameof(q::Quadrature) = q.name
Base.eltype(::Quadrature{T}) where {T} = T
Base.show(io::IO, q::Quadrature) = print(io, q.name == Symbol() ? "Q" : "$(q.name) q", "uadrature rule: $(length(q)) points, eltype = $(eltype(q))")
Base.length(q::Quadrature) = length(q.x)
Base.eachindex(q::Quadrature) = eachindex(q.x)
Base.:(==)(q1::Quadrature, q2::Quadrature) = q1.x == q2.x && q1.w == q2.w
Base.:(≈)(q1::Quadrature, q2::Quadrature) = q1.x ≈ q2.x && q1.w ≈ q2.w
nodes(q::Quadrature, a, b) = (b+a)/2 .+ (b-a)/2 .* q.x # nodes on [a, b]
weights(q::Quadrature, a, b) = (b-a)/2 .* q.w # weights on [a, b]
@inbounds @fastmath function integrate(q::Quadrature, f::AbstractVector, a, b)
    length(f) == length(q) || throw(ArgumentError("$(length(f))-point f is incompatible with $(length(q))-point quadrature rule"))
    return (b-a)/2 * sum(q.w[i] * f[i] for i in eachindex(q)) # integrate samples of f at nodes(q, a, b)
end
@inbounds @fastmath function integrate(q::Quadrature, f::Function, a, b)
    return (b-a)/2 * sum(q.w[i] * f((b+a)/2 + (b-a)/2 * q.x[i]) for i in eachindex(q)) # integrate f(x) on [a, b]
end
(q::Quadrature)(f, a, b) = integrate(q, f, a, b) # calling is equivalent to integrating
