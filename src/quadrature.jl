using FFTW
using QuadGK

struct Quadrature{T, U}
    x::Vector{T} # integration points on [-1, +1]
    w::Vector{T} # integration weights on [-1, +1]
    a::U
    b::U
    name::Symbol

    function Quadrature(x::AbstractVector{T}, w::AbstractVector{T}, interval = (-1.0, 1.0); name = Symbol()) where {T}
        all(-1 .≤ x .≤ 1) || throw(ArgumentError("Quadrature points must be on the canonical interval [-1, 1]"))
        issorted(x) || throw(ArgumentError("Quadrature points must be sorted in ascending order"))
        sum(w) ≈ 2 || throw(ArgumentError("Quadrature weights must sum to 2, but sums to $(sum(w))"))
        length(x) == length(w) || throw(ArgumentError("Quadrature nodes and weights must have same length"))
        new{T, eltype(interval)}(x, w, interval[begin], interval[end], name)
    end
end

function TrapezoidalQuadrature(N::Integer, args...)
    N ≥ 2 || throw(ArgumentError("Trapezoidal rule needs at least 2 points"))
    x = collect(range(-1, 1, length = N))
    w = fill(2/(N-1), N)
    w[begin] = 1/(N-1)
    w[end] = 1/(N-1)
    return Quadrature(x, w, args...; name = Symbol("Trapezoidal"))
end

function TrapezoidalQuadrature(x::AbstractArray, args...)
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
    return Quadrature(x, w, (xmin, xmax), args...; name = Symbol("Trapezoidal"))
end

function SimpsonQuadrature(x::AbstractArray, args...)
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
    return Quadrature(x, w, (xmin, xmax), args...; name = Symbol("Simpson"))
end

SimpsonQuadrature(N::Integer, args...) = SimpsonQuadrature(range(-1, 1, length = N), args...) # reduces to the standard uniform-grid rule (1, 4, 2, 4, 2, ..., 4, 1)

function ClenshawCurtisQuadrature(N::Integer, args...)
    N ≥ 2 || throw(ArgumentError("Clenshaw-Curtis quadrature needs at least 2 points"))
    n = N - 1 # number of FFT points
    x = [-cos(π*m/n) for m in 0:n] # Chebyshev nodes in ascending order
    v = [1 / (1 - 4*min(m, n-m)^2) for m in 0:n-1] # 1/(1-4k^2) and it's mirror image
    w = similar(x)
    w[begin:end-1] .= 2/n .* real(fft(v)) # Discrete Cosine Transform (DCT) for O(N log N) time instead of O(N²)
    w[begin] = w[end] = w[begin]/2 # modify endpoint factors
    return Quadrature(x, w, args...; name = Symbol("Clenshaw-Curtis"))
end

function GaussQuadrature(N::Integer, args...)
    N ≥ 1 || throw(ArgumentError("The number of Gauss quadrature points must be positive"))
    x, w = QuadGK.gauss(N)
    return Quadrature(x, w, args...; name = Symbol("Gauss"))
end

function GaussKronrodQuadrature(N::Integer, args...)
    N ≥ 3 && isodd(N) || throw(ArgumentError("The number of Gauss-Kronrod quadrature points must be at least 3 and odd"))
    x, w = QuadGK.kronrod(div(N-1, 2)) # only returns left half
    x = [x; -reverse(x[begin:end-1])] # add right half (antisymmetric)
    w = [w; reverse(w[begin:end-1])] # add right half (symmetric)
    return Quadrature(x, w, args...; name = Symbol("Gauss-Kronrod"))
end

transform(q::Quadrature, interval) = Quadrature(q.x, q.w, interval; name = q.name) # arrays are reuse and shared!
Base.eltype(::Quadrature{T, U}) where {T, U} = Base.promote_type(T, U)
Base.nameof(q::Quadrature) = q.name
Base.show(io::IO, q::Quadrature) = print(io, q.name == Symbol() ? "Q" : "$(q.name) q", "uadrature rule: on [$(q.a), $(q.b)], $(length(q)) points, eltype = $(eltype(q))")
Base.length(q::Quadrature) = length(q.x)
Base.eachindex(q::Quadrature) = eachindex(q.x)
Base.:(==)(q1::Quadrature, q2::Quadrature) = q1.x == q2.x && q1.w == q2.w && q1.a == q2.a && q1.b == q2.b
Base.:(≈)(q1::Quadrature, q2::Quadrature) = q1.x ≈ q2.x && q1.w ≈ q2.w && q1.a ≈ q2.a && q1.b ≈ q2.b
node(q::Quadrature, i) = (q.b+q.a)/2 + (q.b-q.a)/2 * q.x[i]
nodes(q::Quadrature) = (q.b+q.a)/2 .+ (q.b-q.a)/2 .* q.x
weights(q::Quadrature) = (q.b-q.a)/2 .* q.w # weights for integration on [a, b]
@inbounds @fastmath function integrate(q::Quadrature, f::AbstractVector)
    length(f) == length(q) || throw(ArgumentError("$(length(f))-point f is incompatible with $(length(q))-point quadrature rule"))
    return (q.b-q.a)/2 * sum(q.w[i] * f[i] for i in eachindex(q)) # integrate samples of f on [a, b]
end
@inbounds @fastmath function integrate(q::Quadrature, f::Function)
    return (q.b-q.a)/2 * sum(q.w[i] * f(node(q, i)) for i in eachindex(q)) # integrate f(x) on [a, b]
end
(q::Quadrature)(args...) = integrate(q, args...) # calling is equivalent to integrating
