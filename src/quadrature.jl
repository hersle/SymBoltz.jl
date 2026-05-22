using FFTW
using QuadGK

struct QuadratureRule{T, U}
    x::Vector{T} # integration points on [-1, +1]
    w::Vector{T} # integration weights on [-1, +1]
    a::U
    b::U
    name::Symbol

    function QuadratureRule(x::AbstractVector{T}, w::AbstractVector{T}, interval = (-1.0, 1.0); name = Symbol()) where {T}
        all(-1 .≤ x .≤ 1) || throw(ArgumentError("Quadrature points must be on the canonical interval [-1, 1]"))
        sum(w) ≈ 2 || throw(ArgumentError("Quadrature weights must sum to 2, but sums to $(sum(w))"))
        issorted(x) || throw(ArgumentError("Quadrature points must be sorted in ascending order"))
        length(x) == length(w) || throw(ArgumentError("Quadrature rule must have same number of nodes and weights"))
        new{T, eltype(interval)}(x, w, interval[begin], interval[end], name)
    end
end

function TrapezoidalRule(N::Integer, args...)
    N ≥ 2 || throw(ArgumentError("Trapezoidal rule needs at least 2 points"))
    x = collect(range(-1, 1, length = N))
    w = fill(2/(N-1), N)
    w[begin] = 1/(N-1)
    w[end] = 1/(N-1)
    return QuadratureRule(x, w, args...; name = Symbol("Trapezoidal"))
end

function ClenshawCurtisRule(N::Integer, args...)
    N ≥ 2 || throw(ArgumentError("Clenshaw-Curtis quadrature needs at least 2 points"))
    n = N - 1 # number of FFT points
    x = [-cos(π*m/n) for m in 0:n] # Chebyshev nodes in ascending order
    v = [1 / (1 - 4*min(m, n-m)^2) for m in 0:n-1] # 1/(1-4k^2) and it's mirror image
    w = similar(x)
    w[begin:end-1] .= 2/n .* real(fft(v)) # Discrete Cosine Transform (DCT) for O(N log N) time instead of O(N²)
    w[begin] = w[end] = w[begin]/2 # modify endpoint factors
    return QuadratureRule(x, w, args...; name = Symbol("Clenshaw-Curtis"))
end

function GaussRule(N::Integer, args...)
    N ≥ 1 || throw(ArgumentError("The number of Gauss quadrature points must be positive"))
    x, w = QuadGK.gauss(N)
    return QuadratureRule(x, w, args...; name = Symbol("Gauss"))
end

function GaussKronrodRule(N::Integer, args...)
    N ≥ 3 && isodd(N) || throw(ArgumentError("The number of Gauss-Kronrod quadrature points must be at least 3 and odd"))
    x, w = QuadGK.kronrod(div(N-1, 2)) # only returns left half
    x = [x; -reverse(x[begin:end-1])] # add right half (antisymmetric)
    w = [w; reverse(w[begin:end-1])] # add right half (symmetric)
    return QuadratureRule(x, w, args...; name = Symbol("Gauss-Kronrod"))
end

transform(q::QuadratureRule, interval) = QuadratureRule(q.x, q.w, interval; name = q.name) # arrays are reuse and shared!
Base.eltype(::QuadratureRule{T, U}) where {T, U} = Base.promote_type(T, U)
Base.nameof(q::QuadratureRule) = q.name
Base.show(io::IO, q::QuadratureRule) = print(io, q.name == Symbol() ? "Q" : "$(q.name) q", "uadrature rule: on [$(q.a), $(q.b)], $(length(q)) points, eltype = $(eltype(q))")
Base.length(q::QuadratureRule) = length(q.x)
Base.eachindex(q::QuadratureRule) = eachindex(q.x)
Base.:(==)(q1::QuadratureRule, q2::QuadratureRule) = q1.x == q2.x && q1.w == q2.w && q1.a == q2.a && q1.b == q2.b
Base.:(≈)(q1::QuadratureRule, q2::QuadratureRule) = q1.x ≈ q2.x && q1.w ≈ q2.w && q1.a ≈ q2.a && q1.b ≈ q2.b
node(q::QuadratureRule, i) = (q.b+q.a)/2 + (q.b-q.a)/2 * q.x[i]
nodes(q::QuadratureRule) = (q.b+q.a)/2 .+ (q.b-q.a)/2 .* q.x
weights(q::QuadratureRule) = (q.b-q.a)/2 .* q.w # weights for integration on [a, b]
@inbounds @fastmath function integrate(q::QuadratureRule, f::AbstractVector)
    length(f) == length(q) || throw(ArgumentError("$(length(f))-point f is incompatible with $(length(q))-point quadrature rule"))
    return (q.b-q.a)/2 * sum(q.w[i] * f[i] for i in eachindex(q)) # integrate samples of f on [a, b]
end
@inbounds @fastmath function integrate(q::QuadratureRule, f::Function)
    return (q.b-q.a)/2 * sum(q.w[i] * f(node(q, i)) for i in eachindex(q)) # integrate f(x) on [a, b]
end
(q::QuadratureRule)(args...) = integrate(q, args...) # calling is equivalent to integrating
