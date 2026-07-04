struct Quadrature{T, U}
    x::Vector{T} # integration points on [-1, +1]
    w::Vector{T} # integration weights on [-1, +1]
    a::U
    b::U
    name::Symbol

    function Quadrature(x::AbstractVector{T}, w::AbstractVector{T}, interval = (-1.0, 1.0); name = Symbol()) where {T}
        all(-1 .≤ x .≤ 1) || throw(ArgumentError("Quadrature points must be on the canonical interval [-1, 1]"))
        sum(w) ≈ 2 || throw(ArgumentError("Quadrature weights must sum to 2, but sums to $(sum(w))"))
        issorted(x) || throw(ArgumentError("Quadrature points must be sorted in ascending order"))
        length(x) == length(w) || throw(ArgumentError("Quadrature rule must have same number of nodes and weights"))
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

transform(q::Quadrature, interval) = Quadrature(q.x, q.w, interval; name = q.name) # arrays are reuse and shared!
Base.eltype(::Quadrature{T, U}) where {T, U} = Base.promote_type(T, U)
Base.nameof(q::Quadrature) = q.name
Base.show(io::IO, q::Quadrature) = print(io, q.name == Symbol() ? "Q" : "$(q.name) q", "uadrature rule: on [$(q.a), $(q.b)], $(length(q)) points, eltype = $(eltype(q))")
Base.length(q::Quadrature) = length(q.x)
Base.eachindex(q::Quadrature) = eachindex(q.x)
Base.:(==)(q1::Quadrature, q2::Quadrature) = q1.x == q2.x && q1.w == q2.w && q1.a == q2.a && q1.b == q2.b
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
