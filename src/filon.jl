using Bessels

struct SphericalBesselIntegralCache{T <: AbstractFloat}
    l::Vector{Int}
    x::Vector{T} # uniform grid
    w::Matrix{T} # jₗ(x) at this x-point
    I0::Matrix{T} # ∫dx jₗ(x) x⁰ cumulatively to this x-point
    I1::Matrix{T} # ∫dx jₗ(x) x¹ cumulatively to this x-point
    invdx::T # 1/dx (constant because grid is uniform)
    tmp1::Vector{T} # temporary workspace per l
    tmp2::Vector{T} # temporary workspace per l
    tmp3::Vector{T} # temporary workspace per l
    tmp4::Vector{T} # temporary workspace per l
end

# TODO: adaptively grid size from tolerance parameter? can compare to sin
function SphericalBesselIntegralCache(ls, x1 = 0.0, x2 = 10*ls[end]; mindx_grid = 2π/26, mindx_integral = 2π/2048, method = SimpsonEven())
    Nx = 1 + Int(ceil((x2 - x1) / mindx_grid))
    xs = range(x1, x2, length=Nx)
    dx = step(xs)
    xs = collect(xs)

    mindx_integral ≤ dx ≤ mindx_grid || error("Integration grid must be finer than interpolation grid")

    T = typeof(x1)
    I0s = zeros(T, length(ls), length(xs)) # TODO: add dummy for final+1 lookup?
    I1s = zeros(T, length(ls), length(xs))
    ws = zeros(T, length(ls), length(xs))

    Nsubx = 1 + Int(ceil((xs[2] - xs[1]) / mindx_integral))
    subis = zeros(T, Nsubx)

    @inbounds for il in eachindex(ls)
        l = ls[il]
        for ix in 1:length(xs)-1
            subxs = range(xs[ix], xs[ix+1]; length = Nsubx) # fine integration grid
            subis .= sphericalbesselj.(l, subxs) # jₗ(x) x⁰
            I0s[il, ix+1] = I0s[il, ix] + NumericalIntegration.integrate(subxs, subis, method) # cumulative ∫ jₗ(x) x⁰ dx
            subis .*= subxs # jₗ(x) x¹
            I1s[il, ix+1] = I1s[il, ix] + NumericalIntegration.integrate(subxs, subis, method) # cumulative ∫ jₗ(x) x¹ dx
        end
        ws[il, :] = sphericalbesselj.(l, xs)
    end

    tmp1 = zeros(T, length(ls))
    tmp2 = zeros(T, length(ls))
    tmp3 = zeros(T, length(ls))
    tmp4 = zeros(T, length(ls))
    return SphericalBesselIntegralCache{T}(ls, xs, ws, I0s, I1s, 1.0/dx, tmp1, tmp2, tmp3, tmp4)
end

function Base.show(io::IO, jlint::SphericalBesselIntegralCache)
    print(io, "jₗ(x) integration cache ")
    print(io, "for $(jlint.l[begin]) ≤ l ≤ $(jlint.l[end]) and ")
    print(io, "$(jlint.x[begin]) ≤ x ≤ $(jlint.x[end]) ")
    print(io, "($(Base.format_bytes(Base.summarysize(jlint))))\n")
end

# TODO: Hermite interpolation, have derivative for free: could use both (I0, w) and (I1, w*x)? https://en.wikipedia.org/wiki/Cubic_Hermite_spline
Base.@propagate_inbounds @fastmath function (jlint::SphericalBesselIntegralCache)(out0::AbstractArray{<:AbstractFloat}, out1::AbstractArray{<:AbstractFloat}, x)
    x1 = jlint.x[begin]
    ix₋ = 1 + unsafe_trunc(Int, (x - x1) * jlint.invdx) # O(1) uniform grid lookup (corresponding to left point) # TODO: floor(Int, instead?
    ix₊ = min(ix₋ + 1, length(jlint.x))
    x₋ = jlint.x[ix₋]
    x₊ = jlint.x[ix₊]
    dx = x₊ - x₋
    t = (x - x₋) * jlint.invdx
    h00 = (1 + 2t) * (1 - t)^2 # https://en.wikipedia.org/wiki/Cubic_Hermite_spline
    h01 = t^2 * (3 - 2t)
    h10 = t * (1 - t)^2 * dx # derivatives scaled by dx to map from x to t ∈ [0, 1]
    h11 = t^2 * (t - 1) * dx # derivatives scaled by dx to map from x to t ∈ [0, 1]
    @inbounds #=@simd=# for il in eachindex(out0)
        w₋ = jlint.w[il, ix₋]
        I0₋ = jlint.I0[il, ix₋]
        I1₋ = jlint.I1[il, ix₋]
        I0₋′ = w₋
        I1₋′ = w₋ * x₋
        w₊ = jlint.w[il, ix₊]
        I0₊ = jlint.I0[il, ix₊]
        I1₊ = jlint.I1[il, ix₊]
        I0₊′ = w₊
        I1₊′ = w₊ * x₊
        out0[il] = h00 * I0₋ + h10 * I0₋′ + h01 * I0₊ + h11 * I0₊′
        out1[il] = h00 * I1₋ + h10 * I1₋′ + h01 * I1₊ + h11 * I1₊′
    end
    return nothing
end

@fastmath function integrate(out, jlint::SphericalBesselIntegralCache, xs::AbstractArray, ys::AbstractArray)
    xs[2] > xs[1] || error("x-domain is not strictly increasing")
    xs[begin] ≥ jlint.x[begin] || error("$(xs[begin]) is outside integral cache left bound $(jlint.x[begin])")
    xs[end] ≤ jlint.x[end] || error("$(xs[end]) is outside integral cache right bound $(jlint.x[end])")
    out .= 0.0
    I0₋, I0₊, I1₋, I1₊ = jlint.tmp1, jlint.tmp2, jlint.tmp3, jlint.tmp4 # TODO: thread-unsafe
    i₋ = 1
    x₋ = xs[i₋]
    y₋ = ys[i₋]
    jlint(I0₋, I1₋, x₋)
    @inbounds while i₋ < length(xs) # TODO: write out and unroll loop to look up jlint(x) once per x and l without setindex?
        i₊ = i₋ + 1
        x₊ = xs[i₊]
        y₊ = ys[i₊]
        A = (y₊ - y₋) / (x₊ - x₋) # A in y = Ax + B
        B = y₋ - A*x₋ # B in y = Ax + B
        # TODO: quadratic; C?
        jlint(I0₊, I1₊, x₊) # TODO: write out this to avoid l-indexing?
        #=@simd=# for il in eachindex(out)
            I0 = I0₊[il] - I0₋[il] # ∫dx jₗ(x) x⁰ from x₋ to x₊
            I1 = I1₊[il] - I1₋[il] # ∫dx jₗ(x) x¹ from x₋ to x₊ # TODO: interpolate differences instead?
            out[il] += A*I1 + B*I0
        end
        I0₋, I1₋, I0₊, I1₊ = I0₊, I1₊, I0₋, I1₋ # pointer-like swap for next iteration
        i₋, x₋, y₋ = i₊, x₊, y₊ # pass to next iteration
    end
    return out
end

# Convenience dispatches
function integrate(out, jlint::SphericalBesselIntegralCache, f::Function, x1, x2; kw...)
    xs = range(x1, x2; kw...)
    ys = f.(xs)
    return integrate(out, jlint, xs, ys)
end
function integrate(out, jlint::SphericalBesselIntegralCache, f::Function, x2; kw...)
    x1 = jlint.x[begin]
    return integrate(out, jlint, f, x1, x2; kw...)
end
function integrate(out, jlint::SphericalBesselIntegralCache, f::Function; kw...)
    x2 = jlint.x[end]
    return integrate(out, jlint, f, x2; kw...)
end
function integrate(jlint::SphericalBesselIntegralCache, args...; kw...)
    out = zeros(eltype(jlint.I0), length(jlint.l))
    integrate(out, jlint, args...; kw...)
    return out
end