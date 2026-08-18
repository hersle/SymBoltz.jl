#=
Recurrence for spherical Bessel functions jₗ(k*χ) and hyperspherical Bessel functions Φₗ(χ, k)
inspired by arXiv:1312.2697 and arXiv:1311.0839 and CLASS' hyperspherical.c.

Forward/upward recursion is stable for χ above the turning point where k sinK(K,χ) = √(l(l+1)):
    1. seed Φ₀ = sin(kχ) / (k sinK(K,χ))
    2. seed Φ₁ = (cotK(K,χ) Φ₀ - cos(kχ)/sinK(K,χ)) / √K₁
    3. iterate √Kₗ Φₗ = (2l-1) cotK(K,χ) Φₗ₋₁ - √Kₗ₋₁ Φₗ₋₂ (where √Kₗ = √(k² - K l²))

Backward/downward recursion is stable for χ below the turning point:
    1. seed Φₗₘₐₓ₊₁/Φₗₘₐₓ from a continued fraction expansion
    2. iterate √Kₗ Φₗ₋₁ = (2l+1) cotK(K,χ) Φₗ - √Kₗ₊₁ Φₗ₊₁
    3. rescale all Φₗ depending on Φ₀
=#

"""
    jl_recurrence!(out, lmax, x)

Compute spherical Bessel functions ``jₗ(x)`` for every integer ``l = 0, \\dots, l_\\mathrm{max}`` and save them
in `out[l+1]`. Uses forward/backward recursion above/below the turning point ``√(lmax*(lmax+1))``.
"""
function jl_recurrence!(out::AbstractVector, lmax::Integer, x)
    if x == 0.0
        fill!(out, zero(eltype(x))) # jₗ(0) = 0 for l > 0 (the recursions would divide by 0)
        out[1] = one(eltype(x)) # j₀(0) = 1
    elseif x ≥ √(lmax*(lmax+1)) # "turning point": above it jₗ(x) is oscillatory, below it exponentially suppressed
        jl_forward!(out, lmax, x)
    else
        jl_backward!(out, lmax, x)
    end
    return out
end

function jl_forward!(out::AbstractVector, lmax::Integer, x)
    invx = 1 / x
    out[1] = sin(x) * invx # j₀(x) = sin(x) / x
    out[2] = (out[1] - cos(x)) * invx # j₁(x) = sin(x)/x² - cos(x)/x
    @inbounds for l in 1:lmax-1
        out[l+2] = (2l+1) * invx * out[l+1] - out[l] # jₗ₊₁(x) = (2l+1)/x jₗ(x) - jₗ₋₁(x)
    end
    return out
end

function jl_backward!(out::AbstractVector, lmax::Integer, x)
    invx = 1 / x # precompute; cheaper to multiply by 1/x instead of dividing by x
    jₗ = 2.0^-900 # arbitrary nonzero jₗₘₐₓ seed; renormalized at the end
    jₗ′_jₗ = jl_logderiv(lmax, x) # continued fraction of jₗₘₐₓ′/jₗₘₐₓ
    jₗ₊₁ = jₗ * (lmax*invx - jₗ′_jₗ) # jₗₘₐₓ₊₁ from derivative relation jₗ′ = l/x jₗ - jₗ₊₁
    out[lmax+1] = jₗ
    @inbounds for l in lmax:-1:1
        jₗ₋₁ = (2l+1) * invx * jₗ - jₗ₊₁
        if abs(jₗ₋₁) > 2.0^900 # renormalize previously computed jₗ before overflow
            s = 2.0^-900
            jₗ₋₁ *= s
            jₗ *= s
            @views out[l+1:lmax+1] .*= s
        end
        jₗ₊₁ = jₗ
        jₗ = jₗ₋₁
        out[l] = jₗ
    end
    s = (sin(x) * invx) / out[1]
    out .*= s # renormalize all jₗ to match analytical j₀
    return out
end

# Evaluate jₗ′(x)/jₗ(x) from its continued fraction expansion with the modified Lentz algorithm
function jl_logderiv(l::Integer, x)
    invx = 1 / x
    f = l * invx
    C = f
    D = zero(typeof(x))
    @inbounds for j in 1:1000 # increase max iterations if not converging
        b = (2*(l+j)+1) * invx
        D = 1 / (b - D)
        C = b - 1/C
        CD = C * D
        f *= CD # fₗ = jₗ′/jₗ
        abs(CD - 1) < eps(typeof(x)) && return f # converged
    end
    error("Continued fraction for jₗ′(x)/jₗ(x) did not converge at (l, x) = ($l, $x)")
end

"""
    Φl_recurrence!(out, lmax::Integer, χ, k, K, sqK::AbstractVector, invsqK::AbstractVector)

Compute hyperspherical Bessel functions ``Φₗ(χ, k)`` of a universe with curvature `K` for every integer
``l = 0, \\dots, l_\\mathrm{max}`` and store them in `out[l+1]`. Uses forward/backward recurrence
above/below the turning point, where ``k * sinK(K, χ) ≥ √(lmax*(lmax+1))``.
"""
function Φl_recurrence!(out::AbstractVector, lmax::Integer, χ, k, K, sqK::AbstractVector, invsqK::AbstractVector)
    if χ == 0.0
        fill!(out, zero(eltype(χ))) # Φₗ(0) = 0 for l > 0 (the recursions would divide by 0)
        out[1] = one(eltype(χ)) # Φ₀(0) = 1
    elseif k * sinK(K, χ) ≥ √(lmax*(lmax+1)) # above or below turning point?
        Φl_forward!(out, lmax, χ, k, K, sqK, invsqK)
    else
        Φl_backward!(out, lmax, χ, k, K, sqK, invsqK)
    end
    return out
end

"""
    Φl_recurrence!(out::AbstractMatrix, ls::AbstractArray{<:Integer}, χs::AbstractVector, k, K, sqK::AbstractVector, invsqK::AbstractVector; Φmin = 1e-20)

Same as above for every radial coordinate in `χs` at once, storing ``Φ_l(χ_i)`` in `out[i, l+1]` (so `out` must
be `length(χs)` × ``l_\\mathrm{max}+1``). `χs` must be sorted in descending order, as it is in the line-of-sight
integral, so that the parts above and below the turning point are contiguous views that can be handed to the
(vectorized) forward and backward recursion. Tables for `sqK` and `invsqK` are precomputed unless given explictly.
"""
function Φl_recurrence!(out::AbstractMatrix, ls::AbstractArray{<:Integer}, χs::AbstractVector, k, K, sqK::AbstractVector, invsqK::AbstractVector; Φmin = 1e-20)
    K ≤ 0 || throw(ArgumentError("K must be non-positive (flat or open universe), but got $K"))
    size(out) == (length(χs), length(ls)) || throw(ArgumentError("out size is $(size(out)), but should be $((length(χs), length(ls)))"))
    issorted(χs; rev = true) || throw(ArgumentError("χs is not sorted in descending order"))
    issorted(ls) || throw(ArgumentError("ls is not sorted in descending order"))
    lmax = ls[end]
    turningpoint = √(lmax*(lmax+1))
    ifwd = 1 : searchsortedlast(k .* sinK.(K, χs), turningpoint; rev = true) # LHS monotonically increasing for K < 0
    ibwd = (last(ifwd) + 1) : (iszero(χs[end]) ? length(χs) - 1 : length(χs))
    izero = (last(ibwd) + 1) : length(χs)
    if !isempty(ifwd)
        @views Φl_forward!(out[ifwd, :], ls, χs[ifwd], k, K, sqK, invsqK)
    end
    if !isempty(ibwd)
        @views Φl_backward!(out[ibwd, :], ls, χs[ibwd], k, K, sqK, invsqK; Φmin)
    end
    if !isempty(izero)
        @views out[izero, :] .= 0.0 # Φₗ(0) = 0 for l > 0 (the recursions would divide by 0)
        if ls[1] == 0
            @views out[izero, 1] .= 1.0 # Φ₀(0) = 1, only if l = 0 was requested
        end
    end
    return out
end

# Same, but without sqK and invsqK arrays
function Φl_recurrence!(out::AbstractMatrix, ls::AbstractArray{<:Integer}, χs::AbstractVector, k, K; kwargs...)
    K ≤ 0 || throw(ArgumentError("K must be non-positive (flat or open universe), but got $K"))
    lmax = ls[end]
    sqK = zeros(Float64, lmax+2) # TODO: don't require Float64
    invsqK = zeros(Float64, lmax+2)
    sqrtK_table!(sqK, invsqK, K, k)
    return Φl_recurrence!(out, ls, χs, k, K, sqK, invsqK; kwargs...)
end

function Φl_forward!(out::AbstractVector, lmax::Integer, χ, k, K, sqK::AbstractVector, invsqK::AbstractVector)
    ck = cotK(K, χ)
    sk = sinK(K, χ)
    out[1] = sin(k*χ) / (k * sk) # Φ₀
    out[2] = (ck * out[1] - cos(k*χ) / sk) * invsqK[2] # Φ₁ = (cotK Φ₀ - k cot(kχ) Φ₀) / √K₁
    @inbounds for l in 1:lmax-1
        out[l+2] = ((2l+1) * ck * out[l+1] - sqK[l+1] * out[l]) * invsqK[l+2] # sqK[l+1] = √Kₗ, invsqK[l+2] = 1/√Kₗ₊₁
    end
    return out
end

@inline function Φₗ_save_if_requested!(out, ls, l, Φ, lstoreidx)
    if l == ls[lstoreidx]
        @inbounds for i in 1:size(out, 1)
            out[i, lstoreidx] = Φ[i]
        end
        return lstoreidx + 1
    end
    return lstoreidx
end

# Same for many radial coordinates at once, with out[i, l+1] = Φₗ(χs[i]). The recursion in l is serial, but
# different χ are independent, so the χ-loop is innermost and vectorizes (hence χ is the contiguous dimension).
function Φl_forward!(out::AbstractMatrix, ls::AbstractArray{<:Integer}, χs::AbstractVector, k, K, sqK::AbstractVector, invsqK::AbstractVector)
    ck = cotK.(K, χs) # depends on χ but not on l, so tabulate it once instead of recomputing it every step
    Φₗ₋₁ = similar(χs)
    Φₗ = similar(χs)
    for i in eachindex(χs) # TODO: dot syntax
        sk = sinK(K, χs[i])
        Φₗ₋₁[i] = sin(k*χs[i]) / (k * sk) # Φ₀
        Φₗ[i] = (ck[i] * Φₗ₋₁[i] - cos(k*χs[i]) / sk) * invsqK[2] # Φ₁
    end

    lstoreidx = 1 # ls[lstoreidx] is the next multipole to output
    lstoreidx = Φₗ_save_if_requested!(out, ls, 0, Φₗ₋₁, lstoreidx) # save Φ₀?
    lstoreidx = Φₗ_save_if_requested!(out, ls, 1, Φₗ, lstoreidx) # save Φ₁?

    @inbounds for l in 1:ls[end]-1
        Φₗ₊₁ = Φₗ₋₁ # reuse array name (Φₗ₋₁[i] is not needed anymore writing Φₗ₊₁[i])
        @simd for i in eachindex(χs) # since it does not depend on χ
            Φₗ₊₁[i] = ((2l+1) * ck[i] * Φₗ[i] - sqK[l+1] * Φₗ₋₁[i]) * invsqK[l+2]
        end
        lstoreidx = Φₗ_save_if_requested!(out, ls, l+1, Φₗ₊₁, lstoreidx) # save Φₗ₊₁?
        Φₗ, Φₗ₋₁ = Φₗ₋₁, Φₗ # swap array pointers for next iteration
    end
    return out
end

function Φl_backward!(out::AbstractVector, lmax::Integer, χ, k, K, sqK::AbstractVector, invsqK::AbstractVector)
    ck = cotK(K, χ)
    Φₗ = 2.0^-900 # Φₗₘₐₓ, arbitrary nonzero seed (fixed by the final normalization)
    Φₗ′_Φₗ = Φl_logderiv(lmax, χ, k, K) # continued fraction of Φₗₘₐₓ′/Φₗₘₐₓ
    Φₗ₊₁ = Φₗ * (lmax*ck - Φₗ′_Φₗ) * invsqK[lmax+2] # Φₗₘₐₓ₊₁ from the derivative recursion Φₗ′ = l cotK Φₗ - √Kₗ₊₁ Φₗ₊₁
    out[lmax+1] = Φₗ
    @inbounds for l in lmax:-1:1
        Φₗ₋₁ = ((2l+1) * ck * Φₗ - sqK[l+2] * Φₗ₊₁) * invsqK[l+1] # sqK[l+2] = √Kₗ₊₁, invsqK[l+1] = 1/√Kₗ
        if abs(Φₗ₋₁) > 2.0^900 # renormalize previously computed Φₗ before overflow
            s = 2.0^-900
            Φₗ₋₁ *= s
            Φₗ *= s
            @views out[l+1:lmax+1] .*= s
        end
        Φₗ₊₁ = Φₗ
        Φₗ = Φₗ₋₁
        out[l] = Φₗ
    end
    s = sin(k*χ) / (k * sinK(K, χ)) / out[1]
    out .*= s
    return out
end

# Same for many radial coordinates at once, with out[i, l+1] = Φₗ(χs[i]); see Φl_forward! above. The recurrence
# and the overflow check are split into separate loops so that the recurrence can be @simd-vectorized over χ.
function Φl_backward!(out::AbstractMatrix, ls::AbstractArray{<:Integer}, χs::AbstractVector, k, K, sqK::AbstractVector, invsqK::AbstractVector; Φmin = 0.0)
    lmax = ls[end]
    _out = zeros(lmax + 1) # TODO: take as input argument?
    @inbounds for i in eachindex(χs)
        #println("χ = $(χs[i]), lmax = $lmax")
        Φl_backward!(_out, lmax, χs[i], k, K, sqK, invsqK)

        # Lower lmax to skip the recurrence Φₗ where it is below a tiny threshold
        # Immediately set Φₗ = 0 there
        # Works because Φₗ decreases monotonically as l decreases when χ is below the turning point
        while lmax ≥ 1 && abs(_out[1+lmax]) < Φmin
            _out[1+lmax] = 0.0
            lmax -= 1
        end

        for j in eachindex(ls)
            out[i, j] = _out[1 + ls[j]]
        end
    end
    return out
end

# Evaluate Φₗ'(χ; k)/Φₗ(χ; k) from its continued fraction expansion with the modified Lentz algorithm
function Φl_logderiv(l::Integer, χ, k, K)
    ck = cotK(K, χ)
    σ(m) = √(k^2 - K*m^2) # √Kₘ, computed on the fly since the fraction runs to arbitrarily high m
    nudge(x) = copysign(max(abs(x), 1e-100), x) # keep x not too close to zero
    f = nudge(l * ck) # b₀
    C = f
    D = zero(typeof(χ))
    σj = σ(l+1) # √K_{l+j} at j = 1
    @inbounds for j in 1:1000000 # increase max iterations if not converging
        σj1 = σ(l+j+1) # √K_{l+j+1}
        a = -σj / σj1
        j == 1 && (a *= σj) # the leading √K_{l+1} of the fraction, folded into the first numerator
        b = (2*(l+j)+1) * ck / σj1
        D = 1 / nudge(b + a*D)
        C = nudge(b + a/C)
        CD = C * D
        f *= CD # fₗ = Φₗ′/Φₗ
        abs(CD - 1) < eps(typeof(χ)) && return f # converged
        σj = σj1
    end
    error("Continued fraction for Φₗ′(χ)/Φₗ(χ) did not converge at (l, χ, k, K) = ($l, $χ, $k, $K)")
end

# Tabulate l-dependent coefficients ``√(Kₗ) = √(k² - K l²)`` of the hyperspherical Bessel recurrence and their inverses.
function sqrtK_table!(sqK::AbstractVector, invsqK::AbstractVector, K, k)
    @inbounds for i in eachindex(sqK, invsqK)
        l = i - 1
        Kₗ = k^2 - K*l^2
        sqK[i] = √(Kₗ)
        invsqK[i] = 1 / sqK[i]
    end
    return sqK, invsqK
end
