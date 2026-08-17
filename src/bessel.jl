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
    Φl_recurrence!(out::AbstractMatrix, lmax, χs::AbstractVector, k, K, sqK, invsqK)

Same as above for every radial coordinate in `χs` at once, storing ``Φ_l(χ_i)`` in `out[i, l+1]` (so `out` must
be `length(χs)` × ``l_\\mathrm{max}+1``). `χs` must be sorted in descending order, as it is in the line-of-sight
integral, so that the parts above and below the turning point are contiguous views that can be handed to the
(vectorized) forward and backward recursion.
"""
function Φl_recurrence!(out::AbstractMatrix, lmax::Integer, χs::AbstractVector, k, K, sqK::AbstractVector, invsqK::AbstractVector)
    @assert size(out) == (length(χs), lmax+1) "out is $(size(out)), but should be $((length(χs), lmax+1))"
    @assert issorted(χs; rev = true) "χs must be sorted in descending order"
    ifwd = count(χ -> k * sinK(K, χ) ≥ √(lmax*(lmax+1)), χs) # χs[1:ifwd] are above the turning point,
    izero = count(!=(0.0), χs) # and χs[izero+1:end] are exactly zero
    @views Φl_forward!(out[1:ifwd, :], lmax, χs[1:ifwd], k, K, sqK, invsqK)
    @views Φl_backward!(out[ifwd+1:izero, :], lmax, χs[ifwd+1:izero], k, K, sqK, invsqK)
    @views out[izero+1:end, :] .= 0.0 # Φₗ(0) = 0 for l > 0 (the recursions would divide by 0)
    @views out[izero+1:end, 1] .= 1.0 # Φ₀(0) = 1
    return out
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

# Same for many radial coordinates at once, with out[i, l+1] = Φₗ(χs[i]). The recursion in l is serial, but
# different χ are independent, so the χ-loop is innermost and vectorizes (hence χ is the contiguous dimension).
function Φl_forward!(out::AbstractMatrix, lmax::Integer, χs::AbstractVector, k, K, sqK::AbstractVector, invsqK::AbstractVector)
    ck = cotK.(K, χs) # depends on χ but not on l, so tabulate it once instead of recomputing it every step
    @inbounds for i in eachindex(χs)
        sk = sinK(K, χs[i])
        out[i, 1] = sin(k*χs[i]) / (k * sk) # Φ₀
        out[i, 2] = (ck[i] * out[i, 1] - cos(k*χs[i]) / sk) * invsqK[2] # Φ₁ = (cotK Φ₀ - k cot(kχ) Φ₀) / √K₁
    end
    @inbounds for l in 1:lmax-1
        a, b, c = 2l+1, sqK[l+1], invsqK[l+2] # sqK[l+1] = √Kₗ, invsqK[l+2] = 1/√Kₗ₊₁; hoisted out of the χ-loop
        @simd for i in eachindex(χs) # since it does not depend on χ
            out[i, l+2] = (a * ck[i] * out[i, l+1] - b * out[i, l]) * c
        end
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

# Same for many radial coordinates at once, with out[i, l+1] = Φₗ(χs[i]); see Φl_forward! above. The overflow
# check and rescaling below happen exactly as in the single-χ method, just once per χ instead of once overall.
function Φl_backward!(out::AbstractMatrix, lmax::Integer, χs::AbstractVector, k, K, sqK::AbstractVector, invsqK::AbstractVector)
    ck = cotK.(K, χs)
    Φₗ = similar(ck)
    Φₗ₊₁ = similar(ck)
    @inbounds for i in eachindex(χs)
        Φₗ[i] = 2.0^-900 # Φₗₘₐₓ, arbitrary nonzero seed (fixed by the final normalization)
        Φₗ′_Φₗ = Φl_logderiv(lmax, χs[i], k, K) # continued fraction of Φₗₘₐₓ′/Φₗₘₐₓ
        Φₗ₊₁[i] = Φₗ[i] * (lmax*ck[i] - Φₗ′_Φₗ) * invsqK[lmax+2] # Φₗₘₐₓ₊₁ from the derivative recursion Φₗ′ = l cotK Φₗ - √Kₗ₊₁ Φₗ₊₁
        out[i, lmax+1] = Φₗ[i]
    end
    @inbounds for l in lmax:-1:1
        a, b, c = 2l+1, sqK[l+2], invsqK[l+1] # sqK[l+2] = √Kₗ₊₁, invsqK[l+1] = 1/√Kₗ; hoisted, see Φl_forward!
        for i in eachindex(χs)
            Φₗ₋₁ = (a * ck[i] * Φₗ[i] - b * Φₗ₊₁[i]) * c
            if abs(Φₗ₋₁) > 2.0^900 # renormalize previously computed Φₗ (for this χ) before overflow
                s = 2.0^-900
                Φₗ₋₁ *= s
                Φₗ[i] *= s
                for m in l+1:lmax+1
                    out[i, m] *= s
                end
            end
            Φₗ₊₁[i] = Φₗ₋₁ # old Φₗ₊₁[i] is no longer needed (already used above); overwrite it with the new value
            out[i, l] = Φₗ₋₁
        end
        Φₗ, Φₗ₊₁ = Φₗ₊₁, Φₗ # pointer swap instead of copying every element: what now holds Φₗ₋₁ becomes Φₗ for the next l
    end
    s = Φₗ # reuse as dummy array
    @inbounds for i in eachindex(χs)
        s[i] = sin(k*χs[i]) / (k * sinK(K, χs[i])) / out[i, 1]
    end
    @inbounds for m in 1:lmax+1
        @simd for i in eachindex(χs)
            out[i, m] *= s[i]
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
