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
