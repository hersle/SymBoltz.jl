"""
    massless_neutrinos(g; lmax = 10, interact = false, name = :ν, kwargs...)

Create a particle species for massless neutrinos in the spacetime with metric `g`.
"""
function massless_neutrinos(g; lmax = 10, interact = false, name = :ν, kwargs...)
    description = "Massless neutrinos"
    ν = radiation(g; adiabatic = true, name, description, kwargs...) |> background |> complete

    vars = @variables begin
        F0(τ, k), [description = "Distribution function monopole"]
        F(τ, k)[1:lmax], [description = "Distribution function multipoles"]
        δ(τ, k), [description = "Overdensity (gauge-dependent)"]
        Δ(τ, k), [description = "Overdensity (gauge-independent)"]
        δP(τ, k), [description = "Pressure perturbation"]
        θ(τ, k), [description = "Velocity divergence"]
        u(τ, k), [description = "Velocity"]
        σ(τ, k), [description = "Shear stress"]
        f(τ, k), [description = "Momentum transfer from other species"]
    end
    pars = @parameters begin
        N, [description = "Number of massless neutrino species"]
    end
    eqs = [
        D(F0) ~ -k*F[1] + 4*D(g.Φ)
        D(F[1]) ~ k/3*(F0-2*F[2]+4*g.Ψ) + 4/(3k) * f/(ν.ρ+ν.P) # (ρ+P)θ′ = … + f with θ = 3kF₁/4
        [D(F[l]) ~ k/(2*l+1) * (l*F[l-1] - (l+1)*F[l+1]) for l in 2:lmax-1]...
        D(F[lmax]) ~ k*F[lmax-1] - (lmax+1) / τ * F[lmax]
        δ ~ F0
        Δ ~ δ + 3*g.ℋ*(1+ν.w)*θ/k^2
        δP ~ ν.cₛ² * ν.ρ * δ
        θ ~ 3*k*F[1]/4
        σ ~ F[2]/2
        u ~ θ / k
    ]
    !interact && push!(eqs, f ~ 0)
    ieqs = [
        δ ~ -2 * g.Ψ # adiabatic: δᵢ/(1+wᵢ) == δⱼ/(1+wⱼ) (https://cmb.wintherscoming.no/theory_initial.php#adiabatic)
        θ ~ 1//2 * (k^2*τ) * g.Ψ
        σ ~ 1//15 * (k*τ)^2 * g.Ψ
        F[3] ~ 1//(2*3+1) * k*τ * F[2] # 1/(2l+1) * k*τ * F[l-1] → 0 quickly
        [F[l] ~ 0 for l in 4:lmax]...
    ]
    description = "Massless neutrinos"
    return extend(ν, System(eqs, τ, vars, pars; initialization_eqs = ieqs, name, kwargs...); description)
end

"""
    momentum_quadrature(f, N; x0 = 12, u = x -> 1/(1+x/x0), x = u -> x0*(1-u)/u, dx_du = u -> -x0/u^2, x1 = 0.0, x2 = Inf)

Compute ``N`` dimensionless momentum bins ``xᵢ`` and integral weights ``Wᵢ`` for integrating ``∫dx x²f(x)g(x)`` from ``0`` to ``∞``
against arbitrary weight functions ``g(x)`` with ``N``-point Gaussian quadrature using QuadGK.jl.
The returned weights `Ws` approximates the integral for arbitrary functions `g(x)` with the sum ``sum(Ws .* g.(xs))``.

The keyword arguments specifies an integral substitution ``x(u)`` with derivative ``\\mathrm{d}x/\\mathrm{d}u`` and inverse ``u(x)`` to apply.
The default transformation first maps ``x`` on the infinite domain ``(0, ∞)`` to ``x/x0``, where ``x0`` is a characteristic momentum of ``f``.
It then performs a rational transformation of ``x/x0`` into ``u`` on the finite domain ``(0, 1)`` to make the numerical integral well-defined.
As ``x0 → ∞``, this approaches Gaussian quadrature in ``x`` itself; for the Fermi-Dirac distribution, ``x0 ≈ 12`` is more accurate.
"""
function momentum_quadrature(f, N; x0 = 12, u = x -> 1/(1+x/x0), x = u -> x0*(1-u)/u, dx_du = u -> -x0/u^2, x1 = 0.0, x2 = Inf)
    w(x) = f(x) * x^2 # weight function to integrate against, i.e. want weights for ∫dx*w(x)*g(x) for arbitrary g(x)
    us, Ws = gauss(u -> dx_du(u) * x(u)^2 * f(x(u)), N, u(x1), u(x2)) # get u bins and quadrature weights
    xs = x.(us) # corresponding x values
    return xs, Ws
end

# Gρ(r) = ∫₀¹ds s²√(s²+r²) and GP(r) = ∫₀¹ds s⁴/√(s²+r²) and their derivatives,
# so ∫₀ˣdt t²√(t²+y²) = x⁴Gρ(y/x) and ∫₀ˣdt t⁴/√(t²+y²) = x⁴GP(y/x).
# For r ≫ 1 the closed forms cancel catastrophically, so sum the binomial series in 1/r² instead.
binomial_coeffs(a, p; n = 14) = ntuple(i -> Float64(prod((a - j) / (j + 1) for j in 0:i-2; init = 1//1) / (2(i-1) + p)), n)
const Gρ_series = binomial_coeffs(1//2, 3)
const Gρ′_series = binomial_coeffs(-1//2, 3)
const GP_series = binomial_coeffs(-1//2, 5)
const GP′_series = binomial_coeffs(-3//2, 5)
const GP″_series = 3 .* binomial_coeffs(-5//2, 5) .- GP′_series
Gρ(r) = r > 3 ? r * evalpoly(1/r^2, Gρ_series) : ((2 + r^2) * √(1 + r^2) - r^4 * asinh(1 / max(r, 1e-100))) / 8
Gρ′(r) = r > 3 ? evalpoly(1/r^2, Gρ′_series) : r * (√(1 + r^2) - r^2 * asinh(1 / max(r, 1e-100))) / 2
Gρ″(r) = r > 3 ? evalpoly(1/r^2, GP′_series) / r^3 : ((1 + 3r^2) / √(1 + r^2) - 3r^2 * asinh(1 / max(r, 1e-100))) / 2 # = -GP′(r)/r; needed when ℋ(ρ) is differentiated twice
GP(r) = r > 3 ? evalpoly(1/r^2, GP_series) / r : ((2 - 3r^2) * √(1 + r^2) + 3r^4 * asinh(1 / max(r, 1e-100))) / 8
GP′(r) = r > 3 ? -evalpoly(1/r^2, GP′_series) / r^2 : r * (3r^2 * asinh(1 / max(r, 1e-100)) - (1 + 3r^2) / √(1 + r^2)) / 2
GP″(r) = r > 3 ? evalpoly(1/r^2, GP″_series) / r^3 : (9r^2 * asinh(1 / max(r, 1e-100)) - (1 + 12r^2 + 9r^4) / (1 + r^2)^(3/2)) / 2
@register_symbolic Gρ(r)
@register_symbolic Gρ′(r)
@register_symbolic Gρ″(r)
@register_symbolic GP(r)
@register_symbolic GP′(r)
@register_symbolic GP″(r)
@register_derivative Gρ(r) 1 Gρ′(r)
@register_derivative Gρ′(r) 1 Gρ″(r)
@register_derivative GP(r) 1 GP′(r)
@register_derivative GP′(r) 1 GP″(r)

"""
    massive_neutrinos(g; nx = 4, x0 = 12, lmax = 10, name = :h, kwargs...)

Create a particle species for massive neutrinos in the spacetime with metric `g`.
The momentum integrals use `nx` quadrature points with the momentum scale `x0` (see [`momentum_quadrature`](@ref)).
The background density and pressure integrals are integrated by parts, so they use the same ``x²f₀ \\mathrm{d}\\ln f₀/\\mathrm{d}\\ln x`` kernel
as the perturbations. This makes the discrete background and adiabatic perturbations exactly consistent for any quadrature.
"""
function massive_neutrinos(g; nx = 4, x0 = 12, lmax = 10, name = :h, kwargs...)
    # compute numerical reduced momenta x = q*c / (kB*T) and Gaussian quadrature weights
    # for approximating integrals ∫dx x² f₀(x) g(x) for any g(x) over the infinite domain (0, ∞),
    # but change variables to transform it into a finite domain (0, 1)
    # (see e.g. https://juliamath.github.io/QuadGK.jl/v2.11/quadgk-examples/#Improper-integrals:-Infinite-limits)
    f₀(x) = 1 / (exp(x) + 1) # not exp(E); distribution function is "frozen in"; see e.g. Dodelson exercise 3.9
    dlnf₀_dlnx(x) = -x / (1 + exp(-x))
    x, W = momentum_quadrature(f₀, nx; x0)
    x² = x .^ 2
    dlnf₀ = dlnf₀_dlnx.(x)
    W′ = @. W * dlnf₀ / x # replaces f₀ by f₀′ in the weights, since f₀ dlnf₀/dlnx = xf₀′; integrate background by parts ∫dx f₀v′ = -∫dx f₀′v with v(0) = 0 to get the same f₀′ as the perturbation sources
    ∫dx_x²_f₀(g) = sum(collect(g .* W)) # ≈ ∫dx x²f₀(x)g(x)
    ∫dx_x²_f₀′(g) = sum(collect(g .* W′)) # ≈ ∫dx x²f₀′(x)g(x)

    pars = @parameters begin
        N = 3, [description = "Number of degenerate neutrino masses"]
        m_eV, [description = "Individual neutrino mass (in eV/c²)"] # TODO: only one m?
        m = m_eV * eV/c^2, [description = "Individual neutrino mass (in kg)"]
        T₀, [description = "Temperature today (in K)"]
        y₀ = m*c^2 / (kB*T₀), [description = "Temperature-reduced mass today"]
        Iρ₀ = -∫dx_x²_f₀′(@. x² * Gρ(y₀ / x)), [description = "Density integral today"] # by parts like Iρ; circumvent defining E₀[1:nx] because vector parameter dependencies doesn't work properly with setsym/remake
        Ω₀ = N * 8*Num(π)/3 * 2/(2*Num(π)^2) * (kB*T₀)^4 / (ħ*c)^3 * Iρ₀ / ((H100*g.h*c)^2/GN), [description = "Reduced background density today"]
    end
    vars = @variables begin
        ρ(τ), [description = "Background density"]
        P(τ), [description = "Background pressure"]
        Ω(τ), [description = "Reduced background density"]
        T(τ), [description = "Temperature"]
        y(τ), [description = "Temperature-deuced mass"]
        w(τ), [description = "Equation of state"]
        δP(τ, k), [description = "Pressure perturbation"]
        δ(τ, k), [description = "Overdensity (gauge-dependent)"]
        Δ(τ, k), [description = "Overdensity (gauge-independent)"]
        σ(τ, k), [description = "Shear stress"]
        θ(τ, k), [description = "Velocity divergence"]
        u(τ, k), [description = "Velocity"]
        E(τ)[1:nx], [description = "Dimensionless energies"]
        ψ0(τ, k)[1:nx], [description = "Distribution function monopole"]
        ψ(τ, k)[1:nx,1:lmax], [description = "Distribution function multipoles"]
        In(τ), [description = "Number density integral"]
        Iρ(τ), [description = "Density integral"]
        IP(τ), [description = "Pressure integral"]
        Iδρ(τ, k), [description = "Overdensity integral"]
    end

    eqs = [
        T ~ T₀ / g.a
        y ~ y₀ * g.a
        In ~ -∫dx_x²_f₀′(x / 3) # by parts with v′ = x², v = x³/3
        Iρ ~ -∫dx_x²_f₀′(@. x² * Gρ(y / x)) # by parts with v′ = x²E, v = x⁴Gρ(y/x)
        IP ~ -∫dx_x²_f₀′(@. x² * GP(y / x) / 3) # by parts with v′ = x⁴/(3E), v = x⁴GP(y/x)/3
        ρ ~ N/(π^2) * (kB*T)^4 / (ħ*c)^3 * Iρ / ((H100*g.h*c)^2/GN) # compute g/(2π²ħ³) * ∫dp p² √((pc)² + (mc²)²) / (exp(pc/(kT)) + 1) with dimensionless x = pc/(kT) and degeneracy factor g = 2
        P ~ N/(π^2) * (kB*T)^4 / (ħ*c)^3 * IP / ((H100*g.h*c)^2/GN) # compute g/(2π²ħ³) * ∫dp p² (pc)² / (3√((pc)² + (mc²)²)) / (exp(pc/(kT)) + 1) with dimensionless x = pc/(kT) and degeneracy factor g = 2
        w ~ P / ρ
        Ω ~ 8*Num(π)/3 * ρ

        Iδρ ~ ∫dx_x²_f₀(E .* ψ0)
        δ ~ Iδρ / Iρ
        Δ ~ δ + 3*g.ℋ*(1+w)*θ/k^2
        u ~ ∫dx_x²_f₀(x .* ψ[:,1]) / (Iρ + IP)
        θ ~ u * k
        σ ~ (2//3) * ∫dx_x²_f₀(x² ./ E .* ψ[:,2]) / (Iρ + IP)
        δP ~ P * ∫dx_x²_f₀(x² ./ (3E) .* ψ0) / IP

        [E[i] ~ √(x[i]^2 + y^2) for i in 1:nx]...
        [D(ψ0[i]) ~ -k * x[i]/E[i] * ψ[i,1] - D(g.Φ) * dlnf₀[i] for i in 1:nx]...
        [D(ψ[i,1]) ~ k/3 * x[i]/E[i] * (ψ0[i] - 2*ψ[i,2]) - k/3 * E[i]/x[i] * g.Ψ * dlnf₀[i] for i in 1:nx]...
        [D(ψ[i,l]) ~ k/(2*l+1) * x[i]/E[i] * (l*ψ[i,l-1] - (l+1) * ψ[i,l+1]) for i in 1:nx, l in 2:lmax-1]...
        [D(ψ[i,lmax]) ~ k/(2*lmax+1) * x[i]/E[i] * (lmax*ψ[i,lmax-1] - (lmax+1) * ((2*lmax+1) * E[i]/x[i] * ψ[i,lmax] / (k*τ) - ψ[i,lmax-1])) for i in 1:nx]... # explicitly inserted ψ[lmax+1] to avoid array allocations in newer MTK (see example in https://github.com/SciML/ModelingToolkit.jl/issues/3708)
    ]
    ieqs = [
        [ψ0[i] ~ -1//4 * (-2*g.Ψ) * dlnf₀[i] for i in 1:nx]...
        [ψ[i,1] ~ -1//3 * E[i]/x[i] * (1/2*k*τ*g.Ψ) * dlnf₀[i] for i in 1:nx]...
        [ψ[i,2] ~ -1//2 * (1//15*(k*τ)^2*g.Ψ) * dlnf₀[i] for i in 1:nx]...
        [ψ[i,l] ~ 0 for i in 1:nx, l in 3:lmax]... # TODO: full ICs
    ]
    description = "Massive neutrino"
    pars = [m, m_eV, N, T₀, Ω₀, y₀, T₀, Iρ₀]
    return System(eqs, τ, vars, pars; initialization_eqs = ieqs, name, description, kwargs...)
end
