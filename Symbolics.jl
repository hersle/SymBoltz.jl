using ModelingToolkit
using Symbolics
using LinearAlgebra: diagm # TODO: avoid big dependency?
using Tullio: @einsum # Einstein summation convention # TODO: use Tensorial.@einsum? # TODO: import as einsum
using Base.Iterators
using DifferentialEquations
using Plots

# TODO: use Grassmann.jl? Tensorial.jl? UnitfulBuckinghamPi.jl?
# TODO: solve analytical differential equations symbolically (not possible for now: https://discourse.julialang.org/t/solve-differential-equation-algebraically/106352/6; instead e.g. propose ansatz a^n and solve for n?)
# TODO: write parametrized Tensor type? index T[+1] for upper indices and T[-1] for lower indices?
# TODO: specify tensor symmetries with Tensorial?
# TODO: Spacetime.jl/Metric.jl file with different metrics in different gauges?

const MAXORDER = 1
@variables ϵ # perturbation book-keeping parameter
pertorder(expr, ϵ, n) = substitute(n == 0 ? expr : expand_derivatives((Differential(ϵ^n))(expr)), ϵ => 0; fold=false) 
pertseries(expr, ϵ, n) = sum(pertorder.(expr, ϵ, 0:n) .* ϵ .^ (0:n))
perttrunc(expr) = pertseries(expr, ϵ, MAXORDER) # TODO: // 1 here to enforce rational/fractions?
# TODO: make some kind of O(ϵ^0), O(ϵ^1), ... functions?

# cancel a factor if it multiplies a *whole* expression
# TODO: fix. sometimes it returns Num, sometimes BasicSymbolic, ...
function cancel_factor(expr, factor)
    expr_reduced = expand_derivatives(Differential(factor)(expr))
    return simplify(expand(expr_reduced * factor - expr)) === Num(0) ? expr_reduced : expr
end

function cancel_factors(expr, factors)
    for factor in factors
        expr = cancel_factor(expr, factor)
    end
    return expr
end

@variables t, x, y, z # coordinates # TODO: what about parametrizing with b = ln(a)?
@variables a(t) # background (O(ϵ⁰))
@variables Φ(a,x,y,z), Ψ(a,x,y,z) # TODO: linear perturbations (O(ϵ¹))

X = [t, x, y, z]
∂ = expand_derivatives .∘ Differential.(X)

# calculate metric g_μν -> inverse metric g^μν -> connection Γ^σ_μν -> Ricci tensor Ric_μν -> Ricci scalar R -> Einstein tensor G_μν
g = diagm([-1-2*ϵ*Ψ, a^2*(1+2*ϵ*Φ), a^2*(1+2*ϵ*Φ), a^2*(1+2*ϵ*Φ)]) .|> perttrunc # g_μν
ginv = inv(g) .|> perttrunc # g^μν (give .|> here, but |> with @einsum below)
@einsum Γ[σ,μ,ν] := ginv[σ,ρ] / 2 * (∂[μ](g[ν,ρ]) + ∂[ν](g[ρ,μ]) - ∂[ρ](g[μ,ν])) |> perttrunc |> simplify # Γ^σ_μν
@einsum Ric[μ,ν] := (∂[ρ](Γ[ρ,ν,μ]) - ∂[ν](Γ[ρ,ρ,μ])) / 4 + Γ[ρ,ρ,λ]*Γ[λ,ν,μ] - Γ[ρ,ν,λ]*Γ[λ,ρ,μ] |> perttrunc |> simplify # Ric_μν # /4 since @einsum overcounts first two terms due to λ-summation in last two terms
R = (@einsum _ := ginv[μ,ν] * Ric[μ,ν]) |> perttrunc |> simplify # TODO: why does one term have float?
G = Ric .- g/2 * R .|> perttrunc .|> expand .|> simplify # G_μν # TODO: always do perttrunc > expand > simplify_fractions

# function ∇(T::Tensor{S, Num, R, N}) where {S, R, N} # TODO: specialize at compile time?
function ∇(T) # assume T = T^α1...αr has only upper indices for now # TODO: generalize
    ∇T = zeros(Num, (length(∂), size(T)...))
    for μ in eachindex(∂)
        for α in CartesianIndices(T)
            ∇T[μ,α.I...] = ∂[μ](T[α.I...])
            for r in 1:ndims(T) # rank
                for d in 1:size(T, r)
                    ∇T[μ,α.I...] += Γ[α[r],d,μ] * T[α.I[begin:r-1]..., d, α.I[r+1:end]...]
                end
            end
        end
    end
    ∇T # Tensor{Tuple{∇Tsize...}}(∇T)
end

# TODO: type-based Species/Ingredient/Stuff
function construct_species(w, cs2, subscript)
    # energy density ρ
    ρ0 = Symbolics.variable("ρ$(subscript)0"; T=Symbolics.FnType)(a) # TODO: nicer way to create?
    δ  = Symbolics.variable("δ$(subscript)"; T=Symbolics.FnType)(a, x, y, z)
    ρ1 = δ*ρ0 # TODO: use this as the main variable?
    ρ = ρ0 + ϵ*ρ1 # TODO: more systematic expansion, e.g. create array of [ρ0, ρ1, ρ2, ...]?

    # pressure P
    P = w*ρ0 + ϵ*cs2*ρ1

    # velocity v or u
    vx = Symbolics.variable("v$(subscript)x"; T=Symbolics.FnType)(a, x, y, z)
    vy = Symbolics.variable("v$(subscript)y"; T=Symbolics.FnType)(a, x, y, z)
    vz = Symbolics.variable("v$(subscript)z"; T=Symbolics.FnType)(a, x, y, z)
    u = [0, ϵ*vx/a, ϵ*vy/a, ϵ*vz/a]
    u2 = (@einsum _ := g[μ,ν] * u[μ] * u[ν]) |> perttrunc |> simplify # == g[i,j] * u[i] * u[j] because u[0] == 0
    u = [√(-(1+u2)/g[1,1]), u[2], u[3], u[4]] .|> perttrunc |> simplify # set u[1] from normalization u^2 == -1 of u[i] (u[0] = 0 )
    #u = substitute.(u, 1.0 => 1) # hack

    # energy-momentum tensor T
    @einsum T[μ,ν] := (ρ+P)*u[μ]*u[ν] + P*ginv[μ,ν] # (with high indices)
    T = simplify.(perttrunc.(T))

    ρ, δ, P, vx, vy, vz, T
end

ρr, δr, Pr, vrx, vry, vrz, Tr = construct_species(1/3, 1/3, "r") # TODO: use symbol instead of string?
ρm, δm, Pm, vmx, vmy, vmz, Tm = construct_species(0/1, 0/1, "m")
# TODO: treat cosmological constant as species with ω = -1?
T = Tr + Tm # T^μν
@einsum Tlo[μ,ν] := g[μ,α] * g[ν,β] * T[α,β] # |> perttrunc # T_μν

# ∇_μ T^μν = 0
@einsum tr[ν] := ∇(Tr)[μ,μ,ν] |> perttrunc
@einsum tm[ν] := ∇(Tm)[μ,μ,ν] |> perttrunc

# TODO: ∇_μ (n u^μ) = 0?

# cosmological constant # TODO: treat as species instead?
@variables Λ

# define "full" equations of motion
eoms = [
    # background equations
    pertorder(tr[1], ϵ, 0), # radiation density evolution
    pertorder(tm[1], ϵ, 0), # matter density evolution
    pertorder(G[1,1] + Λ*g[1,1] - 8*π*Tlo[1,1], ϵ, 0), # Friedmann equation

    # TODO: perturbations
    pertorder(tr[1], ϵ, 1), # radiation density evolution
    pertorder(tm[1], ϵ, 1), # matter density evolution
    pertorder(tr[2], ϵ, 1), # radiation velocity evolution
    pertorder(tm[2], ϵ, 1), # matter velocity evolution
    pertorder(G[1,1] + Λ*g[1,1] - 8*π*Tlo[1,1], ϵ, 1), # perturbed Friedmann equation
    pertorder(G[2,3] + Λ*g[2,3] - 8*π*Tlo[2,3], ϵ, 1), # shear stress
] .~ 0

# introduce Hubble parameter H(a) = H0*E(a) and reduced density parameters
@variables E(a) Ωr(a) Ωm(a)
@parameters ΩΛ
@parameters H0 # TODO: turn into normal variables? seems parameters aren't "handled" by differentiation
H = H0 * E
ρcrit = 3/8 * H0^2 / π
eoms = substitute(eoms, Dict(
    pertorder(ρm, ϵ, 0) => Ωm * ρcrit, # TODO: would like ρm[0] => ... instead
    pertorder(ρr, ϵ, 0) => Ωr * ρcrit,
    Λ => ΩΛ * 3*H0^2,
    ∂[1](a) => H*a, # TODO: would need to apply one more time for acceleration equation
), fold=false)

# re-parametrize from t to a
aindep = Symbolics.variable("a")
eoms = substitute.(eoms, a => aindep; fold=false)
Ωr, Ωm, E, Ψ, Φ, δr, δm, vrx, vry, vrz, vmx, vmy, vmz = substitute.([Ωr, Ωm, E, Ψ, Φ, δr, δm, vrx, vry, vrz, vmx, vmy, vmz], a => aindep) # TODO: handle all this more automatically
a = aindep # forget old a(t)

# manipulate background
# take square of E^2 # TODO: improve
eoms[3] = E ~ √(abs(expand(simplify_fractions(Symbolics.solve_for(eoms[3], E^2))))) # TODO: get rid of abs!!!

# TODO: use simplify_fractions instead of simplify to avoid floats!
# TODO: and make it work on equations and not only expressions!
for (i, eom) in enumerate(eoms)
    println(i, " -> ", eom)
end

# Fourier transform perturbations
# TODO: write e.g. δr(a,r) instead of δr(a,x,y,z)?
@parameters i # hacky √(-1), will substitute later # TODO: handle properly!!
@parameters kx ky kz k K
@variables Ψk(a) Φk(a) δrk(a) δmk(a) vrk(a) vmk(a)
kdotx = kx*x + ky*y + kz*z
eoms = expand_derivatives.(substitute(eoms, Dict(
    Ψ => Ψk * exp(i*kdotx),
    Φ => Φk * exp(i*kdotx),
    δr => δrk * exp(i*kdotx),
    δm => δmk * exp(i*kdotx),

    # assume ∇ × v == 0 ⟹ i*k × vk == 0 ⟹ v ∥ k ⟹ vxk == kx/abs(k) * vk
    # also define *real* Fourier coefficients by multiplying vk -> i*vk by i
    vrx => i * kx/k * vrk * exp(i*kdotx),
    vry => i * ky/k * vrk * exp(i*kdotx),
    vrz => i * kz/k * vrk * exp(i*kdotx),
    vmx => i * kx/k * vmk * exp(i*kdotx),
    vmy => i * ky/k * vmk * exp(i*kdotx),
    vmz => i * kz/k * vmk * exp(i*kdotx),
); fold=false))

eoms = substitute.(eoms, Ref(Dict(
    i^2 => -1, # hacky √(-1)
    exp(i*kdotx) => 1, # TODO: avoid?
)); fold=false)
eoms = substitute.(eoms, kx^2 => k^2 - ky^2 - kz^2; fold=false)
eoms = substitute.(eoms, k => K*H0; fold=false)
eoms = substitute.(eoms, Ref(Dict(kx => 1, ky => 1, kz => 1)); fold=false) # TODO: avoid
eoms = substitute.(eoms, i => 1; fold=false) # TODO: avoid
eoms = substitute.(eoms, H0 => 1; fold=false) # TODO: avoid
eoms[9] = Φk ~ Symbolics.solve_for(eoms[9], Φk) # TODO: avoid
eoms = simplify.(eoms)

# eliminate Φ # TODO: avoid
eoms = expand_derivatives.(substitute.(eoms, Φk => eoms[9].rhs))
eoms = expand.(eoms[1:end-1])

# solve
function solve_for_wavenumber(Kval)
    @named sys0 = ODESystem(eoms, a, [Ωr, Ωm, δrk, δmk, vrk, vmk, Ψk], [ΩΛ, K])
    sys = structural_simplify(sys0; simplify=true, allow_symbolic=true, allow_parameter=true, check_consistency=false)
    aini = 1e-7
    Ωr0 = 1e-4
    Ωm0 = 0.3
    ΩΛ0 = 1.0 - Ωr0 - Ωm0
    Ψ0 = -1.0
    y0 = [sys.Ωr => Ωr0 / aini^4, sys.Ωm => Ωm0 / aini^3, sys.ΩΛ => ΩΛ0, sys.Ψk => Ψ0, sys.δmk => -3/2*Ψ0, sys.δrk => -3/2*Ψ0, sys.vmk => -K/(2*√(Ωr0))*Ψ0, sys.vrk => +K/(6*√(Ωr0))*Ψ0, K => Kval] # TODO: set ICs "automatically"
    prob = ODEProblem(sys, y0, (aini, 1.0))
    return solve(prob)
end

sol = solve_for_wavenumber(1.0)

p = plot()
plot!(p, log10.(sol[a]), sol[Ωm] ./ sol[E] .^ 2) # TODO: add observed equations for these?
plot!(p, log10.(sol[a]), sol[Ωr] ./ sol[E] .^ 2)
plot!(p, log10.(sol[a]), sol[Ψk])
display(p)

Ks = 10 .^ range(-2, +2; length=20)
sols = solve_for_wavenumber.(Ks)
Ps = [sol[δmk][end]^2 for sol in sols]

p = plot()
plot!(p, log10.(Ks), log10.(Ps))

# TODO: set background ICs today?

#=
# 9) Boltzmann equations (see Baumann chapter 3)
# TODO: integration: https://docs.sciml.ai/SymbolicNumericIntegration/stable/
#@variables p
#E = p^2 # TODO: mass
#f = 1 / (exp(E/T) - 1) # boson distribution function (+ 1 for fermions)
# TODO: go to Boltzmann equation?
#ρ = g/(2*Num(π))^3 * integral(4*π*p^2 * f(p))
# TODO: start with CDM only?
=#