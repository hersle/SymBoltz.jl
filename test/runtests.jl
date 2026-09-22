using Test
using SymBoltz
using ModelingToolkit
using ForwardDiff
using FiniteDiff
using BenchmarkTools
using Base.Threads
using Statistics
using DelimitedFiles
using StaticArrays
using SciMLBase

lmax = 5
M = ΛCDM(K = nothing; lmax) # flat
pars = parameters_Planck18(M)
prob = CosmologyProblem(M, pars)
prob_dense = CosmologyProblem(M, pars; bgsparse = false, ptsparse = false)
prob_sparse = prob

τ, k, D = SymBoltz.τ, SymBoltz.k, SymBoltz.D

# Must come first because warnings are only given once
@testset "Solve failure warnings" begin
    Ωc0 = prob.bg[1].ps[M.c.Ω₀]
    prob.bg[1].ps[M.c.Ω₀] = NaN # bad
    bgsol = @test_warn "Background solution failed" solvebg(prob.bg[1])
    prob.bg[1].ps[M.c.Ω₀] = Ωc0 # restore good
    bgsol = @test_nowarn solvebg(prob.bg[1])
    bgsols = solvebg(prob)

    @test_warn "Perturbation (mode k = NaN) solution failed" ptsol = solvept(prob.pt, bgsols, [NaN]; thread = false)
    @test_nowarn ptsol = solvept(prob.pt, bgsols, [1.0]; thread = false)
end

@testset "Base.show" begin
    # print to a buffer to keep CI output clean, but still catch errors and warnings
    @test_nowarn sprint(show, prob)
    @test_nowarn sprint(show, prob_dense)
    @test_nowarn sprint(show, prob_sparse)
end

@testset "Solution accessing" begin
    is = [M.g.a, M.g.a, M.g.a, M.g.a]
    τs = [1.0, 2.0, 3.0]
    ks = [1e2, 1e3]
    ni, nt, nk = length(is), length(τs), length(ks)

    # Size of solution output should match input arguments in order
    sol = solve(prob)
    @test size(sol(is[1], τs[1], )) == () # background
    @test size(sol(is, τs[1])) == (ni,)
    @test size(sol(is[1], τs)) == (nt,)
    @test size(sol(is, τs)) == (ni, nt)
    @test_throws "No perturbations" sol(ks, τs, is)
    #@test_throws "below minimum solved time" sol(sol[M.τ][begin]-1, is)
    #@test_throws "above maximum solved time" sol(sol[M.τ][end]+1, is)

    sol = solve(prob, ks)
    @test size(sol(is[1], τs[1], ks[1])) == () # perturbations
    @test size(sol(is, τs[1], ks[1])) == (ni,)
    @test size(sol(is[1], τs, ks[1])) == (nt,)
    @test size(sol(is, τs, ks[1])) == (ni, nt)
    @test size(sol(is[1], τs[1], ks)) == (nk,)
    @test size(sol(is, τs[1], ks)) == (ni, nk)
    @test size(sol(is[1], τs, ks)) == (nt, nk)
    @test size(sol(is, τs, ks)) == (ni, nt, nk)
    @test_throws "below minimum solved wavenumber" sol(is, τs, ks[begin]-1)
    @test_throws "above maximum solved wavenumber" sol(is, τs, ks[end]+1)
    #@test_throws "below minimum solved time" sol(ks[1], sol[M.τ][begin]-1, is)
    #@test_throws "above maximum solved time" sol(ks[1], sol[M.τ][end]+1, is)

    # TODO: also test array indexing
end

@testset "Backwards integration from today" begin
    sol = solve(prob)
    @test length(sol.bg) == 2 && sol.bg[end].t[end] < sol.bg[end].t[begin] # χ and κ are integrated backwards in the last stage
    @test sol.bg[end].t[begin] == sol.bg[1].t[end] # from where the first stage terminates
    @test sol[M.χ][end] == 0.0
    @test sol[M.b.κ][end] == 0.0
    @test sol[M.χ] ≈ sol[M.τ][end] .- sol[M.τ]
end

@testset "Accessing derivative variables" begin
    ks = 1e3
    sol = solve(prob, ks)
    τ0 = sol[M.τ][end]

    # derivatives are not interpolated; add them to the model as equations like dx ~ D(x) instead
    @test_throws "not present in the system" sol(D(M.g.a), τ0)
    @test_throws "not present in the system" sol(D(M.g.Φ), τ0, ks)
end

@testset "Solution interpolation" begin
    ks = 10 .^ range(-2, 4, length=100)
    sol = solve(prob, ks)
    ks = range(extrema(ks)..., length=500)
    τs = range(extrema(sol[M.τ])..., length=500)
    is = [M.g.a, M.G.ρ, M.g.Φ, M.g.Ψ]
    @test sol(is, τs, ks; smart = true) == sol(is, τs, ks; smart = false)
end

@testset "Spherical Bessel function" begin
    l = 0:1000
    jlfast = zeros(size(l))
    jlslow = zeros(size(l))

    for x in 0.0:0.001:10.0 # near 0 is most sketchy
        # Test jₗ(x)
        SymBoltz.jl!(jlfast, l, x) # "unsafe" implementation
        SymBoltz.jlsafe!(jlslow, l, x) # safe implementation
        @test jlfast ≈ jlslow
    end
end

@testset "Spherical Bessel function chain rule" begin
    x = 0.0:0.1:10.0

    # Test jl(l, x) chain rule
    crazy(l, x) = sin(7*SymBoltz.jl(l, x^2)) # crazy composite function involving jl
    for l in 1:500
        dcrazy_fd(l, x) = FiniteDiff.finite_difference_derivative(x -> crazy(l, x), x)
        dcrazy_ad(l, x) = ForwardDiff.derivative(x -> crazy(l, x), x)
        @test all(isapprox.(dcrazy_ad.(l, x), dcrazy_fd.(l, x); atol = 1e-6))
    end
end

@testset "Spherical Bessel function cache" begin
    ls = 10:10:100
    i5 = 0
    i10 = 1
    jl_lin = SphericalBesselCache(ls; dx = 2π/150, hermite = false)
    jl_her = SphericalBesselCache(ls; dx = 2π/15, hermite = true)
    for (jl, atol) in [(jl_lin, 1e-5), (jl_her, 1e-5)]
        @test_throws BoundsError jl(i5, 0.0) # not cached
        @test_throws BoundsError jl(i10, -1.0)
        @test_throws BoundsError jl(i10, jl.x[end] + 1.0)
        @test isapprox(jl(i10, 0.0), SymBoltz.sphericalbesselj(10, 0.0); atol = 1e-16)
        @test isapprox(jl(i10, jl.x[end]), SymBoltz.sphericalbesselj(10, jl.x[end]); atol = 1e-16)
        @test isapprox(jl(i10, 123.456), SymBoltz.sphericalbesselj(10, 123.456); atol)

        xs = range(jl.x[begin], jl.x[end], step=0.001)
        is = eachindex(ls)
        @test all(isapprox.(jl.(is', xs), SymBoltz.jl.(ls', xs); atol))
        @test (@ballocated $jl($i10, π)) == 0 # non-allocating

        j10(x) = jl(i10, x)
        @test isfinite(ForwardDiff.derivative(j10, π))
    end
end

@testset "Extend array" begin
    @test_throws "Cannot extend empty array" SymBoltz.extend_array([], 0)
    @test SymBoltz.extend_array(1.0:1.0:3.0, 0) == 1.0:1.0:3.0
    @test SymBoltz.extend_array(1.0:1.0:3.0, 4) == 1.0:0.2:3.0
end

@testset "Timeseries" begin
    sol = solve(prob)
    τs = SymBoltz.timeseries(sol; Nextra=1) # naive implementation could transform endpoints slightly through exp(log(τ))
    zs = sol(M.g.z, τs)
    as = sol(M.g.a, τs)
    @test τs[end-1] != τs[end] # ensure callback does not duplicate last point
    @test isapprox(as[end], 1.0; atol = 1e-12) && as[end] >= 1.0 # a(τ₀) ≈ 1.0, but not less, so root finding algorithms on time series work with different signs today
    @test isapprox(zs[end], 0.0; atol = 1e-12) && zs[end] <= 0.0 # z(τ₀) ≈ 0.0, but not more, so root finding algorithms on time series work with different signs today

    # Invert z to τ with root finding and check we get the same τ
    @test all(isapprox.(τs, SymBoltz.timeseries(sol, M.g.z, zs); atol = 1e-12))
    @test all(isapprox.(τs, SymBoltz.timeseries(sol, M.g.a, as); atol = 1e-12))

    # Invert z and ż to τ with Hermite spline and check we get the same τ
    @test all(isapprox.(τs, SymBoltz.timeseries(sol, M.g.z, M.g.ż, zs); atol = 1e-6)) # TODO: make more reliable
end

@testset "Source grid" begin
    τs = [1.0, 2.0]
    ks = [1.0, 10.0, 100.0]

    # scalar S: returns matrix of scalar sources
    Ss = source_grid(prob, M.τ + M.k, τs, ks)
    @test Ss isa Matrix{Float64}
    @test size(Ss) == (length(τs), length(ks))
    @test isequal(Ss, τs .+ transpose(ks))

    # vector S: returns matrix of vector sources
    Ss = source_grid(prob, [M.τ + M.k, M.τ * M.k], τs, ks)
    @test Ss isa Matrix{Vector{Float64}}
    @test size(Ss) == (length(τs), length(ks))
    @test isequal(getindex.(Ss, 1), τs .+ transpose(ks))
    @test isequal(getindex.(Ss, 2), τs .* transpose(ks))

    # same with interpolation with Cubic splines
    kinterp = CubicSplineInterpolator(range(1.0, 200.0, length=5))
    Ss = source_grid(prob, M.τ + M.k, τs, ks, kinterp)
    @test size(Ss) == (length(τs), length(ks))
    @test all(isapprox.(Ss, τs .+ transpose(ks)))

    # SVector S: returns matrix of SVectors
    Ss = source_grid(prob, SVector(M.τ + M.k, M.τ * M.k), τs, ks)
    @test Ss isa Matrix{SVector{2, Float64}}
    @test size(Ss) == (length(τs), length(ks))
    @test isequal(getindex.(Ss, 1), τs .+ transpose(ks))
    @test isequal(getindex.(Ss, 2), τs .* transpose(ks))

    # source_grid_chebyshev: scalar S
    kinterp = ChebyshevInterpolator(extrema(ks)..., 1)
    Ss = source_grid(prob, M.τ + M.k, τs, ks, kinterp)
    @test Ss isa Matrix{Float64}
    @test size(Ss) == (length(τs), length(ks))
    @test Ss ≈ τs .+ transpose(ks)

    # source_grid_chebyshev: SVector S
    Ss = source_grid(prob, SVector(M.τ + M.k, M.τ * M.k), τs, ks, kinterp)
    @test Ss isa Matrix{SVector{2, Float64}}
    @test size(Ss) == (length(τs), length(ks))
    @test getindex.(Ss, 1) ≈ τs .+ transpose(ks)
    @test getindex.(Ss, 2) ≈ τs .* transpose(ks)
end

@testset "Initial conditions" begin
    τini = prob.bg[1].tspan[1]
    ks = [1e2, 1e3]
    sol = solve(prob, ks)

    # Check that a ≈ √(Ωᵣ₀) * t
    Ωγ0 = M.γ.Ω₀
    Ων0 = M.ν.Ω₀
    Ωh0 = M.h.Ω₀ / M.h.Iρ₀ * 7π^4/120
    Ωr0 = Ωγ0 + Ων0 + Ωh0
    @test isapprox(sol[M.g.a][begin], sol[√(Ωr0)*M.τ][begin]; atol = 1e-10)

    # Check that τ ≈ 1 / g.ℋ
    @test isapprox(sol[M.τ][begin], sol[1/M.g.ℋ][begin]; atol = 1e-10)

    # Check that Fₗ(0) ∝ kˡ
    Fls = sol([M.γ.F0; collect(M.γ.F)], τini, ks)
    @test all(isapprox.(Fls[:,1] ./ Fls[:,2], map(l -> (ks[1]/ks[2])^l, 0:size(Fls)[1]-1))[1:4])

    # Check initial ratio of metric potentials
    @test all(isapprox.(sol(M.g.Φ / M.g.Ψ, τini, ks), sol((1+2/5*M.fν), τini); atol = 1e-4))

    # Check initial adiabatic perturbations
    species = [M.c, M.b, M.γ, M.ν, M.h]
    y0s = sol([s.δ/(1+s.w) for s in species], τini, ks) # should be equal for all species
    y1s = sol([s.u/M.k for s in species], τini, ks) # should be equal for all species
    @test isapprox(minimum(y0s), maximum(y0s); rtol = 1e-3)
    @test isapprox(minimum(y1s), maximum(y1s); rtol = 1e-3)
    y2s = sol([s.σ/M.k^2 for s in [M.ν, M.h]], τini, ks) # should be equal for massless and massive neutrinos
    @test isapprox(minimum(y2s), maximum(y2s); rtol = 1e-3)

    # Perturbations span the same τ as the background
    sol = solve(prob, [1e0, 1e1])
    @test all([ptsol.t[begin] == sol[M.τ][begin] && ptsol.t[end] == sol[M.τ][end] for ptsol in sol.pts])
end

@testset "Automatic background/thermodynamics splining" begin
    sol = solve(prob, 1.0) # solve with one perturbation mode to activate splining
    τs = SymBoltz.timeseries.(sol, log10(M.g.a), range(-8, 0, length=100)) # TODO a => as syntax
    function checkvar(var, atol, rtol)
        vals1 = sol(var, τs) # from background
        vals2 = sol(var, τs, 1.0) # from splined perturbations
        return all(isapprox.(vals1, vals2; atol, rtol))
    end
    @test checkvar(M.g.a, 1e-6, 0)
    @test checkvar(M.b.κ̇, 0, 1e-2)
    @test checkvar(M.b.κ, 0, 1e-4)
    @test checkvar(M.b.v, 1e-3, 0)
    @test checkvar(M.b.v̇, 0, 1e1) # TODO: improve
    @test checkvar(M.b.cₛ², 1e-4, 0)
    @test checkvar(M.b.T, 0, 1e-5)
    @test checkvar(M.b.Xe, 1e-5, 0)
end

@testset "Solve background+perturbations together (without splining background)" begin
    prob_nospline_dense = CosmologyProblem(M, pars; spline = false, ptsparse = false)
    prob_nospline_sparse = CosmologyProblem(M, pars; spline = false, ptsparse = true)
    @test_nowarn sprint(show, prob_nospline_dense)
    @test_nowarn sprint(show, prob_nospline_sparse)
    ks = [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
    @test issuccess(solve(prob_nospline_dense, ks))
    @test issuccess(solve(prob_nospline_sparse, ks))
    # TODO: compare number of timesteps with different tolerances? looks like runtime difference is proportional to difference in number of steps
    # TODO: although not splining is slower for normal values, can be opposite be true for AD?
end

# TODO: optionally spline observables, too, with cubic hermite splines using analytic derivatives, to save computations in perturbations (e.g. visibility function and derivatives for CMB)
@testset "Do not spline observed variables" begin
    @test_throws "not an unknown" SymBoltz.mtkcompile_spline(M, [M.g.H])
end

@testset "Primordial power spectrum pivot scale" begin
    h = pars[M.g.h]
    k = 0.05 * (L100/h) # 0.05/Mpc (≠ 0.05/(Mpc/h)) in units of H₀/c
    sol = solve(prob, k)
    @test sol[M.I.kpivot] ≈ k

    ks = 1.0:100.0
    sol = solve(prob, ks)
    P1 = sol(M.I.P, sol[M.τ][begin], ks)
    P2 = spectrum_primordial(ks, sol)
    @test all(isapprox.(P1, P2))
end

@testset "Time and optical depth today" begin
    ks = 1.0
    prob = CosmologyProblem(M, pars) # recreate since solution usually modifies problem parameters
    sol = solve(prob, ks)
    τ0 = sol[M.τ][end]
    @test sol(M.g.a, τ0) ≈ sol(M.g.a, τ0, ks) ≈ 1.0
    @test sol(M.χ, τ0) == sol(M.χ, τ0, ks) == 0.0
    @test sol(M.b.κ, τ0) == sol(M.b.κ, τ0, ks) == 0.0
end

@testset "Equal parameters in background and perturbation solutions" begin
    sol = solve(prob, [1.0, 10.0, 100.0])
    pars = [ # choose lots of background parameters that should be equal in perturbations
        M.g.h,
        M.c.Ω₀,
        M.b.Ω₀, M.b.YHe, M.b.fHe,
        M.γ.Ω₀, M.γ.T₀,
        M.ν.Ω₀, M.ν.T₀, M.ν.N,
        M.h.Ω₀, M.h.T₀, M.h.m, M.h.y₀, M.h.Iρ₀,
        M.Λ.Ω₀,
        M.I.As, M.I.kpivot, M.I.ns
    ]
    @test allequal([extrema(sol[M.τ]); map(pt -> extrema(pt.t), sol.pts)]) # background and perturbation should have equal timespans
    @test all(allequal([sol.bg[end].ps[par]; map(pt -> pt.ps[par], sol.pts)]) for par in pars)
end

@testset "Success checking" begin
    @test issuccess(solve(prob, 1.0))
    @test !issuccess(solve(prob, 0.0))
end

@testset "Consistent AD and FD derivatives of matter power spectrum" begin
    k = 10 .^ range(0, 3; length = 20)
    diffpars = [M.c.Ω₀, M.b.Ω₀] # TODO: h, ...
    probf = remake_function(prob, diffpars)
    function logP(logθ)
        θ = exp.(logθ)
        prob′ = probf(θ)
        P = spectrum_matter(prob′, k)
        return log.(P)
    end
    logθ = [log(pars[par]) for par in diffpars]
    ∂logP_∂logθ_ad = ForwardDiff.jacobian(logP, logθ)
    ∂logP_∂logθ_fd = FiniteDiff.finite_difference_jacobian(logP, logθ, Val{:central}; relstep = 1e-3)
    @test all(isapprox.(∂logP_∂logθ_ad, ∂logP_∂logθ_fd; atol = 1e-3))

    #= for debug plotting
    using CairoMakie
    fig = Figure()
    ax = Axis(fig[1, 1])
    for i in eachindex(diffpars)
        color = Makie.wong_colors()[i]
        alpha = 0.6
        lines!(ax, log10.(k), ∂logP_∂logθ_ad[:, i]; color, alpha, linestyle = :solid)
        lines!(ax, log10.(k), ∂logP_∂logθ_fd[:, i]; color, alpha, linestyle = :dash)
    end
    fig
    =#
end

@testset "Consistent AD and FD derivatives of CMB power spectrum" begin
    l = 25:25:1000
    jl = SphericalBesselCache(l)
    diffpars = [M.c.Ω₀, M.b.Ω₀] # TODO: h, ...
    probf = remake_function(prob, diffpars)
    function logDlTT(logθ)
        θ = exp.(logθ)
        prob′ = probf(θ)
        DlTT = spectrum_cmb(:TT, prob′, jl; normalization = :Dl)
        return log.(DlTT)
    end
    logθ = [log(pars[par]) for par in diffpars]
    ∂logDlTT_∂logθ_ad = ForwardDiff.jacobian(logDlTT, logθ)
    ∂logDlTT_∂logθ_fd = FiniteDiff.finite_difference_jacobian(logDlTT, logθ, Val{:central}; relstep = 1e-3) # 1e-4 screws up at small l
    @test all(isapprox.(∂logDlTT_∂logθ_ad, ∂logDlTT_∂logθ_fd; atol = 1e0)) # TODO: fix and decrease tolerance!!!

    #= for debug plotting
    using CairoMakie
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "l")
    for i in eachindex(diffpars)
        color = Makie.wong_colors()[i]
        alpha = 0.6
        label = replace("∂(log(Dₗ)) / ∂(log($(diffpars[i])))", "₊" => ".")
        lines!(ax, l, ∂logDlTT_∂logθ_ad[:, i]; color, alpha, linestyle = :solid, label = "$label (AD)")
        lines!(ax, l, ∂logDlTT_∂logθ_fd[:, i]; color, alpha, linestyle = :dash, label = "$label (FD)")
    end
    axislegend(ax; backgroundcolor = :transparent)
    fig
    =#
end

@testset "Parameter updater for dependent parameters" begin
    prob0 = CosmologyProblem(M, merge(pars, Dict(M.γ.T₀ => NaN)))
    Ω0total = M.γ.Ω₀ + M.ν.Ω₀ + M.c.Ω₀ + M.b.Ω₀ + M.h.Ω₀ + M.Λ.Ω₀
    getter = SymBoltz.getsym(prob0, [M.γ.T₀, M.γ.Ω₀, Ω0total]) # TODO: define Ω0total in model?
    @test all(isnan.(getter(prob0)))

    probf = remake_function(prob0, [M.γ.T₀])
    prob1 = probf([2.73])
    vals = getter(prob1)
    @test vals[1] == 2.73
    @test isfinite(vals[2])
    @test vals[3] ≈ 1.0
    @test issuccess(solve(prob1, 1.0))
end

@testset "Parameter updater and remake" begin
    stages(prob) = [prob.bg..., prob.pt]
    sameprob(prob1, prob2) = all(isequal(s1.u0, s2.u0) && isequal(s1.p.tunable, s2.p.tunable) for (s1, s2) in zip(stages(prob1), stages(prob2))) # u0 and p contain NaN
    Ω0total = M.γ.Ω₀ + M.ν.Ω₀ + M.h.Ω₀ + M.b.Ω₀ + M.c.Ω₀ + M.Λ.Ω₀
    ks = 10 .^ range(0, 3, length=10)
    prob_copy = deepcopy(prob)

    # update one parameter in all possible ways
    prob1 = CosmologyProblem(M, merge(pars, Dict(M.c.Ω₀ => 0.3)))
    newprobs1 = [
        remake_function(prob, M.c.Ω₀)(0.3),
        remake_function(prob, [M.c.Ω₀])([0.3]),
        remake_function(prob, (M.c.Ω₀,))((0.3,)),
        remake(prob, M.c.Ω₀ => 0.3),
        remake(prob, [M.c.Ω₀ => 0.3]),
        remake(prob, Dict(M.c.Ω₀ => 0.3)),
    ]

    # update several parameters in all possible ways
    prob2 = CosmologyProblem(M, merge(pars, Dict(M.c.Ω₀ => 0.3, M.g.h => 0.7, M.I.ns => 0.95)))
    newprobs2 = [
        remake_function(prob, [M.c.Ω₀, M.g.h, M.I.ns])([0.3, 0.7, 0.95]),
        remake_function(prob, [M.I.ns, M.g.h, M.c.Ω₀])([0.95, 0.7, 0.3]),
        remake_function(prob, (M.c.Ω₀, M.g.h, M.I.ns))((0.3, 0.7, 0.95)),
        remake(prob, [M.c.Ω₀ => 0.3, M.g.h => 0.7, M.I.ns => 0.95]),
        remake(prob, Dict(M.c.Ω₀ => 0.3, M.g.h => 0.7, M.I.ns => 0.95)),
        remake(remake(remake(prob, M.c.Ω₀ => 0.3), M.g.h => 0.7), M.I.ns => 0.95),
        remake_function(prob1, [M.g.h, M.I.ns])([0.7, 0.95]),
    ]

    for (freshprob, newprobs) in [(prob1, newprobs1), (prob2, newprobs2)]
        for newprob in newprobs
            @test sameprob(newprob, freshprob) # numerical values equal those of a problem created from scratch
            @test newprob.pt.ps[M.c.Ω₀] == 0.3
            @test all(stage.ps[M.c.Ω₀] == 0.3 for stage in newprob.bg if M.c.Ω₀ in SymBoltz.problem_symbols(stage))
            @test newprob.bg[end].ps[M.Λ.Ω₀] == newprob.pt.ps[M.Λ.Ω₀] == freshprob.pt.ps[M.Λ.Ω₀] != prob.pt.ps[M.Λ.Ω₀] # dependent parameter is updated
            @test newprob.bg[end].ps[Ω0total] ≈ newprob.pt.ps[Ω0total] ≈ 1.0
        end
    end
    @test sameprob(prob, prob_copy) # original problem is unchanged
    @test !sameprob(prob, prob1) && !sameprob(prob, prob2) && !sameprob(prob1, prob2)

    # only parameters given to the problem can be updated
    @test_throws "Cannot update" remake(prob, M.Λ.Ω₀ => 0.7) # dependent parameter
    @test_throws "Cannot update" remake_function(prob, [M.c.Ω₀, M.Λ.Ω₀])
    @test_throws "Cannot update" remake(prob, M.g.a => 1.0) # variable
    @test_throws "Cannot update" remake(prob, SymBoltz.k => 1.0) # wavenumber

    # differentiation
    function Pk(Ωc0)
        newprob = remake_function(prob, M.c.Ω₀)(Ωc0)
        return spectrum_matter(newprob, ks)
    end
    isnonzero(x) = isfinite(x) && !iszero(x)
    @test all(isnonzero.(Pk(0.3)))
    @test all(isnonzero.(Pk(0.3)))
    @test all(isnonzero.(ForwardDiff.derivative(Pk, 0.3)))
    @test all(isnonzero.(ForwardDiff.derivative(Pk, 0.3))) # twice (successive calls failed with earlier bug)
end

@testset "Dedicated background/perturbation solvers" begin
    bgsols = solvebg(prob) # TODO: @inferred
    @test bgsols isa Tuple{Vararg{SymBoltz.ODESolution}}

    ks = 1.0:1.0:10.0
    ptsol = solvept(prob.pt, bgsols, ks) # TODO: @inferred
    @test ptsol isa Vector{<:SymBoltz.ODESolution}

    # custom output_func for e.g. source function
    getS = SymBoltz.getsym(prob.pt, M.ST)
    τi, τ0 = extrema(bgsols[end].t)
    τs = range(τi, τ0, length = 768)
    ks = range(1.0, 1000.0, length = 1000)
    Ss = solvept(prob.pt, bgsols, ks; saveat = τs, output_func = (ptsol, _) -> getS(ptsol))
    Ss = stack(Ss)
    @test size(Ss) == (length(τs), length(ks))
end

@testset "Background differentiation test" begin
    diffpars = [M.g.h, M.c.Ω₀, M.b.Ω₀, M.γ.T₀, M.ν.N, M.h.m_eV, M.b.YHe, M.I.ln_As1e10, M.I.ns]
    probf = remake_function(prob, diffpars)
    τ0(θ) = solve(probf(θ))[M.τ][end]
    θ0 = [pars[par] for par in diffpars]
    dτ0_ad = ForwardDiff.gradient(τ0, θ0)
    dτ0_fd = FiniteDiff.finite_difference_gradient(τ0, θ0)
    @test all(isapprox.(dτ0_ad, dτ0_fd; atol = 1e-2))
    @test all(isapprox.(dτ0_ad[end-2:end], 0.0; atol = 1e-10))
    @test all(isapprox.(dτ0_fd[end-2:end], 0.0; atol = 1e-2))
end

@testset "Stability of different RECFAST models" begin
    M1 = ΛCDM(K = nothing; Hswitch = 0)
    M2 = ΛCDM(K = nothing; Heswitch = 0)
    M3 = ΛCDM(K = nothing; reionization = false)
    for M in [M1, M2, M3]
        prob = CosmologyProblem(M, pars)
        sol = solve(prob)
        @test issuccess(sol)
        @test all(sol[M.b.rec.XH⁺] .≤ 1 + 1e-5)
        @test all(sol[M.b.rec.XHe⁺] .≤ 1 + 1e-5)
        @test all(sol[M.b.rec.XHe⁺⁺/M.b.fHe] .≤ 1 + 1e-5)
    end
end

using QuasiMonteCarlo
function stability(M::System, ks, vary::Dict, nsamples; verbose = false, error = false, kwargs...)
    prob0 = CosmologyProblem(M, Dict(keys(vary) .=> NaN))
    pars = collect(keys(vary))
    probf = remake_function(prob0, pars)
    lo = [bound[1] for bound in values(vary)] # lower corner of parameter space
    hi = [bound[2] for bound in values(vary)] # uppper corner of parameter space
    samples = QuasiMonteCarlo.sample(nsamples, lo, hi, LatinHypercubeSample())
    nsuccess = 0
    if verbose
        println("Varying ", keys(vary))
        println("Solving for wavenumbers ", ks)
    end
    for sample in eachcol(samples)
        prob = probf(sample)
        sol = solve(prob, ks; kwargs...)
        if issuccess(sol)
            nsuccess += 1
        else
            solve(prob, ks; verbose=true, kwargs...) # solve again with verbose output for debugging
            error && Base.error("FAIL: ", sample)
        end
        verbose && println(issuccess(sol) ? "PASS" : "FAIL", ": ", sample)
    end
    return nsuccess / nsamples
end
vary = Dict(par => (0.5val, 1.5val) for (par, val) in pars) # ± 50% around fiducial values
ks = [1e0, 1e1, 1e2, 1e3]
@testset "Stability of problems throughout parameter space with Latin hypercube sampling" begin
    @test stability(M, ks, vary, 100; error = true) == 1.0 # 100%

    M1 = ΛCDM(K = nothing, Hswitch = 0; lmax)
    @test stability(M1, ks, vary, 100; error = true) == 1.0

    M2 = ΛCDM(K = nothing, Heswitch = 0; lmax)
    @test stability(M2, ks, vary, 100; error = true) == 1.0

    M3 = ΛCDM(K = nothing, reionization = false; lmax)
    @test stability(M3, ks, vary, 100; error = true) == 1.0
end

using SpecialFunctions: zeta as ζ
@testset "Momentum quadrature strategy" begin
    f(x) = 1 / (exp(x) + 1)
    for N in 1:5
        xs, Ws = SymBoltz.momentum_quadrature(f, 4)
        num(n) = sum(Ws .* xs .^ (n-2)) # numerical quadrature of ∫dx x^n/(exp(x)+1) from 0 to ∞
        anal(n) = factorial(n) * (1 - 1/2^n) * ζ(n+1) # <3 analytical expression for ∫dx x^n/(exp(x)+1) from 0 to ∞ (https://math.stackexchange.com/a/4111560)
        for n in 2:8
            @test isapprox(num(n), anal(n); rtol = 10.0^(-6+n-N))
        end
    end
end

@testset "CMB spectra" begin
    # Without l-interpolation
    ls = range(2, 2500; length = 200)
    jl = SphericalBesselCache(ls)
    DlsTT = spectrum_cmb(:TT, prob, jl; normalization = :Dl)
    DlsEE = spectrum_cmb(:EE, prob, jl; normalization = :Dl)
    @test all(isfinite.(DlsTT))
    @test all(isfinite.(DlsEE))

    # With l-interpolation using fallback cubic splines vs Chebyshev
    Dls = DlsTT # reference from above
    jl_cubic = SphericalBesselCache(range(2, 2500; length = 60))
    jl_cheb = SphericalBesselCache(ChebyshevInterpolator(2, 2500, 60))
    jl_chebint = SphericalBesselCache(ChebyshevIntegerInterpolator(2, 2500, 60))
    Dls_cubic = spectrum_cmb(:TT, prob, jl_cubic, ls; normalization = :Dl)
    Dls_cheb = spectrum_cmb(:TT, prob, jl_cheb, ls; normalization = :Dl)
    Dls_chebint = spectrum_cmb(:TT, prob, jl_chebint, ls; normalization = :Dl)
    @test isapprox(Dls_cubic, Dls; rtol = 1e-1)
    @test isapprox(Dls_cheb, Dls; rtol = 1e-4)
    @test isapprox(Dls_chebint, Dls; rtol = 1e-4)

    # Error with bad input
    @test_throws "outside the l-range" spectrum_cmb(:TT, prob, jl, 1:3000; normalization = :Dl)
end

@testset "Toggle threading" begin
    @test length(unique(fetch.(map(i -> SymBoltz.@spawnif(threadid(), true), 1:10)))) > 1
    @test only(unique(fetch.(map(i -> SymBoltz.@spawnif(threadid(), false), 1:10)))) == 1
end

@testset "Sparse Jacobian" begin
    # sparse background should work for ΛCDM, but since it is a small system the dense version should be a bit faster
    prob_sparse_bg = CosmologyProblem(M, pars; pt = false, bgjac = true, bgsparse = true)
    @test all(SymBoltz.issparse, prob_sparse_bg.bg)
    @test issuccess(solve(prob_sparse_bg))

    # with ΛCDM model
    k = [1e0, 1e1, 1e2, 1e3]
    sol = solve(prob_sparse, k; bgalg = SymBoltz.Rodas4P(linsolve = SymBoltz.LUFactorization()), ptalg = SymBoltz.KenCarp4(linsolve = SymBoltz.PureKLUFactorization()))
    @test issuccess(sol)

    M2 = RMΛ()
    pars2 = Dict(M2.m.Ω₀ => 0.3, M2.r.Ω₀ => 1e-5, M2.g.h => NaN, M2.r.T₀ => NaN)
    prob2 = CosmologyProblem(M2, pars2; bgsparse = true, ptsparse = true) # sparse background and perturbations
    @test all(SymBoltz.issparse, prob2.bg) && SymBoltz.issparse(prob2.pt)
    bgopts = (alg = SymBoltz.Rodas4P(linsolve = SymBoltz.PureKLUFactorization()),) # extra options still override the named ones
    ptopts = (alg = SymBoltz.KenCarp4(linsolve = SymBoltz.PureKLUFactorization()),)
    sol = solve(prob2, k; bgopts, ptopts)
    @test issuccess(sol)
end

@testset "Is-in-place and specialization level" begin
    ks = 10 .^ range(0, 3, length=5)
    for iip in (true, false), specialize in (SciMLBase.AutoSpecialize, SciMLBase.FullSpecialize)
        prob = CosmologyProblem(M, pars; iip, specialize)
        @test all(isinplace(stage) == iip for stage in prob.bg)
        @test isinplace(prob.pt) == iip
        @test all(isinplace(stage.f) == iip for stage in prob.bg)
        @test isinplace(prob.pt.f) == iip
        @test all(typeof(stage.f).parameters[2] == specialize for stage in prob.bg)
        @test typeof(prob.pt.f).parameters[2] == specialize
        @test issuccess(solve(prob, ks))

        # iip/specialize should propagate to problems created from an existing one
        probf = remake_function(prob, [M.c.Ω₀])
        prob2 = probf([0.3])
        @test all(isinplace(stage) == iip for stage in prob2.bg)
        @test isinplace(prob2.pt) == iip
        @test all(isinplace(stage.f) == iip for stage in prob.bg)
        @test isinplace(prob.pt.f) == iip
        @test all(typeof(stage.f).parameters[2] == specialize for stage in prob2.bg)
        @test typeof(prob2.pt.f).parameters[2] == specialize
        @test issuccess(solve(prob2, ks))
    end
end

@testset "Check compatibility between dense/sparse Jacobian and (non)linear solver" begin
    @test !issuccess(solve(prob_dense; bgopts = (alg = SymBoltz.Tsit5(), maxiters = 5))) # alg without linsolve
    @test !issuccess(solve(prob_sparse; bgopts = (alg = SymBoltz.Tsit5(), maxiters = 5)))
    @test issuccess(solve(prob_dense, 1.0)) # should automatically find compatible linsolves
    @test issuccess(solve(prob_sparse, 1.0))
    @test issuccess(solve(prob_dense, 1.0; bgopts = (alg = SymBoltz.Rodas5P(),), ptopts = (alg = SymBoltz.Rodas5P(),))) # should automatically find compatible linsolves
    @test issuccess(solve(prob_sparse, 1.0; bgopts = (alg = SymBoltz.Rodas5P(),), ptopts = (alg = SymBoltz.Rodas5P(),)))
    @test_throws "dense Jacobian must be solved with dense" solve(prob_dense; bgopts = (alg = SymBoltz.Rodas5P(linsolve = SymBoltz.PureKLUFactorization()),)) # has dense background
    @test_throws "sparse Jacobian must be solved with sparse" solve(prob_sparse, 1.0; ptopts = (alg = SymBoltz.Rodas5P(linsolve = SymBoltz.RFLUFactorization()),)) # has sparse perturbations
    @test issuccess(solve(prob_dense, 1.0; bgalg = SymBoltz.default_bgalg(prob_dense), ptalg = SymBoltz.default_ptalg(prob_dense; accuracy = 0)))
    @test issuccess(solve(prob_dense, 1.0; bgalg = SymBoltz.default_bgalg(prob_dense), ptalg = SymBoltz.default_ptalg(prob_dense; accuracy = 1)))
    @test issuccess(solve(prob_dense, 1.0; bgalg = SymBoltz.default_bgalg(prob_dense), ptalg = SymBoltz.default_ptalg(prob_dense; accuracy = 2)))
    @test issuccess(solve(prob_sparse, 1.0; bgalg = SymBoltz.default_bgalg(prob_sparse), ptalg = SymBoltz.default_ptalg(prob_sparse; accuracy = 0)))
    @test issuccess(solve(prob_sparse, 1.0; bgalg = SymBoltz.default_bgalg(prob_sparse), ptalg = SymBoltz.default_ptalg(prob_sparse; accuracy = 1)))
    @test issuccess(solve(prob_sparse, 1.0; bgalg = SymBoltz.default_bgalg(prob_sparse), ptalg = SymBoltz.default_ptalg(prob_sparse; accuracy = 2)))
end

@testset "Matter power spectrum with different arguments" begin
    modes = [:m, :c, :b, :cb, :cbh, :h]
    ks = [1e-1, 1e0, 1e1, 1e2]
    τs = [1.5, 3.0]
    sol = solve(prob, ks)
    @test size(spectrum_matter(modes, prob, ks, τs)) == (6, 2, 4) # general form
    @test size(spectrum_matter(modes, sol,  ks, τs)) == (6, 2, 4)
    @test size(spectrum_matter(modes, prob, ks)) == (6, 4) # omit τ; should use τ0
    @test size(spectrum_matter(modes, sol,  ks)) == (6, 4)
    @test size(spectrum_matter(prob, ks, τs)) == (2, 4) # omit modes; should use :m
    @test size(spectrum_matter(sol,  ks, τs)) == (2, 4)
    @test size(spectrum_matter(prob, ks)) == (4,) # omit modes and τ; should use :m and τ0
    @test size(spectrum_matter(sol,  ks)) == (4,)
end

@testset "Matter power spectrum converged to 0.1%" begin
    k = 10 .^ range(-1, 4, length=100)
    @time P0 = spectrum_matter(prob, k; bgalg = SymBoltz.default_bgalg(prob; stiff=true), bgabstol = 1e-10, bgreltol = 1e-10, ptalg = SymBoltz.default_ptalg(prob; accuracy=2), ptabstol = 1e-10, ptreltol = 1e-10)
    @time P  = spectrum_matter(prob, k)
    errs = abs.(P./P0 .- 1)
    @test all(errs .< 1e-3)
end

@testset "Zero allocations in ODE functions" begin
    for prob in [prob_dense, prob_sparse]
        sol = solve(prob, 1.0)
        for (subname, subsol) in [[(Symbol(:bg, i), bgsol) for (i, bgsol) in enumerate(sol.bg)]; (:pt, sol.pts[1])]
            subprob = subsol.prob
            u0 = subsol.u[begin]
            p = subprob.p
            t = subprob.tspan[begin]
            fout = similar(u0)
            Jout = isnothing(subprob.f.jac_prototype) ? zeros(length(u0), length(u0)) : subprob.f.jac_prototype
            fform = hasproperty(subprob.f, :f) ? "analytical" : "numerical"
            Jform = hasproperty(subprob.f, :jac) ? "analytical" : "numerical"
            Jform *= SymBoltz.issparse(subprob) ? "+sparse" : "+dense"
            println("Checking allocations for $subname with $fform f with output type $(typeof(fout)) and size $(size(fout))")
            @test (@ballocated $(subprob.f)($fout, $u0, $p, $t)) == 0
            println("Checking allocations for $subname with $Jform J with output type $(typeof(Jout)) and size $(size(Jout))")
            @test (@ballocated $(subprob.f.jac)($Jout, $u0, $p, $t)) == 0
        end
    end
end

@testset "Find and classify inner variables" begin
    # Parameters
    @test SymBoltz.isbackground(M.h.Ω₀)
    @test SymBoltz.isperturbation(M.h.Ω₀)
    @test string(only(SymBoltz.basevars(M.h.Ω₀))) == "h₊Ω₀"

    # Background variables (are also perturbation variables)
    @test SymBoltz.isbackground(M.h.ρ) && SymBoltz.isperturbation(M.h.ρ) && string(only(SymBoltz.basevars(M.h.ρ))) == "h₊ρ(τ)"
    @test SymBoltz.isbackground(D(M.h.ρ)) && SymBoltz.isperturbation(D(M.h.ρ)) && string(only(SymBoltz.basevars(D(M.h.ρ)))) == "h₊ρ(τ)" # differentiated
    @test SymBoltz.isbackground(M.h.E) && SymBoltz.isperturbation(M.h.E) && string(only(SymBoltz.basevars(M.h.E))) == "h₊E(τ)" # indexable
    @test SymBoltz.isbackground(M.h.E[1]) && SymBoltz.isperturbation(M.h.E[1]) && string(only(SymBoltz.basevars(M.h.E[1]))) == "h₊E(τ)" # indexed
    @test SymBoltz.isbackground(D(M.h.E)) && SymBoltz.isperturbation(D(M.h.E)) && string(only(SymBoltz.basevars(D(M.h.E)))) == "h₊E(τ)" # differentiated+indexable
    @test SymBoltz.isbackground(D(M.h.E[1])) && SymBoltz.isperturbation(D(M.h.E[1])) && string(only(SymBoltz.basevars(D(M.h.E[1])))) == "h₊E(τ)" # differentiated+indexed

    # Perturbation variables (are not background variables)
    @test !SymBoltz.isbackground(M.h.δ) && SymBoltz.isperturbation(M.h.δ) && string(only(SymBoltz.basevars(M.h.δ))) == "h₊δ(τ, k)"
    @test !SymBoltz.isbackground(D(M.h.δ)) && SymBoltz.isperturbation(D(M.h.δ)) && string(only(SymBoltz.basevars(D(M.h.δ)))) == "h₊δ(τ, k)" # differentiated
    @test !SymBoltz.isbackground(M.h.ψ) && SymBoltz.isperturbation(M.h.ψ) && string(only(SymBoltz.basevars(M.h.ψ))) == "h₊ψ(τ, k)" # indexable
    @test !SymBoltz.isbackground(M.h.ψ[1,1]) && SymBoltz.isperturbation(M.h.ψ[1,1]) && string(only(SymBoltz.basevars(M.h.ψ[1,1]))) == "h₊ψ(τ, k)" # indexed
    @test !SymBoltz.isbackground(D(M.h.ψ)) && SymBoltz.isperturbation(D(M.h.ψ)) && string(only(SymBoltz.basevars(D(M.h.ψ)))) == "h₊ψ(τ, k)" # differentiated+indexable
    @test !SymBoltz.isbackground(D(M.h.ψ[1,1])) && SymBoltz.isperturbation(D(M.h.ψ[1,1])) && string(only(SymBoltz.basevars(D(M.h.ψ[1,1])))) == "h₊ψ(τ, k)" # differentiated+indexed
end

@testset "Remove background initial conditions" begin
    @test isempty(SymBoltz.remove_background_initial_conditions!([D(M.g.a) ~ M.g.a/M.τ])) # should remove
    @test !isempty(SymBoltz.remove_background_initial_conditions!([M.g.Ψ ~ 20M.C / (15+4M.fν)])) # should keep
end

@testset "Split off a closed subsystem" begin
    @independent_variables t
    Dt = Differential(t)
    @variables x(t) y(t) z(t) w(t)
    @parameters c
    eqs = [Dt(x) ~ -x, Dt(y) ~ -y + x, z ~ c*y] # y depends on x, and the algebraic z depends on y
    sys = System(eqs, t, [x, y, z], [c]; initial_conditions = [x => 1.0, y => 1.0], name = :Chain)

    # x is closed on its own, and x and y are closed together
    @test isequal(unknowns(SymBoltz.split_system(sys, [x])), [x])
    @test isequal(unknowns(SymBoltz.split_system(sys, [x, y])), [x, y])

    # but y cannot be integrated without x, nor z without both
    @test_throws "depend on x(t)" SymBoltz.split_system(sys, [y])
    @test_throws "depend on x(t), y(t)" SymBoltz.split_system(sys, [z])
    @test_throws "w(t) is not an unknown or parameter" SymBoltz.split_system(sys, [w])

    # the system must be flattened first
    syssub = System(eqs, t, [x, y, z], [c]; systems = [System([Dt(z) ~ -z], t, [z], []; name = :sub)], name = :Chain)
    @test_throws "flatten it first" SymBoltz.split_system(syssub, [x])
end

@testset "Shooting method" begin
    M = BDΛCDM()

    # 1) unspecified ΩΛ0, constrained ℋ = 1 today
    pars1 = merge(parameters_Planck18(M), Dict(M.G.ω => 100.0, M.G.ϕini => 0.95, M.G.ϕ̇ini => 0.0))
    prob1 = CosmologyProblem(M, pars1, Dict(M.Λ.Ω₀ => 0.5), [M.g.ℋ ~ 1])
    sol1 = solve(prob1)
    @test issuccess(sol1)
    @test sol1[M.g.ℋ][end] ≈ 1.0 atol=1e-4
    @test sol1.bg[1][D(M.G.ϕ)][begin] == 0.0 # TODO: sol1[D(M.G.ϕ)] when later stages recognize splined derivatives

    # 1) same, but with bracketing root-finder
    prob1_bracket = CosmologyProblem(M, pars1, Dict(M.Λ.Ω₀ => (0.5, 1.0)), [M.g.ℋ ~ 1])
    sol1_bracket = solve(prob1_bracket)
    @test issuccess(sol1_bracket) && sol1_bracket[M.g.ℋ][end] ≈ 1.0 && sol1_bracket.bg[1][D(M.G.ϕ)][begin] == 0.0
    @test sol1_bracket[M.Λ.Ω₀] ≈ sol1[M.Λ.Ω₀] atol=1e-4

    # 2) unspecified ΩΛ0 and ϕini
    pars2 = merge(parameters_Planck18(M), Dict(M.G.ω => 100.0, M.G.ϕ̇ini => 0.0))
    prob2 = CosmologyProblem(M, pars2, Dict(M.G.ϕini => 1-1/(1+M.G.ω/5), M.Λ.Ω₀ => 0.5), [M.g.ℋ ~ 1, M.G.G ~ 1]) # ω-dependent ϕini ≈ 0.95
    sol2 = solve(prob2)
    @test issuccess(sol2)
    @test isapprox(sol2[M.g.ℋ][end], 1.0; atol = 1e-5)
    @test isapprox(sol2[M.G.G][end], 1.0; atol = 1e-5)
    @test sol2.bg[1][D(M.G.ϕ)][begin] == 0.0

    # start shooting with valid but bad initial guess
    prob_bad = CosmologyProblem(M, pars1, Dict(M.Λ.Ω₀ => 5.0), [M.g.ℋ ~ 1])
    sol_bad = solve(prob_bad)
    @test issuccess(sol_bad)
    @test sol_bad[M.g.ℋ][end] ≈ 1.0 atol=1e-4

    # initial shooting guess in bad/unstable region
    prob_stupid = CosmologyProblem(M, pars1, Dict(M.Λ.Ω₀ => -1.0), [M.g.ℋ ~ 1])
    @test_throws "Shooting failed to converge" solve(prob_stupid)

    # bracketing method with both initial guesses of the same sign
    prob_stupid = CosmologyProblem(M, pars1, Dict(M.Λ.Ω₀ => (0.0, 0.5)), [M.g.ℋ ~ 1])
    @test_throws "Shooting failed to converge" solve(prob_stupid)

    # illegal input
    @test_throws "Only parameters can be specified" CosmologyProblem(M, merge(pars1, Dict(D(M.G.ϕ) => 0.0)), Dict(M.Λ.Ω₀ => 0.5), [M.g.ℋ ~ 1])
    @test_throws "Cannot update" remake(prob1, Dict(M.G.ϕ => 0.9))
    @test_throws "Got 2 shooting parameters" CosmologyProblem(M, pars2, Dict(M.G.ϕ => 0.95, M.Λ.Ω₀ => 0.5), [M.g.ℋ ~ 1])
    @test_throws "Shooting with multiple parameters requires scalar guesses, but got interval guesses" CosmologyProblem(M, pars2, Dict(M.G.ϕ => (0.5, 1.5), M.Λ.Ω₀ => (0.5, 1.0)), [M.g.ℋ ~ 1, M.G.G ~ 1])
    @test_throws "requires nonbracketing" solve(prob1; shootalg = SymBoltz.default_shootalg(prob1_bracket))
    @test_throws "requires nonbracketing" solve(prob2; shootalg = SymBoltz.default_shootalg(prob1_bracket))
    @test_throws "requires bracketing" solve(prob1_bracket; shootalg = SymBoltz.default_shootalg(prob1))

    # test that Base.show works for different shooting guess/condition combinations
    @test_nowarn sprint(show, prob1)
    @test_nowarn sprint(show, sol1)
    @test_nowarn sprint(show, prob1_bracket)
    @test_nowarn sprint(show, sol1_bracket)
    @test_nowarn sprint(show, prob2)
    @test_nowarn sprint(show, sol2)
    @test_nowarn sprint(show, prob_bad)
    @test_nowarn sprint(show, sol_bad)

    # shooting in model
    vars = @variables a(τ) ℋ(τ)
    pars = @parameters Ωr0 Ωm0 ΩΛ0 [shoot=true]
    eqs = [ℋ ~ √(Ωr0/a^4 + Ωm0/a^3 + ΩΛ0) * a, D(a) ~ a*ℋ]
    initialization_eqs = [a^2*ℋ^2 ~ a^2/τ^2] # this form avoids 1/a and is more stable
    guesses = Dict(ΩΛ0 => 1 - Ωr0 - Ωm0, a => √(Ωr0) * τ)
    constraints = [ℋ ~ 1]
    @named M = System(eqs, τ, vars, pars; initialization_eqs, guesses, constraints)
    @test Set(keys(SymBoltz.shootvars(M))) == Set(ΩΛ0)
    pars = Dict(Ωr0 => 1e-5, Ωm0 => 0.3)
    prob = CosmologyProblem(M, pars)
    @test isnothing(prob.pt)
    sol = solve(prob)
    @test sol[ℋ][end] ≈ 1 && sol[Ωr0 + Ωm0 + ΩΛ0] ≈ 1

    # a guess passed to CosmologyProblem must override the one declared in the model
    @test CosmologyProblem(M, pars, Dict(ΩΛ0 => 0.123)).shoot[ΩΛ0] == 0.123
    # ... and must actually be used: a guess in a bad region makes the shooting fail,
    # while the model's own guess would have converged
    @test_throws "Shooting failed to converge" solve(CosmologyProblem(M, pars, Dict(ΩΛ0 => -1e3)))
    # an interval guess selects a bracketing solver, so it must reach the solver too
    @test issuccess(solve(CosmologyProblem(M, pars, Dict(ΩΛ0 => (0.5, 1.0)))))

    # shooting with numerical continuity equations
    vars = @variables a(τ) ℋ(τ) ρ(τ) ρr(τ) ρm(τ) ρΛ(τ)
    pars = @parameters ρri ρmi ρΛi
    eqs = [
        ℋ ~ √(8π/3*ρ) * a
        ρ ~ ρr + ρm + ρΛ
        D(a) ~ ℋ*a
        D(ρr) ~ -4ℋ*ρr
        D(ρm) ~ -3ℋ*ρm
        ρΛ ~ ρΛi
    ]
    initial_conditions = [
        ρr => ρri
        ρm => ρmi
    ]
    initialization_eqs = [
        ℋ ~ 1 / τ
    ]
    guesses = [
        a => √(Ωr0) * τ
    ]
    @named M = System(eqs, τ, vars, pars; initial_conditions, initialization_eqs, guesses)
    p = Dict(
        ρri => NaN,
        ρmi => NaN,
        ρΛi => NaN,
    )
    prob = CosmologyProblem(M, p, Dict(ρri => 1e-5/1e-8^4, ρmi => 0.3/1e-8^3, ρΛi => 0.7), [8π/3*ρr ~ 1e-5, 8π/3*ρm ~ 0.3, 8π/3*ρΛ ~ (1-0.3-1e-5)])
    sol = solve(prob; verbose = true)
    @test isapprox(sol[a/τ][begin], sol[√(8π/3*ρr*a^4)][begin]; atol = 1e-6)

    # TODO: require all shooting variables to be @parameters
    M = ΛCDM(; Λ = SymBoltz.cosmological_constant(SymBoltz.metric(); analytical = false))
    p = parameters_Planck18(M)
    shoot_guesses = Dict(M.Λ.ρ => 0.1)
    shoot_conditions = [M.g.ℋ ~ 1]
    @test_throws "must be declared with @parameters" CosmologyProblem(M, p, shoot_guesses, shoot_conditions)
end

@testset "Underdetermined/overdetermined initialization" begin
    # Underdetermined
    M2 = flatten(M)
    ieqs = ModelingToolkit.get_initialization_eqs(M2)
    deleteat!(ieqs, findfirst(eq -> isequal(eq.lhs, D(M2.a)), ieqs)) # delete ℋ = 1/τ
    @test_throws ModelingToolkit.StateSelection.ExtraVariablesSystemException CosmologyProblem(M2, pars)

    # Overdetermined
    M2 = flatten(M)
    ieqs = ModelingToolkit.get_initialization_eqs(M2)
    eq = ieqs[findfirst(eq -> isequal(eq.lhs, M2.Ψ), ieqs)] # IC for Ψ
    push!(ieqs, eq.lhs ~ 2eq.rhs) # overconstrain
    @test CosmologyProblem(M2, pars; pt = false) isa CosmologyProblem
    @test_throws ModelingToolkit.StateSelection.ExtraEquationsSystemException CosmologyProblem(M2, pars)
end

@testset "Expand" begin
    @variables w(τ) ρ(τ) ℋ(τ)
    eqs = [
        D(ρ) ~ -3ℋ*(1+w)ρ
        w ~ 1//3
        ℋ  ~ 123
    ]
    @test isequal(expandeq(eqs, D(ρ); protect = Set(ℋ)), -4ℋ*ρ)
end

@testset "High lmax" begin
    M = ΛCDM(lmax = 32)
    prob = CosmologyProblem(M, pars)
    @test issuccess(solve(prob, [1e-1, 1e0, 1e1, 1e2, 1e3]))
end

@testset "CLASS comparison" begin
    @test endswith(pwd(), "test") # should be in test/ subdirectory
    if false # switch to true to generate CLASS output
        using CLASS
        prob_class = CLASSProblem(
            "output" => "mPk, tCl, pCl, lCl",
            "ic" => "ad",
            "modes" => "s",
            "gauge" => "newtonian",
            "h" => pars[M.g.h],
            "T_cmb" => pars[M.γ.T₀],
            "l_max_g" => lmax,
            "l_max_pol_g" => lmax,
            "Omega_b" => pars[M.b.Ω₀],
            "YHe" => pars[M.b.YHe],
            "recombination" => "recfast",
            "recfast_Hswitch" => 1,
            "recfast_Heswitch" => 6,
            "reio_parametrization" => "reio_camb",
            "Omega_cdm" => pars[M.c.Ω₀],
            "N_ur" => pars[M.ν.N],
            "N_ncdm" => 1,
            "deg_ncdm" => pars[M.h.N],
            "m_ncdm" => pars[M.h.m_eV],
            "T_ncdm" => (4/11)^(1/3),
            "l_max_ur" => lmax,
            "l_max_ncdm" => lmax,
            "ln_A_s_1e10" => pars[M.I.ln_As1e10],
            "n_s" => pars[M.I.ns],
            "Omega_Lambda" => 0.0, # determine automatically
            "w0_fld" => -1.0,
            "wa_fld" => 0.0,
            "Omega_k" => 0.0,
            "Omega_scf" => 0.0,
            "Omega_dcdmdr" => 0.0,
            "tight_coupling_approximation" => 5, # compromise_CLASS; cannot turn off, and more accurate second_order_CLASS is incompatible with newtonian gauge
            "tight_coupling_trigger_tau_c_over_tau_h" => 1e-2, # cannot turn off
            "tight_coupling_trigger_tau_c_over_tau_k" => 1e-3, # cannot turn off
            "radiation_streaming_approximation" => 3, # turn off
            "ur_fluid_approximation" => 3, # turn off
            "ncdm_fluid_approximation" => 3, # turn off
        )
        sol_class = solve(prob_class)
        ks_class = sol_class[:pk][!, "k (h/Mpc)"] * L100
        Pks_class = sol_class[:pk][!, "P (Mpc/h)^3"] / L100^3
        ls_class = sol_class[:cl][!, "l"]
        DlTTs_class = sol_class[:cl][!, "TT"]
        DlEEs_class = sol_class[:cl][!, "EE"]
        Dlϕϕs_class = sol_class[:cl][!, "phiphi"]
        open("./class_Pk.dat", "w") do f # tests run from test/ directory
            writedlm(f, [ks_class Pks_class])
        end
        open("./class_Cl.dat", "w") do f
            writedlm(f, [ls_class DlTTs_class DlEEs_class Dlϕϕs_class])
        end
    end

    # Matter power spectrum
    Pk_class = readdlm("./class_Pk.dat")
    ks_class, Pks_class = Pk_class[:, 1], Pk_class[:, 2]
    ks = ks_class # solve at same wavenumbers as CLASS
    Pks = spectrum_matter(prob, ks)
    @test isapprox(Pks, Pks_class; rtol = 1e-3)

    # CMB power spectrum
    Cl_class = readdlm("./class_Cl.dat")
    ls_class, DlTTs_class, DlEEs_class, Dlϕϕs_class = Cl_class[:, 1], Cl_class[:, 2], Cl_class[:, 3], Cl_class[:, 4]
    ls = unique(Int.(round.(exp.(range(log(ls_class[begin]), log(ls_class[end]), length=200)))))
    jl = SphericalBesselCache(ls)
    Dls = spectrum_cmb([:TT, :EE, :ψψ], prob, jl, ls_class; normalization = :Dl)
    @test isapprox(Dls[:, 1], DlTTs_class; rtol = 2e-3)
    @test isapprox(Dls[:, 2], DlEEs_class; rtol = 2e-3)
    @test isapprox(Dls[ls_class .< 11, 3], Dlϕϕs_class[ls_class .< 11]; rtol = 1e-2) # full line-of-sight integration below l_limber (l = 2 has higher error and breaks isapprox(...; rtol = 2e-3) for all l)
    @test isapprox(Dls[ls_class .≥ 11, 3], Dlϕϕs_class[ls_class .≥ 11]; rtol = 1e-3) # Limber approximation enabled
end

@testset "Error if nonfinite error message" begin
    # Vector input of scalars
    x = rand(4)
    @test_nowarn SymBoltz.error_if_nonfinite(x)
    x[3] = NaN
    @test_throws "NaN at CartesianIndex(3" SymBoltz.error_if_nonfinite(x) # julia<1.13 says CartesianIndex(3,); julia≥1.13 says CartesianIndex(3)

    # Matrix input of SVector
    x = rand(SVector{2, Float64}, 4, 5)
    @test_nowarn SymBoltz.error_if_nonfinite(x)
    x[3] = x[3] .+ NaN
    @test_throws "[NaN, NaN] at CartesianIndex(3, 1)" SymBoltz.error_if_nonfinite(x)
end

@testset "Spherical Bessel function and CMB power spectrum with non-integer ℓ" begin
    ls = loggrid(2, 2500; length = 100)
    jl = SphericalBesselCache(ls)
    xs = range(jl.x[begin], jl.x[end]; length = 1000)
    @test all(isfinite, jl.(transpose(eachindex(jl.l)), xs))
    #plot(); for i in eachindex(jl.l) plot!(x -> jl(i, x), xlims = (0, 10), label = "l = $(jl.l[i])") end; plot!()

    ls_all = 2:2500
    Dls = spectrum_cmb(:TT, prob, jl, ls_all; normalization = :Dl)
    @test all(isfinite, Dls)
    #plot(ls_all, Dls; xscale = :log10)
end

@testset "Interpolation" begin
    x = range(0.0, 10.0; length=20)
    interp = CubicSplineInterpolator(x)
    @test eltype(interp) == eltype(x)
    @test issorted(interp)
    x′ = range(x[begin], x[end]; length = 1000)
    y′ = interpolate(interp, sin.(x), x′)
    @test all(interpolate(x, sin.(x), x′) .== y′) # should fall exactly back to cubic spline interpolation
    @test isapprox(y′, sin.(x′); atol = 1e-1)

    # same with integer-only x
    x = 0:10
    interp = CubicSplineInterpolator(x)
    @test eltype(interp) == eltype(x)
    @test issorted(interp)
    y′ = interpolate(interp, sin.(x), x′)
    @test all(interpolate(x, sin.(x), x′) .== y′) # should fall exactly back to cubic spline interpolation
    @test isapprox(y′, sin.(x′); atol = 1e-0)

    x = range(0.0, 10.0; length=20)
    interp = ChebyshevInterpolator(x[begin], x[end], 20)
    @test eltype(interp) == eltype(x)
    @test issorted(interp)
    y′ = interpolate(interp, sin.(interp), x′)
    @test isapprox(y′, sin.(x′); atol = 1e-10) # more accurate than cubic splines
    @test isapprox(interp.ws, SymBoltz.baryweights(interp.xs); atol = 1e-12)

    xbreak = (0.0, 5.0, 10.0)
    interp = PiecewiseChebyshevInterpolator(xbreak, (10, 20))
    @test eltype(interp) == eltype(xbreak)
    @test issorted(interp)
    y′ = interpolate(interp, sin.(interp), x′)
    @test isapprox(y′[x′ .≤ 5.0], sin.(x′[x′ .≤ 5.0]); atol = 1e-4) # lower order, less accurate
    @test isapprox(y′[x′ .≥ 5.0], sin.(x′[x′ .≥ 5.0]); atol = 1e-12) # higher order, more accurate

    interp = ChebyshevIntegerInterpolator(0, 100, 22)
    @test eltype(interp) <: Integer
    @test issorted(interp.xs)
    @test all(isinteger, interp.xs)
    @test allunique(interp.xs)
    @test extrema(interp) == (0, 100)
    x′ = range(interp[begin], interp[end]; length = 1000)
    y′ = interpolate(interp, sin.(π/30 .* interp), x′)
    @test isapprox(y′, sin.(π/30 .* x′); atol = 1e-10)
    @test_throws "collide" ChebyshevIntegerInterpolator(0, 100, 23)
end

@testset "Model with logarithmic scale factor as independent variable" begin
    @independent_variables b
    D = Differential(b)
    pars = @parameters Ωr0 Ωm0 ΩΛ0
    vars = @variables a(b) ρ(b) ρr(b) ρm(b) H(b) ℋ(b) τ(b) Φ(b, k) Ψ(b, k) δρ(b, k) δr(b, k) θr(b, k) δm(b, k) θm(b, k)
    eqs = [
        # background (Friedmann equation)
        a ~ exp(b)
        ρr ~ 3/8π * Ωr0/a^4
        ρm ~ 3/8π * Ωm0/a^3
        ρ ~ ρr + ρm + 3/8π * ΩΛ0
        H ~ √(8π/3 * ρ)
        ℋ ~ a * H
        D(τ) ~ 1 / ℋ
        # perturbations (Newtonian gauge, no anisotropic stress)
        Ψ ~ Φ
        δρ ~ ρr*δr + ρm*δm
        D(Φ) ~ (-4π/3*a^2/ℋ*δρ - k^2/(3ℋ)*Φ - ℋ*Ψ) / ℋ
        D(δr) ~ -4/3*θr/ℋ + 4*D(Φ)
        D(θr) ~ (k^2*δr/4 + k^2*Ψ) / ℋ
        D(δm) ~ -θm/ℋ + 3*D(Φ)
        D(θm) ~ -θm + k^2*Ψ/ℋ
    ]
    initial_conditions = [
        τ => 2 * (√(Ωr0 + Ωm0*a) - √(Ωr0)) / Ωm0
        Φ => 20/15
        δr => -2*Φ
        δm => -3/2*Φ
        θr => 1/2*k^2*τ*Φ
        θm => 1/2*k^2*τ*Φ
        ΩΛ0 => 1 - Ωm0 - Ωr0
    ]
    M = complete(System(eqs, b, vars, [pars; k]; initial_conditions, name = :RMΛ))
    p = Dict(M.Ωr0 => 1e-4, M.Ωm0 => 0.3)

    # use timespan that goes past today, but terminate when a = 1 with a callback
    prob = CosmologyProblem(M, p; tspan = (-9.0, 0.1), terminate = M.a ~ 1)
    ks = 10.0 .^ 0:0.5:3
    sol = solve(prob, ks)
    @test issuccess(sol)
    @test sol[M.a][end] ≈ 1.0
    @test sol[M.H][end] ≈ 1.0

    # integrate exactly to a = 1 with no termination condition
    prob = CosmologyProblem(M, p; tspan = (-9.0, 0.0), terminate = nothing)
    sol = solve(prob, ks)
    @test issuccess(sol)
    @test sol[M.a][end] ≈ 1.0
    @test sol[M.H][end] ≈ 1.0
end

@testset "Interacting background integrated backwards with b = log(a)" begin
    # 1) Simplified interacting models
    function interacting_model(f; name = :QΛCDM)
        @independent_variables b # = log(a)
        D = Differential(b)
        pars = @parameters Ωr0 Ωb0 Ωc0 ΩΛ0 w0 wa α
        vars = @variables a(b) ρ(b) ρr(b) ρb(b) ρc(b) ρΛ(b) H(b) ℋ(b) Q(b)
        eqs = [
            a ~ exp(b)
            ρr ~ 3/8π * Ωr0/a^4
            ρb ~ 3/8π * Ωb0/a^3
            ρ ~ ρr + ρb + ρc + ρΛ
            H ~ √(8π/3 * ρ)
            ℋ ~ a * H
            Q ~ α * f(ρc, ρΛ) # = a*Q/ℋ in the equations above
            D(ρc) ~ -3ρc + Q
            D(ρΛ) ~ -3*(1+w0+wa*(1-a))*ρΛ - Q
        ]
        initial_conditions = [
            ρc => 3/8π * Ωc0 # today # TODO: override when solving forward?
            ρΛ => 3/8π * ΩΛ0 # today # TODO: override when solving forward?
            ΩΛ0 => 1 - Ωr0 - Ωb0 - Ωc0
        ]
        return complete(System(eqs, b, vars, pars; initial_conditions, name))
    end

    M1 = interacting_model((ρc, ρΛ) -> ρc)
    M2 = interacting_model((ρc, ρΛ) -> ρΛ)
    M3 = interacting_model((ρc, ρΛ) -> ρc + ρΛ)
    M4 = interacting_model((ρc, ρΛ) -> ρc*ρΛ / (ρc+ρΛ))
    p = Dict(
        M1.Ωr0 => 1e-5,
        M1.Ωb0 => 0.05,
        M1.Ωc0 => 0.30,
        M1.w0 => -1.1,
        M1.wa => 0.2,
        M1.α => 50.0,
    )
    prob1 = CosmologyProblem(M1, p; tspan = (0, -8), terminate = nothing)
    prob2 = CosmologyProblem(M2, p; tspan = (0, -8), terminate = nothing)
    prob3 = CosmologyProblem(M3, p; tspan = (0, -8), terminate = nothing)
    prob4 = CosmologyProblem(M4, p; tspan = (0, -8), terminate = nothing)
    @test issuccess(solve(prob1))
    @test issuccess(solve(prob2))
    @test issuccess(solve(prob3))
    @test issuccess(solve(prob4))

    # 2) Full interacting models (in own module to not leak globals)
    Base.include(Module(), "interacting.jl") # do not pollute global namespace
end

@testset "Automatic background stages" begin
    @independent_variables t
    D = Differential(t)
    @variables x(t) y(t) [backwards = true] z(t) w(t)

    # starting backwards gives 2 stages ([y], [x, z, w]), while starting forwards would give 3 ([x, w], [y], [z])
    sys = System([D(x) ~ -x, D(y) ~ -y, D(z) ~ y - z, D(w) ~ x], t; name = :sys)
    @test isequal(SymBoltz.split_stages(sys), (([y], [x, z, w]), (true, false)))

    # mutually dependent variables are integrated together
    sys = System([D(x) ~ z, D(y) ~ -y, D(z) ~ x], t; name = :sys)
    @test isequal(SymBoltz.split_stages(sys), (([x, z], [y]), (false, true)))
    sys = System([D(x) ~ y, D(y) ~ x], t; name = :sys)
    @test_throws "must be integrated in the same direction" SymBoltz.split_stages(sys)
end

@testset "Solving model with mixed forward-backward direction in background" begin
    @independent_variables b
    D = Differential(b)
    pars = @parameters Ωr0 Ωm0 ΩΛ0 k
    vars = @variables a(b) ρ(b) ρr(b) [backwards = true] ρm(b) [backwards = true] ρΛ(b) H(b) ℋ(b) τ(b) Φ(b,k) δρ(b,k) δr(b,k) θr(b,k) δm(b,k) θm(b,k)
    eqs = [
        # background equations
        a ~ exp(b)
        ρ ~ ρr + ρm + ρΛ
        H ~ √(8π/3 * ρ)
        ℋ ~ a * H
        D(τ) ~ 1 / ℋ
        D(ρr) ~ -4*ρr
        D(ρm) ~ -3*ρm
        ρΛ ~ 3/8π * ΩΛ0
        # perturbation equations
        δρ ~ ρr*δr + ρm*δm
        D(Φ) ~ (-4π/3*a^2/ℋ*δρ - k^2/(3ℋ)*Φ - ℋ*Φ) / ℋ
        D(δr) ~ -4/3*θr/ℋ + 4*D(Φ)
        D(θr) ~ k^2 * (δr/4 + Φ) / ℋ
        D(δm) ~ -θm/ℋ + 3*D(Φ)
        D(θm) ~ -θm + k^2*Φ/ℋ
    ]
    initial_conditions = [
        ΩΛ0 => 1 - Ωr0 - Ωm0
        ρr => 3/8π * Ωr0 / a^4
        ρm => 3/8π * Ωm0 / a^3
        τ => a / √(Ωr0)
        Φ => 20/15
        δr => -2*Φ
        δm => -3/2*Φ
        θr => 1/2*k^2*τ*Φ
        θm => 1/2*k^2*τ*Φ
    ]
    M = complete(System(eqs, b, vars, pars; initial_conditions, name = :RMΛ))
    p = Dict(M.Ωr0 => 1e-4, M.Ωm0 => 0.3)
    prob = CosmologyProblem(M, p; bg = ([ρr, ρm], [τ]), tspan = (-8.0, 0.0), terminate = nothing)
    ks = 10.0 .^ (0:3)
    sol = solve(prob, ks)
    @test issuccess(sol)

    # the same stages are detected automatically from the dependencies and backwards metadata
    probauto = CosmologyProblem(M, p; tspan = (-8.0, 0.0), terminate = nothing)
    @test probauto.bg[1].tspan == reverse(probauto.tspan) # the first stage is backwards
    @test probauto.bg[2].tspan == probauto.tspan # the second stage is forwards
    @test issetequal(unknowns(probauto.bg[1].f.sys), [ρr, ρm])
    @test solve(probauto, ks)(M.δm, 0.0, ks) ≈ sol(M.δm, 0.0, ks)

    # τ cannot be split off before the densities it depends on
    @test_throws "depend on ρm(b), ρr(b)" CosmologyProblem(M, p; bg = ([τ], [ρr, ρm]), tspan = (-8.0, 0.0), terminate = nothing)

    # the variables of one stage must be declared with the same direction
    @test_throws "must be integrated in the same direction" CosmologyProblem(M, p; bg = ([ρr, ρm, τ],), tspan = (-8.0, 0.0), terminate = nothing)

    # the backward solve hits its boundary condition ρ(a=1) = 3/8π*Ω₀ exactly
    @test sol(M.ρr, 0.0) == 3/8π * p[M.Ωr0]
    @test sol(M.ρm, 0.0) == 3/8π * p[M.Ωm0]

    # variables of earlier stages are evaluated from their spline in the last stage
    bs = range(-8.0, 0.0, length = 5)
    @test sol(M.ρr, bs) ≈ sol.bg[1](bs; idxs = M.ρr).u rtol = 1e-4
    @test sol([M.ρr, M.τ], bs) == [sol(M.ρr, bs)'; sol(M.τ, bs)']

    # the stages can be solved one by one
    bgsol1 = solvebg(prob.bg[1])
    bgsol2 = solvebg(SymBoltz.setupbg(prob.bg[2], (bgsol1,)))
    @test solvept(prob.pt, (bgsol1, bgsol2), ks)[end].u[end] ≈ sol.pts[end].u[end]

    # ForwardDiff should differentiate through the whole backward-forward-perturbation chain
    diffpars = [M.Ωr0, M.Ωm0]
    probf = remake_function(prob, diffpars)
    tol = 1e-10
    function output(θ)
        solθ = solve(probf(θ), ks; bgopts = (reltol = tol, abstol = tol), ptopts = (reltol = tol, abstol = tol)) # not sol, which would overwrite the outer one
        return [solθ(M.τ, 0.0); vec(solθ(M.Φ, 0.0, ks)); vec(solθ(M.δm, 0.0, ks))]
    end
    θ0 = [p[par] for par in diffpars]
    J_ad = ForwardDiff.jacobian(output, θ0)
    J_fd = FiniteDiff.finite_difference_jacobian(output, θ0, Val{:central}; relstep = 1e-2, absstep = 1e-6) # absstep keeps Ωr0 = 1e-4 positive; smaller steps are dominated by ODE solver noise
    @test all(isapprox.(J_ad, J_fd; rtol = 1e-3))
end

@testset "Solving model with background integrated only backwards" begin
    @independent_variables b
    D = Differential(b)
    pars = @parameters Ωr0 Ωm0 ΩΛ0 h As ns k
    vars = @variables a(b) ρ(b) ρr(b) [backwards = true] ρm(b) [backwards = true] ρΛ(b) H(b) ℋ(b) Φ(b,k) δρ(b,k) δr(b,k) θr(b,k) δm(b,k) θm(b,k) Δm(b,k)
    eqs = [
        # background equations
        a ~ exp(b)
        ρ ~ ρr + ρm + ρΛ
        H ~ √(8π/3 * ρ)
        ℋ ~ a * H
        D(ρr) ~ -4*ρr
        D(ρm) ~ -3*ρm
        ρΛ ~ 3/8π * ΩΛ0
        # perturbation equations
        δρ ~ ρr*δr + ρm*δm
        D(Φ) ~ (-4π/3*a^2/ℋ*δρ - k^2/(3ℋ)*Φ - ℋ*Φ) / ℋ
        D(δr) ~ -4/3*θr/ℋ + 4*D(Φ)
        D(θr) ~ k^2 * (δr/4 + Φ) / ℋ
        D(δm) ~ -θm/ℋ + 3*D(Φ)
        D(θm) ~ -θm + k^2*Φ/ℋ
        Δm ~ δm + 3ℋ*θm/k^2 # gauge-invariant overdensity
    ]
    initial_conditions = [
        ΩΛ0 => 1 - Ωr0 - Ωm0
        ρr => 3/8π * Ωr0 / a^4
        ρm => 3/8π * Ωm0 / a^3
        Φ => 20/15
        δr => -2*Φ
        δm => -3/2*Φ
        θr => 1/2*k^2*a/√(Ωr0)*Φ # τ ≈ a/√(Ωr0) early in radiation domination
        θm => 1/2*k^2*a/√(Ωr0)*Φ
    ]
    M = complete(System(eqs, b, vars, pars; initial_conditions, name = :RMΛ))
    p = Dict(M.Ωr0 => 1e-4, M.Ωm0 => 0.3, M.h => 0.7, M.As => 2e-9, M.ns => 0.96)

    # the default event a ~ 1 is already satisfied where the backwards stage starts, so it never triggers and the whole span is integrated
    probdefault = CosmologyProblem(M, p; tspan = (-8.0, 0.0))
    @test solve(probdefault).bg[1].t[end] == -8.0
    prob = CosmologyProblem(M, p; tspan = (-8.0, 0.0), terminate = nothing)
    @test length(prob.bg) == 1
    ks = 10.0 .^ (0:3)
    sol = solve(prob, ks)
    @test issuccess(sol) && length(sol.bg) == 1

    # the backwards first stage can terminate at an event, which shrinks the span of the later perturbations
    probterm = CosmologyProblem(M, p; tspan = (-8.0, 0.0), terminate = M.a ~ 1e-3)
    solterm = solve(probterm, ks)
    @test solterm.bg[1].t[end] ≈ log(1e-3) && solterm.pts[1].t[begin] ≈ log(1e-3) && solterm.pts[1].t[end] ≈ 0.0

    # splitting all background unknowns off into a first stage leaves the last with none, but gives the same result
    probbg = CosmologyProblem(M, p; bg = ([ρr, ρm], []), tspan = (-8.0, 0.0), terminate = nothing)
    @test isempty(unknowns(probbg.bg[end].f.sys))
    @test solve(probbg, ks)(M.δm, 0.0, ks) ≈ sol(M.δm, 0.0, ks)

    # the backward solve hits its boundary conditions today exactly
    @test sol(M.ρr, 0.0) == 3/8π * p[M.Ωr0]
    @test sol(M.H, 0.0) ≈ 1.0 # since Ωr0 + Ωm0 + ΩΛ0 = 1

    # the stages can be solved one by one
    bgsols = solvebg(prob)
    @test solvept(prob.pt, bgsols, ks)[end].u[end] ≈ sol.pts[end].u[end]

    # the matter power spectrum works without thermodynamics, but CMB spectra need source functions
    @test spectrum_matter(prob, ks) ≈ spectrum_matter(sol, ks) ≈ spectrum_primordial(ks, sol) .* sol(M.Δm, 0.0, ks) .^ 2
    @test_throws Exception spectrum_cmb(:TT, prob, SphericalBesselCache(25:25:100))

    # ForwardDiff should differentiate through the backward-perturbation chain
    diffpars = [M.Ωr0, M.Ωm0]
    probf = remake_function(prob, diffpars)
    tol = 1e-10
    function output(θ)
        solθ = solve(probf(θ), ks; bgopts = (reltol = tol, abstol = tol), ptopts = (reltol = tol, abstol = tol)) # not sol, which would overwrite the outer one
        return [vec(solθ(M.Φ, 0.0, ks)); vec(solθ(M.δm, 0.0, ks))]
    end
    θ0 = [p[par] for par in diffpars]
    J_ad = ForwardDiff.jacobian(output, θ0)
    J_fd = FiniteDiff.finite_difference_jacobian(output, θ0, Val{:central}; relstep = 1e-2, absstep = 1e-6) # absstep keeps Ωr0 = 1e-4 positive; smaller steps are dominated by ODE solver noise
    @test all(isapprox.(J_ad, J_fd; rtol = 1e-3))
end
