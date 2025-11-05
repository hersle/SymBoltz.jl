using Test, SymBoltz, Unitful, UnitfulAstro, ModelingToolkit, ForwardDiff, FiniteDiff, BenchmarkTools

M = ΛCDM(K = nothing) # flat
pars = SymBoltz.parameters_Planck18(M)
prob = CosmologyProblem(M, pars)

# Must come first because warnings are only given once
@testset "Solve failure warnings" begin
    Ωc0 = prob.bg.ps[M.c.Ω₀]
    prob.bg.ps[M.c.Ω₀] = NaN # bad
    bgsol = @test_warn "Background solution failed" solvebg(prob.bg)
    prob.bg.ps[M.c.Ω₀] = Ωc0 # restore good
    bgsol = @test_nowarn solvebg(prob.bg)

    @test_warn "Perturbation (mode k = NaN) solution failed" ptsol = solvept(prob.pt, bgsol, [NaN], prob.bgspline; thread = false)
    @test_nowarn ptsol = solvept(prob.pt, bgsol, [1.0], prob.bgspline; thread = false)
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

@testset "Accessing derivative variables" begin
    ks = 1.0 / u"Mpc"
    sol = solve(prob, ks)
    D = Differential(M.τ)
    τ0 = sol[M.τ0]
    @test isapprox(sol(D(M.g.a), τ0), 1.0; atol = 1e-4)
    @test isapprox(sol(D(M.g.z), τ0), -1.0; atol = 1e-4)
    @test_throws "Could not express derivative" sol(D(D(M.g.z)), τ0)
    sol(D(M.g.Φ), τ0, ks)
    @test isapprox(sol(D(M.g.Φ), τ0, ks), sol(D(M.g.Ψ), τ0, ks); atol = 2e-5)
    @test_throws "Could not express derivative" sol(D(D(M.g.Φ)), τ0, ks)
    @test_throws "Could not express derivative" sol(D(D(M.g.Ψ)), τ0, ks)
end

@testset "Solution interpolation" begin
    ks = 10 .^ range(-5, 1, length=100) / u"Mpc"
    sol = solve(prob, ks)
    ks = range(extrema(ks)..., length=500)
    τs = range(extrema(sol.bg.t)..., length=500)
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
    jl = SphericalBesselCache(ls)

    @test_throws BoundsError jl(5, 0.0) # not cached
    @test_throws BoundsError jl(10, -1.0)
    @test_throws BoundsError jl(10, jl.x[end] + 1.0)
    @test jl(10, 0.0) == SymBoltz.sphericalbesselj(10, 0.0)
    @test jl(10, jl.x[end]) == SymBoltz.sphericalbesselj(10, jl.x[end])
    @test isapprox(jl(10, 123.456), SymBoltz.sphericalbesselj(10, 123.456); atol = 1e-3)

    t1 = @belapsed $jl(10, π)
    t2 = @belapsed SymBoltz.sphericalbesselj(10, π)
    @test t1 < t2 # faster
    @test (@ballocated $jl(10, π)) == 0 # non-allocating

    j10(x) = jl(10, x)
    @test isfinite(ForwardDiff.derivative(j10, π))
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

@testset "Initial conditions" begin
    τini = prob.bg.tspan[1]
    ks = [1e-1, 1e0] / u"Mpc"
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
    @test all(Fls[:,1] ./ Fls[:,2] .≈ map(l -> (ks[1]/ks[2])^l, 0:size(Fls)[1]-1))

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

@testset "Whole model without autosplining" begin
    @test mtkcompile(M) isa Any
end

# TODO: optionally spline observables, too, with cubic hermite splines using analytic derivatives, to save computations in perturbations (e.g. visibility function and derivatives for CMB)
@testset "Do not spline observed variables" begin
    @test_throws "not an unknown" SymBoltz.mtkcompile_spline(M, [M.g.H])
end

@testset "Wavenumber units and primordial power spectrum pivot scale" begin
    h = pars[M.g.h]
    k = 0.05 / u"Mpc" # ≠ 0.05/(Mpc/h)
    sol = solve(prob, k)
    @test sol[M.I.kpivot] ≈ SymBoltz.k_dimensionless(k, sol.bg)

    ks = 1.0:100.0
    sol = solve(prob, ks)
    P1 = sol(M.I.P, sol.bg.t[begin], ks)
    P2 = spectrum_primordial(ks, sol)
    @test all(isapprox.(P1, P2))

    ks = collect(1.0:100.0) / u"Mpc"
    sol = solve(prob, ks)
    P1 = sol(M.I.P, sol.bg.t[begin], ks)
    P2 = spectrum_primordial(SymBoltz.k_dimensionless.(ks, h), sol)
    @test all(isapprox.(P1, P2))
end

@testset "Time and optical depth today" begin
    ks = 1.0
    prob = CosmologyProblem(M, pars) # recreate since solution usually modifies problem parameters
    sol = solve(prob, ks)
    τ0 = sol[M.τ0]
    @test sol(M.g.a, τ0) ≈ sol(M.g.a, τ0, ks) ≈ 1.0
    @test sol(M.χ, τ0) == sol(M.χ, τ0, ks) == 0.0
    @test sol(M.b.κ, τ0) == sol(M.b.κ, τ0, ks) == 0.0
end

@testset "Equal parameters in background and perturbation solutions" begin
    sol = solve(prob, [1.0, 10.0, 100.0])
    pars = [ # choose lots of background parameters that should be equal in perturbations
        M.τ0, M.g.h,
        M.c.Ω₀,
        M.b.Ω₀, M.b.YHe, M.b.fHe, M.b.κ0,
        M.γ.Ω₀, M.γ.T₀,
        M.ν.Ω₀, M.ν.T₀, M.ν.Neff,
        M.h.Ω₀, M.h.T₀, M.h.m, M.h.y₀, M.h.Iρ₀,
        M.Λ.Ω₀,
        M.I.As, M.I.kpivot, M.I.ns
    ]
    @test allequal([extrema(sol.bg.t); map(pt -> extrema(pt.t), sol.pts)]) # background and perturbation should have equal timespans
    @test all(allequal([sol.bg.ps[par]; map(pt -> pt.ps[par], sol.pts)]) for par in pars)
end

@testset "Success checking" begin
    @test issuccess(solve(prob, 1.0))
    @test !issuccess(solve(prob, 0.0))
end

@testset "Consistent AD and FD derivatives of power spectra" begin
    k = 10 .^ range(-3, 0; length = 20) / u"Mpc"
    diffpars = [M.c.Ω₀, M.b.Ω₀] # TODO: h, ...
    probgen = SymBoltz.parameter_updater(prob, diffpars)
    function logP(logθ)
        θ = exp.(logθ)
        prob′ = probgen(θ)
        P = spectrum_matter(prob′, k)
        return log.(P / u"Mpc^3")
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
        lines!(ax, log10.(k*u"Mpc"), ∂logP_∂logθ_ad[:, i]; color, alpha, linestyle = :solid)
        lines!(ax, log10.(k*u"Mpc"), ∂logP_∂logθ_fd[:, i]; color, alpha, linestyle = :dash)
    end
    fig
    =#

    l = 25:25:1000
    jl = SphericalBesselCache(l)
    function logDlTT(logθ)
        θ = exp.(logθ)
        prob′ = probgen(θ)
        DlTT = spectrum_cmb(:TT, prob′, jl; normalization = :Dl)
        return log.(DlTT)
    end
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

    probgen = SymBoltz.parameter_updater(prob0, M.γ.T₀)
    prob1 = probgen(2.73)
    vals = getter(prob1)
    @test vals[1] == 2.73
    @test isfinite(vals[2])
    @test vals[3] == 1.0
    @test issuccess(solve(prob1, 1.0))
end

@testset "Parameter updater and remake" begin
    probgen = SymBoltz.parameter_updater(prob, [M.c.Ω₀])

    newprob = probgen([0.3])
    @test newprob.bg.ps[M.c.Ω₀] == newprob.pt.ps[M.c.Ω₀] == 0.3
    @test newprob.bg.ps[M.γ.Ω₀ + M.ν.Ω₀ + M.h.Ω₀ + M.b.Ω₀ + M.c.Ω₀ + M.Λ.Ω₀] == newprob.pt.ps[M.γ.Ω₀ + M.ν.Ω₀ + M.h.Ω₀ + M.b.Ω₀ + M.c.Ω₀ + M.Λ.Ω₀] ≈ 1.0

    ks = 10 .^ range(0, 3, length=10) # faster than with u"Mpc" # TODO: investigate further: Unitful is very slow with autodiff?
    sol = solve(newprob, ks)
    @test all(map(SymBoltz.successful_retcode, sol.pts))

    function Pk(Ωc0)
        newprob = probgen([Ωc0])
        return spectrum_matter(newprob, ks)
    end
    isnonzero(x) = isfinite(x) && !iszero(x)
    @test all(isnonzero.(Pk(0.3)))
    @test all(isnonzero.(Pk(0.3)))
    @test all(isnonzero.(ForwardDiff.derivative(Pk, 0.3)))
    @test all(isnonzero.(ForwardDiff.derivative(Pk, 0.3))) # twice (successive calls failed with earlier bug)
end

@testset "Dedicated background/perturbation solvers" begin
    bgsol = solvebg(prob.bg) # TODO: @inferred
    @test bgsol isa SymBoltz.ODESolution

    ks = 1.0:1.0:10.0
    ptsol = solvept(prob.pt, bgsol, ks, prob.bgspline) # TODO: @inferred
    @test ptsol isa SymBoltz.EnsembleSolution

    # custom output_func for e.g. source function
    getS = SymBoltz.getsym(prob.pt, M.ST0)
    τi = bgsol.t[begin]
    τ0 = bgsol.t[end]
    τs = range(τi, τ0, length = 768)
    ks = range(1.0, 1000.0, length = 1000)
    ptsol = solvept(prob.pt, bgsol, ks, prob.bgspline; saveat = τs, output_func = (ptsol, _) -> (getS(ptsol), false))
    Ss = stack(ptsol.u)
    @test size(Ss) == (length(τs), length(ks))
end

# TODO: test for spectrum_cmb etc.
# TODO: define closure based on this?
@testset "Type stability" begin
    getτ0 = SymBoltz.getsym(prob.bg, M.τ0)
    sol = solve(prob; save_everystep = false, save_start = true, save_end = true)
    τi = sol.bg.t[begin]
    τ0 = @inferred getτ0(sol.bg)
    @test sol.bg.t[end] == τ0

    ls = 20:20:500
    jl = SphericalBesselCache(ls)
    kτ0s_coarse, kτ0s_fine = SymBoltz.cmb_kτ0s(ls[begin], ls[end])
    ks_coarse, ks_fine = kτ0s_coarse / τ0, kτ0s_fine / τ0
    ks_fine = clamp.(ks_fine, ks_coarse[begin], ks_coarse[end]) # numerics can change endpoint slightly
    τs = range(τi, τ0, length = 768)
    #sol = solve(prob, ks_coarse)
    #Ss = sol(ks_fine, τs, M.ST0)
    ptopts = (alg = SymBoltz.DEFAULT_PTALG, reltol = 1e-8, abstol = 1e-8)
    sol = solve(prob, ks_coarse; ptopts = (ptopts..., saveat = τs))
    Sgetter = SymBoltz.getsym(prob.pt, M.ST0)
    Ss = @inferred source_grid(sol, [M.ST0], τs) # TODO: save allocation time with out-of-place version?
    @test Ss == source_grid(prob, [M.ST0], τs, ks_coarse; ptopts)
    Ss = @inferred source_grid(Ss, ks_coarse, ks_fine) # TODO: save allocation time with out-of-place version?
    Ss = @view Ss[1, :, :]
    Θ0s = @inferred los_integrate(Ss, ls, τs, ks_fine, jl) # TODO: sequential along τ? # TODO: cache kτ0 and x=τ/τ0 (only depends on l)
    P0s = @inferred spectrum_primordial(ks_fine, pars[M.g.h], prob.bg.ps[M.I.As])
    Cls = @inferred spectrum_cmb(Θ0s, Θ0s, P0s, ls, ks_fine)
end

@testset "Background differentiation test" begin
    diffpars = [M.g.h, M.c.Ω₀, M.b.Ω₀, M.γ.T₀, M.ν.Neff, M.h.m_eV, M.b.YHe, M.I.ln_As1e10, M.I.ns]
    probgen = SymBoltz.parameter_updater(prob, diffpars)
    getτ0 = SymBoltz.getsym(prob, M.τ0)
    τ0(θ) = getτ0(solve(probgen(θ)))
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
function stability(M::System, ks, vary::Dict, nsamples; verbose = false, kwargs...)
    prob0 = CosmologyProblem(M, Dict(keys(vary) .=> NaN))
    pars = collect(keys(vary))
    probgen = parameter_updater(prob0, pars)
    lo = [bound[1] for bound in values(vary)] # lower corner of parameter space
    hi = [bound[2] for bound in values(vary)] # uppper corner of parameter space
    samples = QuasiMonteCarlo.sample(nsamples, lo, hi, LatinHypercubeSample())
    nsuccess = 0
    verbose && println("Varying ", keys(vary))
    for sample in eachcol(samples)
        prob = probgen(sample)
        sol = solve(prob, ks; verbose, kwargs...)
        if issuccess(sol)
            nsuccess += 1
        end
        verbose && println(issuccess(sol) ? "PASS" : "FAIL", ": ", sample)
    end
    return nsuccess / nsamples
end
vary = Dict(par => (0.5val, 1.5val) for (par, val) in pars) # ± 50% around fiducial values
ks = [1e0, 1e1, 1e2, 1e3]
@testset "Stability of problems throughout parameter space with Latin hypercube sampling" begin
    @test stability(M, ks, vary, 100; verbose = true) == 1.0 # 100%

    M1 = ΛCDM(K = nothing; Hswitch = 0)
    @test stability(M1, ks, vary, 100; verbose = true) == 1.0

    M2 = ΛCDM(K = nothing; Heswitch = 0)
    @test stability(M2, ks, vary, 100; verbose = true) == 1.0

    M3 = ΛCDM(K = nothing; reionization = false)
    @test stability(M3, ks, vary, 100; verbose = true) == 1.0
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
    jl = SphericalBesselCache(20:20:3000)
    DlTT = spectrum_cmb(:TT, prob, jl; normalization = :Dl)
    DlEE = spectrum_cmb(:EE, prob, jl; normalization = :Dl)
end
