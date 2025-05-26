using Test, SymBoltz, Unitful, UnitfulAstro, ModelingToolkit, ForwardDiff, FiniteDiff

M = ΛCDM(K = nothing) # flat
pars = SymBoltz.parameters_Planck18(M)
prob = CosmologyProblem(M, pars)

@testset "Solution accessing" begin
    ks = [1e2, 1e3]
    τs = [1.0, 2.0, 3.0]
    is = [M.g.a, M.g.a, M.g.a, M.g.a]
    nk, nt, ni = length(ks), length(τs), length(is)

    # Size of solution output should match input arguments in order
    sol = solve(prob)
    @test size(sol(τs[1], is[1])) == () # background
    @test size(sol(τs[1], is)) == (ni,)
    @test size(sol(τs, is[1])) == (nt,)
    @test size(sol(τs, is)) == (nt, ni)
    @test_throws "No perturbations" sol(ks, τs, is)
    #@test_throws "below minimum solved time" sol(sol[M.τ][begin]-1, is)
    #@test_throws "above maximum solved time" sol(sol[M.τ][end]+1, is)

    sol = solve(prob, ks)
    @test size(sol(ks[1], τs[1], is[1])) == () # perturbations
    @test size(sol(ks[1], τs[1], is)) == (ni,)
    @test size(sol(ks[1], τs, is[1])) == (nt,)
    @test size(sol(ks[1], τs, is)) == (nt, ni)
    @test size(sol(ks, τs[1], is[1])) == (nk,)
    @test size(sol(ks, τs[1], is)) == (nk, ni)
    @test size(sol(ks, τs, is[1])) == (nk, nt)
    @test size(sol(ks, τs, is)) == (nk, nt, ni)
    @test_throws "below minimum solved wavenumber" sol(ks[begin]-1, τs, is)
    @test_throws "above maximum solved wavenumber" sol(ks[end]+1, τs, is)
    #@test_throws "below minimum solved time" sol(ks[1], sol[M.τ][begin]-1, is)
    #@test_throws "above maximum solved time" sol(ks[1], sol[M.τ][end]+1, is)

    # TODO: also test array indexing
end

@testset "Accessing derivative variables" begin
    ks = 1.0 / u"Mpc"
    sol = solve(prob, ks)
    D = Differential(M.τ)
    τ0 = sol[M.τ0]
    @test isapprox(sol(τ0, D(M.g.a)), 1.0; atol = 1e-5)
    @test isapprox(sol(τ0, D(M.g.z)), -1.0; atol = 1e-5)
    @test_throws "Could not express derivative" sol(τ0, D(D(M.g.z)))
    sol(ks, τ0, D(M.g.Φ))
    @test isapprox(sol(ks, τ0, D(M.g.Φ)), sol(ks, τ0, D(M.g.Ψ)); atol = 1e-5)
    @test_throws "Could not express derivative" sol(ks, τ0, D(D(M.g.Φ)))
    @test_throws "Could not express derivative" sol(ks, τ0, D(D(M.g.Ψ)))
end

@testset "Solution interpolation" begin
    ks = 10 .^ range(-5, 1, length=100) / u"Mpc"
    sol = solve(prob, ks)
    ks = range(extrema(ks)..., length=500)
    τs = range(extrema(sol.bg.t)..., length=500)
    is = [M.g.a, M.G.ρ, M.g.Φ, M.g.Ψ]
    @test sol(ks, τs, is; smart = true) == sol(ks, τs, is; smart = false)
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

        # Test jₗ(x) / x²
        SymBoltz.jl_x2!(jlfast, l, x) # unsafe implementation
        if x == 0.0
            jlslow[3] = 1/15 # l = 2
            jlslow[4:1000] .= 0.0 # l ≥ 3
        else
            SymBoltz.jlsafe!(jlslow, l, x) # safe implementation
            jlslow ./= x^2
            @test jlfast[1:2] ≈ jlslow[1:2] # test l = 0 and 1 only for x > 0 (diverges for x = 0)
        end
        @test jlfast[3:end] ≈ jlslow[3:end] # l ≥ 2
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
    crazy(l, x) = sin(7*SymBoltz.jl(l, x^2)) # crazy composite function involving jl

    # Test (jl/x^2)(l, x) chain rule
    crazy(l, x) = sin(7*SymBoltz.jl_x2(l, x^2)) # plot(x -> crazy(5, x)) looks very cool
    for l in 2:500
        dcrazy_fd(l, x) = FiniteDiff.finite_difference_derivative(x -> crazy(l, x), x)
        dcrazy_ad(l, x) = ForwardDiff.derivative(x -> crazy(l, x), x)
        @test all(isapprox.(dcrazy_ad.(l, x), dcrazy_fd.(l, x); atol = 1e-6))
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
    zs = sol(τs, M.g.z)
    as = sol(τs, M.g.a)
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
    ks = [1e-1, 1e0] / u"Mpc"
    sol = solve(prob, ks)

    # Check that Fₗ(0) ∝ kˡ
    Fls = sol(ks, 0.0, collect(M.γ.F))
    @test all(Fls[1,:] ./ Fls[2,:] .≈ map(l -> (ks[1]/ks[2])^l, 0:length(M.γ.F)-1))

    # Check initial ratio of metric potentials
    sol(ks, 0.0, M.g.Φ / M.g.Ψ) .≈ sol(0.0, (1+2/5*M.fν))
end

@testset "Automatic background/thermodynamics splining" begin
    sol = solve(prob, 1.0) # solve with one perturbation mode to activate splining
    τs = SymBoltz.timeseries.(sol, log10(M.g.a), range(-8, 0, length=100)) # TODO a => as syntax
    tests = [
        (M.g.a, 1e-6, 0)
        (M.b.rec.κ̇, 0, 1e-2)
        (M.b.rec.κ, 0, 1e-4)
        (M.b.rec.v, 1e-3, 0)
        (M.b.rec.v̇, 0, 1e-0) # TODO: improve
        (M.b.rec.cₛ², 1e-4, 0)
        (M.b.rec.Tb, 0, 1e-5)
        (M.b.rec.Xe, 1e-6, 0)
    ]
    for (var, atol, rtol) in tests
        vals1 = sol(τs, var) # from background
        vals2 = sol(1.0, τs, var) # from splined perturbations
        @test all(isapprox.(vals1, vals2; atol, rtol))
    end
end

@testset "Whole model without autosplining" begin
    @test structural_simplify(M) isa Any
end

#=
@testset "Spline observed variables" begin
    prob1 = prob
    prob2 = CosmologyProblem(M, pars; spline = [unknowns(prob1.bg.f.sys); M.b.rec.v])
    @test haskey(prob2.var2spl, M.b.rec.v)

    obs2 = Dict(string(eq.lhs) => string(eq.rhs) for eq in observed(prob2.pt.f.sys)) # for easy lookup based on LHS
    @test obs2["b₊rec₊v(τ)"] == "SymBoltz.value(b₊rec₊v_spline, τ)"
    @test obs2["b₊rec₊vˍτ(τ)"] == "SymBoltz.derivative(b₊rec₊v_spline, τ, 1)"

    sol1 = solve(prob1, 1.0)
    sol2 = solve(prob2, 1.0)
    τs = sol1[M.τ]
    v1s = sol1(1.0, τs, M.b.rec.v)
    v2s = sol2(1.0, τs, M.b.rec.v)
    @test all(isapprox.(v1s, v2s; atol = 1e-3))
end
=#
@testset "Do not spline observed variables" begin
    @test_throws "not an unknown" SymBoltz.structural_simplify_spline(M, [M.g.E])
end

@testset "Primordial power spectrum pivot scale" begin
    h = pars[M.g.h]
    k = 0.05 / u"Mpc" # ≠ 0.05/(Mpc/h)
    sol = solve(prob, k)
    @test sol[M.I.kpivot] ≈ SymBoltz.k_dimensionless(k, h)
end

@testset "Time and optical depth today" begin
    ks = 1.0
    prob = CosmologyProblem(M, pars) # recreate since solution usually modifies problem parameters
    sol = solve(prob, ks)
    τ0 = sol[M.τ0]
    @test sol(τ0, M.g.a) ≈ sol(ks, τ0, M.g.a) ≈ 1.0
    @test sol(τ0, M.χ) == sol(ks, τ0, M.χ) == 0.0
    @test sol(τ0, M.b.rec.κ) == sol(ks, τ0, M.b.rec.κ) == 0.0
end

@testset "Equal parameters in background and perturbation solutions" begin
    sol = solve(prob, [1.0, 10.0, 100.0])
    pars = [ # choose lots of background parameters that should be equal in perturbations
        M.τ0, M.g.h,
        M.c.Ω₀,
        M.b.Ω₀, M.b.rec.Yp, M.b.rec.fHe, M.b.rec.re1z, M.b.rec.re2z, M.b.rec.κ0,
        M.γ.Ω₀, M.γ.T₀,
        M.ν.Ω₀, M.ν.T₀, M.ν.Neff,
        M.h.Ω₀, M.h.T₀, M.h.m, M.h.y₀, M.h.Iρ₀,
        M.Λ.Ω₀,
        M.I.As, M.I.kpivot, M.I.ns
    ]
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
    function logDlTT(logθ)
        θ = exp.(logθ)
        prob′ = probgen(θ)
        DlTT = spectrum_cmb(:TT, prob′, l; normalization = :Dl)
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
    prob1 = probgen(2.7)
    vals = getter(prob1)
    @test vals[1] == 2.7
    @test isfinite(vals[2])
    @test vals[3] == 1.0
    @test issuccess(solve(prob, 1.0))
end

@testset "Parameter updater and remake" begin
    probgen = SymBoltz.parameter_updater(prob, [M.c.Ω₀])

    newprob = probgen([0.3])
    @test newprob.bg.ps[M.c.Ω₀] == newprob.pt.ps[M.c.Ω₀] == 0.3
    @test newprob.bg.ps[M.γ.Ω₀ + M.ν.Ω₀ + M.h.Ω₀ + M.b.Ω₀ + M.c.Ω₀ + M.Λ.Ω₀] == newprob.pt.ps[M.γ.Ω₀ + M.ν.Ω₀ + M.h.Ω₀ + M.b.Ω₀ + M.c.Ω₀ + M.Λ.Ω₀] == 1.0

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

# TODO: test for spectrum_cmb etc.
# TODO: define closure based on this?
@testset "Type stability" begin
    getτ0 = SymBoltz.getsym(prob.bg, M.τ0)
    sol = solve(prob; save_everystep = false, save_start = false, save_end = true)
    τ0 = @inferred getτ0(sol.bg)
    @test sol.bg.t[begin] == sol.bg.t[end] == τ0

    ls = 20:20:500
    kτ0s_coarse, kτ0s_fine = SymBoltz.cmb_kτ0s(ls[begin], ls[end])
    ks_coarse, ks_fine = kτ0s_coarse / τ0, kτ0s_fine / τ0
    ks_fine = clamp.(ks_fine, ks_coarse[begin], ks_coarse[end]) # numerics can change endpoint slightly
    τs = range(0.0, τ0, length = 768)
    #sol = solve(prob, ks_coarse)
    #Ss = sol(ks_fine, τs, M.ST0)
    sol = solve(prob, ks_coarse; ptopts = (reltol = 1e-8, saveat = τs,))
    Sgetter = SymBoltz.getsym(prob.pt, M.ST0)
    Ss = @inferred source_grid(sol, Sgetter, sol.ks, ks_fine, τs) # TODO: could save allocation time with out-of-place version
    Θ0s = @inferred los_integrate(Ss, ls, ks_fine, τs) # TODO: sequential along τ? # TODO: cache kτ0 and x=τ/τ0 (only depends on l)
    P0s = @inferred spectrum_primordial(ks_fine, pars[M.g.h], pars[M.I.As])
    Cls = @inferred spectrum_cmb(Θ0s, Θ0s, P0s, ls, ks_fine)
end
