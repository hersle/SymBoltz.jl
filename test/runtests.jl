using Test, SymBoltz, Unitful, UnitfulAstro, ModelingToolkit, ForwardDiff, FiniteDiff

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

    @test_warn "Perturbation (mode k = NaN) solution failed" ptsol = solvept(prob.pt, bgsol, [NaN], prob.var2spl; thread = false)
    @test_nowarn ptsol = solvept(prob.pt, bgsol, [1.0], prob.var2spl; thread = false)
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
    @test isapprox(sol(D(M.g.Φ), τ0, ks), sol(D(M.g.Ψ), τ0, ks); atol = 1e-5)
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

    # Test (jl/x^2)(l, x) chain rule
    crazy_x2(l, x) = sin(7*SymBoltz.jl_x2(l, x^2)) # plot(x -> crazy_x2(5, x)) looks very cool
    for l in 2:500
        dcrazy_x2_fd(l, x) = FiniteDiff.finite_difference_derivative(x -> crazy_x2(l, x), x)
        dcrazy_x2_ad(l, x) = ForwardDiff.derivative(x -> crazy_x2(l, x), x)
        @test all(isapprox.(dcrazy_x2_ad.(l, x), dcrazy_x2_fd.(l, x); atol = 1e-6))
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

@testset "Initial conditions" begin
    ks = [1e-1, 1e0] / u"Mpc"
    sol = solve(prob, ks)

    # Check that Fₗ(0) ∝ kˡ
    Fls = sol([M.γ.F0; collect(M.γ.F)], 0.0, ks)
    @test all(Fls[:,1] ./ Fls[:,2] .≈ map(l -> (ks[1]/ks[2])^l, 0:size(Fls)[1]-1))

    # Check initial ratio of metric potentials
    @test all(isapprox.(sol(M.g.Φ / M.g.Ψ, 0.0, ks), sol((1+2/5*M.fν), 0.0); atol = 1e-5))
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
    @test checkvar(M.b.rec.κ̇, 0, 1e-2)
    @test checkvar(M.b.rec.κ, 0, 1e-4)
    @test checkvar(M.b.rec.v, 1e-3, 0)
    @test checkvar(M.b.rec.v̇, 0, 1e1) # TODO: improve
    @test checkvar(M.b.rec.cₛ², 1e-4, 0)
    @test checkvar(M.b.rec.Tb, 0, 1e-5)
    @test checkvar(M.b.rec.Xe, 1e-5, 0)
end

@testset "Whole model without autosplining" begin
    @test mtkcompile(M) isa Any
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
    @test_throws "not an unknown" SymBoltz.mtkcompile_spline(M, [M.g.E])
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
    @test sol(M.g.a, τ0) ≈ sol(M.g.a, τ0, ks) ≈ 1.0
    @test sol(M.χ, τ0) == sol(M.χ, τ0, ks) == 0.0
    @test sol(M.b.rec.κ, τ0) == sol(M.b.rec.κ, τ0, ks) == 0.0
end

@testset "Equal parameters in background and perturbation solutions" begin
    sol = solve(prob, [1.0, 10.0, 100.0])
    pars = [ # choose lots of background parameters that should be equal in perturbations
        M.τ0, M.g.h,
        M.c.Ω₀,
        M.b.Ω₀, M.b.rec.YHe, M.b.rec.fHe, M.b.rec.re1z, M.b.rec.re2z, M.b.rec.κ0,
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

@testset "Dedicated background/perturbation solvers" begin
    bgsol = solvebg(prob.bg) # TODO: @inferred
    @test bgsol isa SymBoltz.ODESolution

    ks = 1.0:1.0:10.0
    ptsol = solvept(prob.pt, bgsol, ks, prob.var2spl) # TODO: @inferred
    @test ptsol isa SymBoltz.EnsembleSolution

    # custom output_func for e.g. source function
    getS = SymBoltz.getsym(prob.pt, M.ST0)
    τ0 = bgsol.t[end]
    τs = range(0.0, τ0, length = 768)
    ks = range(1.0, 1000.0, length = 1000)
    ptsol = solvept(prob.pt, bgsol, ks, prob.var2spl; saveat = τs, output_func = (ptsol, _) -> (getS(ptsol), false))
    Ss = stack(ptsol.u)
    @test size(Ss) == (length(τs), length(ks))
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
    sol = solve(prob, ks_coarse; ptopts = (alg = SymBoltz.KenCarp4(), reltol = 1e-8, abstol = 1e-8, saveat = τs,))
    Sgetter = SymBoltz.getsym(prob.pt, M.ST0)
    Ss = @inferred source_grid(sol, [M.ST0], τs) # TODO: save allocation time with out-of-place version?
    @test Ss == source_grid(prob, [M.ST0], τs, ks_coarse)
    Ss = @inferred source_grid(Ss, ks_coarse, ks_fine) # TODO: save allocation time with out-of-place version?
    Ss = @view Ss[1, :, :]
    Θ0s = @inferred los_integrate(Ss, ls, τs, ks_fine) # TODO: sequential along τ? # TODO: cache kτ0 and x=τ/τ0 (only depends on l)
    P0s = @inferred spectrum_primordial(ks_fine, pars[M.g.h], prob.bg.ps[M.I.As])
    Cls = @inferred spectrum_cmb(Θ0s, Θ0s, P0s, ls, ks_fine)
end

@testset "Background differentiation test" begin
    diffpars = [M.g.h, M.c.Ω₀, M.b.Ω₀, M.γ.T₀, M.ν.Neff, M.h.m_eV, M.b.rec.YHe, M.I.ln_As1e10, M.I.ns]
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

# TODO: why does this fail with ν = nothing??
@testset "RECFAST flags" begin
    for Heflag in [0, 1, 2, 3, 6]
        M = ΛCDM(h = nothing, K = nothing; Heflag)
        pars = parameters_Planck18(M)
        prob = CosmologyProblem(M, pars)
        #bgopts = (alg = SymBoltz.Rodas4P(), reltol = 1e-10, abstol = 1e-10)
        sol = solve(prob)
        println(Heflag, ": ", issuccess(sol))
        #@test all(sol[M.b.rec.XH⁺] .≤ 1.0)
        #@test all(sol[M.b.rec.XHe⁺] .≤ 1.0)
        #@test all(sol[M.b.rec.XHe⁺⁺/M.b.rec.fHe] .≤ 1.0)
    end
end
