using Test, SymBoltz, Unitful, UnitfulAstro, ModelingToolkit

M = ΛCDM(K = nothing) # flat
pars = SymBoltz.parameters_Planck18(M)
prob = CosmologyProblem(M, pars)

@testset "Solution accessing" begin
    ks = [1e2, 1e3]
    ts = [1.0, 2.0, 3.0]
    is = [M.g.a, M.g.a, M.g.a, M.g.a]
    nk, nt, ni = length(ks), length(ts), length(is)

    # Size of solution output should match input arguments in order
    sol = solve(prob)
    @test size(sol(ts[1], is[1])) == () # background
    @test size(sol(ts[1], is)) == (ni,)
    @test size(sol(ts, is[1])) == (nt,)
    @test size(sol(ts, is)) == (nt, ni)
    @test_throws "No perturbations" sol(ks, ts, is)
    #@test_throws "below minimum solved time" sol(sol[M.t][begin]-1, is)
    #@test_throws "above maximum solved time" sol(sol[M.t][end]+1, is)

    sol = solve(prob, ks)
    @test size(sol(ks[1], ts[1], is[1])) == () # perturbations
    @test size(sol(ks[1], ts[1], is)) == (ni,)
    @test size(sol(ks[1], ts, is[1])) == (nt,)
    @test size(sol(ks[1], ts, is)) == (nt, ni)
    @test size(sol(ks, ts[1], is[1])) == (nk,)
    @test size(sol(ks, ts[1], is)) == (nk, ni)
    @test size(sol(ks, ts, is[1])) == (nk, nt)
    @test size(sol(ks, ts, is)) == (nk, nt, ni)
    @test_throws "below minimum solved wavenumber" sol(ks[begin]-1, ts, is)
    @test_throws "above maximum solved wavenumber" sol(ks[end]+1, ts, is)
    #@test_throws "below minimum solved time" sol(ks[1], sol[M.t][begin]-1, is)
    #@test_throws "above maximum solved time" sol(ks[1], sol[M.t][end]+1, is)

    # TODO: also test array indexing
end

@testset "Accessing derivative variables" begin
    ks = 1.0 / u"Mpc"
    sol = solve(prob, ks)
    D = Differential(M.t)
    t0 = time_today(sol)
    @test isapprox(sol(t0, D(M.g.a)), 1.0; atol = 1e-5)
    @test isapprox(sol(t0, D(M.g.z)), -1.0; atol = 1e-5)
    @test_throws "Could not express derivative" sol(t0, D(D(M.g.z)))
    sol(ks, t0, D(M.g.Φ))
    @test isapprox(sol(ks, t0, D(M.g.Φ)), sol(ks, t0, D(M.g.Ψ)); atol = 1e-5)
    @test_throws "Could not express derivative" sol(ks, t0, D(D(M.g.Φ)))
    @test_throws "Could not express derivative" sol(ks, t0, D(D(M.g.Ψ)))
end

@testset "Solution interpolation" begin
    ks = 10 .^ range(-5, 1, length=100) / u"Mpc"
    sol = solve(prob, ks)
    ks = range(extrema(ks)..., length=500)
    ts = range(extrema(sol.bg.t)..., length=500)
    is = [M.g.a, M.G.ρ, M.g.Φ, M.g.Ψ]
    @test sol(ks, ts, is; smart = true) == sol(ks, ts, is; smart = false)
end

@testset "Spherical Bessel function" begin
    l = 0:1000
    jlfast = zeros(size(l))
    jlslow = zeros(size(l))
    for x in 0.0:0.001:10.0 # near 0 is most sketchy
        SymBoltz.sphericalbesseljfast!(jlfast, l, x) # unsafe implementation
        SymBoltz.sphericalbesseljslow!(jlslow, l, x) # safe implementation
        @test jlfast ≈ jlslow
    end
end

@testset "Extend array" begin
    @test_throws "Cannot extend empty array" SymBoltz.extend_array([], 0)
    @test SymBoltz.extend_array(1.0:1.0:3.0, 0) == 1.0:1.0:3.0
    @test SymBoltz.extend_array(1.0:1.0:3.0, 4) == 1.0:0.2:3.0
end

@testset "Timeseries" begin
    sol = solve(prob)
    ts = SymBoltz.timeseries(sol; Nextra=1) # naive implementation could transform endpoints slightly through exp(log(t))
    @test sol(ts, M.g.a) isa AbstractArray # should be in bounds

    # First compute z at known t
    ts = range(extrema(ts)..., length=500)
    zs = sol(ts, M.g.z)

    # 1) Invert z to t with root finding
    ts1 = SymBoltz.timeseries(sol, M.g.z, zs)
    @test all(isapprox.(ts, ts1; atol = 1e-12))

    # 2) Invert z and ż to t with Hermite spline
    ts2 = SymBoltz.timeseries(sol, M.g.z, M.g.ż, zs)
    @test all(isapprox.(ts, ts2; atol = 1e-6))
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
    ts = SymBoltz.timeseries.(sol, log10(M.g.a), range(-8, 0, length=100)) # TODO a => as syntax
    tests = [
        (M.g.a, 1e-6, 0)
        (M.b.rec.τ̇, 0, 1e-2)
        (M.b.rec.τ, 0, 1e-4)
        (M.b.rec.v, 1e-3, 0)
        (M.b.rec.v̇, 0, 1e-0) # TODO: improve
        (M.b.rec.cₛ², 1e-4, 0)
        (M.b.rec.Tb, 0, 1e-5)
        (M.b.rec.Xe, 1e-6, 0)
    ]
    for (var, atol, rtol) in tests
        vals1 = sol(ts, var) # from background
        vals2 = sol(1.0, ts, var) # from splined perturbations
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
    @test obs2["b₊rec₊v(t)"] == "SymBoltz.value(b₊rec₊v_spline, t)"
    @test obs2["b₊rec₊vˍt(t)"] == "SymBoltz.derivative(b₊rec₊v_spline, t, 1)"

    sol1 = solve(prob1, 1.0)
    sol2 = solve(prob2, 1.0)
    ts = sol1[M.t]
    v1s = sol1(1.0, ts, M.b.rec.v)
    v2s = sol2(1.0, ts, M.b.rec.v)
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
    sol = solve(prob)
    t0 = time_today(sol)
    @test sol(t0, M.g.a) ≈ 1.0
    @test sol(t0, M.b.rec.τ) == 0.0
end

@testset "Success checking" begin
    @test issuccess(solve(prob, 1.0))
    @test !issuccess(solve(prob, 0.0))
end
