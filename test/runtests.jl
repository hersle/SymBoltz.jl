using Test, SymBoltz, Unitful, UnitfulAstro

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

    ts1 = SymBoltz.timeseries.(sol, M.g.z, [0, 9, 99]) # invert z to t with root finding
    ts2 = SymBoltz.timeseries(sol, M.g.z, [0, 9, 99]; Nextra=9) # invert z to t with splining (need points for accuracy)
    @test all(ts1 .≈ ts2)
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
    ts = SymBoltz.timeseries.(sol, log10(M.g.a), range(-8, 0, length=100))
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
    @test_nowarn structural_simplify(M)
end
