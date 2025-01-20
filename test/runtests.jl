using Test
using SymBoltz

M = ΛCDM()
pars = SymBoltz.parameters_Planck18(M)
prob = CosmologyProblem(M, pars, Dict(M.Λ.Ω₀ => 0.5), [M.g.ℰ ~ 1])

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

