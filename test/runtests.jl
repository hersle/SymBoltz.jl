using Test
using SymBoltz

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
