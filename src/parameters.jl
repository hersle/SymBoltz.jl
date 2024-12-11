function parameters_Planck18(M::CosmologyModel)
    # https://arxiv.org/pdf/1807.06209, Table 5
    h = 0.6736
    params = Dict(
        M.g.h => h,
        M.γ.T₀ => 2.7255,
        M.c.Ω₀ => 0.1200 / h^2,
        M.b.Ω₀ => 0.0224 / h^2,
        M.b.rec.Yp => 0.2454,
        M.ν.Neff => 2.99
    )
    have(M.sys, :h) && push!(params, M.h.m => 0.06 * eV/c^2)
    return params
end
