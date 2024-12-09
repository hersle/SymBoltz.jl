function parameters_Planck18(M::CosmologyModel)
    # https://arxiv.org/pdf/1807.06209, Table 5
    h = 0.6736
    params = Dict(
        M.g.h => h,
        M.γ.T₀ => 2.7255,
        M.c.Ω₀ => 0.1200 / h^2,
        M.b.Ω₀ => 0.0224 / h^2,
        M.b.rec.Yp => 0.2454,
        M.ν.Neff => 2.99 # TODO: vs 3.046 or 3.044?

        # TODO: how to handle
        # 1) backwards with default for Λ.Ω₀? --> skip E = 1
        # 2) forwards with no default?        --> shooting
        #M.Λ.Ω => 1 - sum(s.Ω₀ for s in [M.γ, M.ν, M.h, M.c, M.b])
    )
    have(M.sys, :h) && push!(params, M.h.m => 0.06 * eV/c^2) # TODO: ?
    return params
end
