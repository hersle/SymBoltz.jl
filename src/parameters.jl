function parameters_Planck18(M::System)
    # https://arxiv.org/pdf/1807.06209#table.5
    h = 0.6736
    params = Dict(
        M.g.h => h,
        M.γ.T₀ => 2.7255,
        M.c.Ω₀ => 0.1200 / h^2,
        M.b.Ω₀ => 0.0224 / h^2,
        M.b.rec.YHe => 0.2454,
    )
    have(M, :ν) && push!(params, M.ν.Neff => 2.99)
    have(M, :h) && push!(params, M.h.m_eV => 0.06)
    have(M, :I) && merge!(params, Dict(
        M.I.ln_As1e10 => log(2.099e-9*1e10),
        M.I.ns => 0.965
    ))
    return params
end
