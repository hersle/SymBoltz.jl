function parameters_Planck18(M::System)
    # https://arxiv.org/pdf/1807.06209#table.5
    h = 0.6736
    Neff = 3.046 # total effective number of neutrino species
    Nh = 3 # number of degenerate massive neutrinos
    params = Dict(
        M.g.h => h,
        M.γ.T₀ => 2.7255,
        M.c.Ω₀ => 0.1200 / h^2,
        M.b.Ω₀ => 0.0224 / h^2,
        M.b.YHe => 0.2454,
    )
    have(M, :h) && merge!(params, Dict(
        M.h.m_eV => 0.06 / Nh, # mass per neutrino
        M.h.N => Nh,
    ))
    have(M, :ν) && push!(params, M.ν.Neff => Neff - (have(M, :h) ? Nh : 0)) # each massive neutrino contributes 1 to Neff (since Tν = (4/11)^(1/3)*Tγ) # TODO: rename M.ν.Neff to M.ν.N, similar to CLASS' Nur name
    have(M, :I) && merge!(params, Dict(
        M.I.ln_As1e10 => log(2.099e-9*1e10),
        M.I.ns => 0.965
    ))
    return params
end
