"""
    photons(g; polarization = true, lmax = 10, name = :γ, kwargs...)

Create a particle species for photons in the spacetime with metric `g`.
"""
function photons(g; polarization = true, lmax = 10, name = :γ, kwargs...)
    lmax >= 3 || error("Need lmax >= 3")
    description = "Photons"
    γ = radiation(g; adiabatic = true, name, description, kwargs...) |> background |> complete # prevent namespacing in extension below

    vars = @variables begin
        F0(τ, k), [description = "Distribution function monopole"]
        F(τ, k)[1:lmax], [description = "Distribution function multipoles"]
        Θ0(τ, k), [description = "Temperature perturbation monopole"]
        Θ(τ, k)[1:lmax], [description = "Temperature perturbation multipoles"]
        δ(τ, k), [description = "Overdensity"]
        θ(τ, k), [description = "Velocity divergence"]
        u(τ, k), [description = "Velocity"]
        σ(τ, k), [description = "Shears tress"]
        κ̇(τ), [description = "Optical depth derivative"]
        θb(τ, k), [description = "Baryon velocity divergence"]
        Π(τ, k), [description = "Anisotropic stress perturbation"]
        Π̇(τ, k), [description = "Anisotropic stress perturbation derivative"]
        G0(τ, k), [description = "Polarization component 0"]
        G(τ, k)[1:lmax], [description = "Polarization component"]
    end
    eqs = [
        # Parameter equations
        γ.Ω₀ ~ π^2/15 * (kB*γ.T₀)^4 / (ħ^3*c^5) * 8π*GN / (3*(H100*g.h)^2)

        # Bertschinger & Ma (64) with anₑσₜ -> -κ̇
        D(F0) ~ -k*F[1] + 4*D(g.Φ)
        D(F[1]) ~ k/3*(F0-2*F[2]+4*g.Ψ) - 4//3 * κ̇/k * (θb - θ) # D(θ) ~ -κ̇ (θb-θγ)
        [D(F[l]) ~ k/(2l+1) * (l*F[l-1] - (l+1)*F[l+1]) + κ̇ * (F[l] - δkron(l,2)//10*Π) for l in 2:lmax-1]...
        D(F[lmax]) ~ k*F[lmax-1] - (lmax+1) / τ * F[lmax] + κ̇ * F[lmax] # τ ≈ 1/ℋ
        δ ~ F0
        θ ~ 3*k*F[1]/4
        σ ~ F[2]/2
        u ~ θ / k
        Π ~ F[2] + G0 + G[2]
        Π̇ ~ D(Π)
        Θ0 ~ F0/4
        [Θ[l] ~ F[l]/4 for l in 1:lmax]...
    ]
    ics = [
        F0 ~ -2*g.Ψ # Dodelson (7.89) # TODO: derive automatically
        F[1] ~ 2//3 * k*τ*g.Ψ # Dodelson (7.95)
        F[2] ~ (polarization ? -8//15 : -20//45) * k/κ̇ * F[1] # depends on whether polarization is included
        [F[l] ~ -l//(2*l+1) * k/κ̇ * F[l-1] for l in 3:lmax]...
    ]
    if polarization
        append!(eqs, [
            D(G0) ~ k * (-G[1]) + κ̇ * (G0 - Π/2)
            D(G[1]) ~ k/(2*1+1) * (1*G0 - 2*G[2]) + κ̇ * G[1]
            [D(G[l]) ~ k/(2l+1) * (l*G[l-1] - (l+1)*G[l+1]) + κ̇ * (G[l] - δkron(l,2)//10*Π) for l in 2:lmax-1]...
            D(G[lmax]) ~ k*G[lmax-1] - (lmax+1) / τ * G[lmax] + κ̇ * G[lmax]
        ])
        append!(ics, [
            G0 ~ 5//16 * F[2],
            G[1] ~ -1//16 * k/κ̇ * F[2],
            G[2] ~ 1//16 * F[2],
            [G[l] ~ -l//(2l+1) * k/κ̇ * G[l-1] for l in 3:lmax]...
        ])
    else
        append!(eqs, [collect(G .~ 0)...]) # pin to zero
    end
    description = "Photon radiation"
    return extend(γ, System(eqs, τ, vars, []; initialization_eqs=ics, name, kwargs...); description)
end
