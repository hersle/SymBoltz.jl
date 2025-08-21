"""
    species_constant_eos(g, _w, ẇ = 0, _σ = 0; analytical = true, θinteract = false, adiabatic = false, name = :s, kwargs...)

Create a symbolic component for a particle species with equation of state `w ~ P/ρ` in the spacetime with the metric `g`.
"""
function species_constant_eos(g, _w, ẇ = 0, _σ = 0; analytical = true, θinteract = false, adiabatic = false, name = :s, kwargs...)
    @assert ẇ == 0 && _σ == 0 # TODO: relax (need to include in ICs)
    if analytical
        pars = @parameters begin
            Ω₀, [description = "Reduced background density today"]
        end
    else
        pars = []
    end
    vars = @variables begin
        w(τ), [description = "Equation of state"]
        ρ(τ), [description = "Background density"]
        P(τ), [description = "Background pressure"]
        Ω(τ), [description = "Reduced background density"]
        cₛ²(τ), [description = "Speed of sound squared"]
        δ(τ, k), [description = "Overdensity"]
        θ(τ, k), [description = "Velocity divergence"]
        θinteraction(τ, k), [description = "Velocity divergence interaction"]
        σ(τ, k), [description = "Shear stress"]
        u(τ, k), [description = "Velocity"]
        u̇(τ, k), [description = "Velocity derivative"]
    end
    eqs = [
        w ~ _w
        P ~ w * ρ
        analytical ? (ρ ~ 3/(8*Num(π))*Ω₀ * g.a^(-3*(1+w))) : (D(ρ) ~ -3 * g.ℰ * (ρ + P)) # alternative derivative: D(ρ) ~ -3 * g.ℰ * (ρ + P)
        Ω ~ 8*Num(π)/3 * ρ

        D(δ) ~ -(1+w)*(θ-3*D(g.Φ)) - 3*g.ℰ*(cₛ²-w)*δ # Bertschinger & Ma (30) with Φ -> -Φ; or Baumann (4.4.173) with Φ -> -Φ
        D(θ) ~ -g.ℰ*(1-3*w)*θ - ẇ/(1+w)*θ + cₛ²/(1+w)*k^2*δ - k^2*σ + k^2*g.Ψ + θinteraction # Bertschinger & Ma (30) with θ = kv # TODO: θinteraction mult by ℋ?
        u ~ θ / k
        u̇ ~ D(u)
        σ ~ _σ
    ]
    adiabatic && push!(eqs, cₛ² ~ w)
    ics = [
        δ ~ -3//2 * (1+w) * g.Ψ # adiabatic: δᵢ/(1+wᵢ) == δⱼ/(1+wⱼ) (https://cmb.wintherscoming.no/theory_initial.php#adiabatic) # TODO: match CLASS with higher-order (for photons)? https://github.com/lesgourg/class_public/blob/22b49c0af22458a1d8fdf0dd85b5f0840202551b/source/perturbations.c#L5631-L5632
        θ ~ 1//2 * (k^2/g.ℰ) * g.Ψ # τ ≈ 1/ℰ # TODO: include σ ≠ 0 # solve u′ + ℋ(1-3w)u = w/(1+w)*kδ + kΨ with Ψ=const, IC for δ, Φ=-Ψ, ℋ=H₀√(Ωᵣ₀)/a after converting ′ -> d/da by gathering terms with u′ and u in one derivative using the trick to multiply by exp(X(a)) such that X′(a) will "match" the terms in front of u
    ]
    !θinteract && push!(eqs, (θinteraction ~ 0))
    return System(eqs, τ, vars, [pars; k]; initialization_eqs=ics, name, kwargs...)
end

"""
    matter(g; name = :m, kwargs...)

Create a particle species for matter (with equation of state `w ~ 0`) in the spacetime with metric `g`.
"""
function matter(g; name = :m, kwargs...)
    description = "Matter"
    return species_constant_eos(g, 0; name, description, kwargs...)
end

"""
    radiation(g; name = :r, kwargs...)

Create a particle species for radiation (with equation of state `w ~ 1/3`) in the spacetime with metric `g`.
"""
function radiation(g; name = :r, kwargs...)
    r = species_constant_eos(g, 1//3; name, kwargs...) |> complete
    pars = @parameters begin
        T₀, [description = "Temperature today (in K)"]
    end
    vars = @variables begin
        T(τ), [description = "Temperature"] # TODO: define in constant_eos? https://physics.stackexchange.com/questions/650508/whats-the-relation-between-temperature-and-scale-factor-for-arbitrary-eos-1
    end
    eqs = [T ~ T₀ / g.a]
    description = "Radiation"
    return extend(r, System(eqs, τ, vars, pars; name); description)
end

"""
    cosmological_constant(g; name = :Λ, kwargs...)

Create a particle species for the cosmological constant (with equation of state `w ~ -1`) in the spacetime with metric `g`.
"""
function cosmological_constant(g; name = :Λ, analytical = true, kwargs...)
    description = "Cosmological constant"
    Λ = species_constant_eos(g, -1; name, analytical, adiabatic = true, description, kwargs...) |> background |> complete # discard ill-defined perturbations
    vars = @variables begin
        δ(τ, k), [description = "Overdensity"]
        θ(τ, k), [description = "Velocity divergence"]
        σ(τ, k), [description = "Shear stress"]
    end
    eqs = [δ ~ 0, θ ~ 0, σ ~ 0]
    description = "Cosmological constant"
    return extend(Λ, System(eqs, τ, vars, []; name); description) # manually set perturbations to zero
end

"""
    w0wa(g; kwargs...)

Create a particle species for the w₀-wₐ dark energy (CPL) parametrization in the spacetime with metric `g`.
"""
function w0wa(g; name = :X, analytical = false, kwargs...)
    pars = @parameters begin
        w0, [description = "Equation of state today"]
        wa, [description = "Equation of state evolution"]
        cₛ², [description = "Rest-frame speed of sound squared"]
    end
    vars = @variables begin
        ρ(τ), [description = "Background density"]
        P(τ), [description = "Background pressure"]
        w(τ), [description = "Equation of state"]
        ẇ(τ), [description = "Equation of state derivative"]
        cₐ²(τ), [description = "Adiabatic speed of sound squared"]
        δ(τ, k), [description = "Overdensity"]
        θ(τ, k), [description = "Velocity divergence"]
        σ(τ, k), [description = "Shear stress"]
    end
    # TODO: generate equations with a generic species_eos function
    eqs = [
        w ~ w0 + wa * (1 - g.a)
        ẇ ~ D(w)
        P ~ w * ρ
    ]
    if analytical
        append!(pars, @parameters Ω₀)
        push!(eqs, ρ ~ 3/(8*Num(π))*Ω₀ * abs(g.a)^(-3 * (1 + w0 + wa)) * exp(-3 * wa * (1 - g.a))) # energy density # TODO: get rid of abs
    else
        push!(eqs, D(ρ) ~ -3 * g.ℰ * ρ * (1 + w))
    end
    append!(eqs, [
        # Following https://arxiv.org/pdf/1002.1311 section II
        cₐ² ~ w - ẇ/(3*g.ℰ*(1+w))
        D(δ) ~ -(1+w)*(θ-3*D(g.Φ)) - 3*g.ℰ*(cₛ²-w)*δ - 9*(g.ℰ/k)^2*(1+w)*(cₛ²-cₐ²)*θ
        D(θ) ~ -g.ℰ*(1-3*cₛ²)*θ + cₛ²/(1+w)*k^2*δ - k^2*σ + k^2*g.Ψ
        σ ~ 0
    ])
    ics = [
        δ ~ -3//2 * (1+w) * g.Ψ # adiabatic ICs, see e.g. https://arxiv.org/abs/1004.5509 eq. (3.17)
        θ ~ 1//2 * (k^2/g.ℰ) * g.Ψ # τ ≈ 1/ℰ; adiabatic ICs, see e.g. https://arxiv.org/abs/1004.5509 eq. (3.18)
    ]
    description = "w₀wₐ (CPL) dark energy"
    return System(eqs, τ, vars, pars; initialization_eqs=ics, name, description, kwargs...)
end

"""
    photons(g; polarization = true, lmax = 6, name = :γ, kwargs...)

Create a particle species for photons in the spacetime with metric `g`.
"""
function photons(g; polarization = true, lmax = 6, name = :γ, kwargs...)
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
        D(F[lmax]) ~ k*F[lmax-1] - (lmax+1) * g.ℰ * F[lmax] + κ̇ * F[lmax] # τ ≈ 1/ℰ
        δ ~ F0
        θ ~ 3*k*F[1]/4
        σ ~ F[2]/2
        Π ~ F[2] + G0 + G[2]
        Π̇ ~ D(Π)
        Θ0 ~ F0/4
        [Θ[l] ~ F[l]/4 for l in 1:lmax]...
    ]
    ics = [
        F0 ~ -2*g.Ψ # Dodelson (7.89) # TODO: derive automatically
        F[1] ~ 2//3 * k/g.ℰ*g.Ψ # Dodelson (7.95)
        F[2] ~ (polarization ? -8//15 : -20//45) * k/κ̇ * F[1] # depends on whether polarization is included
        [F[l] ~ -l//(2*l+1) * k/κ̇ * F[l-1] for l in 3:lmax]...
    ]
    if polarization
        append!(eqs, [
            D(G0) ~ k * (-G[1]) + κ̇ * (G0 - Π/2)
            D(G[1]) ~ k/(2*1+1) * (1*G0 - 2*G[2]) + κ̇ * G[1]
            [D(G[l]) ~ k/(2l+1) * (l*G[l-1] - (l+1)*G[l+1]) + κ̇ * (G[l] - δkron(l,2)//10*Π) for l in 2:lmax-1]...
            D(G[lmax]) ~ k*G[lmax-1] - (lmax+1) * g.ℰ * G[lmax] + κ̇ * G[lmax]
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

"""
    massless_neutrinos(g; lmax = 6, name = :ν, kwargs...)

Create a particle species for massless neutrinos in the spacetime with metric `g`.
"""
function massless_neutrinos(g; lmax = 6, name = :ν, kwargs...)
    description = "Massless neutrinos"
    ν = radiation(g; adiabatic = true, name, description, kwargs...) |> background |> complete

    vars = @variables begin
        F0(τ, k), [description = "Distribution function monopole"]
        F(τ, k)[1:lmax+1], [description = "Distribution function multipoles"]
        δ(τ, k), [description = "Overdensity"]
        θ(τ, k), [description = "Velocity divergence"]
        σ(τ, k), [description = "Shear stress"]
    end
    pars = @parameters begin
        Neff, [description = "Effective number of neutrino species"] # TODO: massless vs. massive?
    end
    eqs = [
        D(F0) ~ -k*F[1] + 4*D(g.Φ)
        D(F[1]) ~ k/3*(F0-2*F[2]+4*g.Ψ)
        [D(F[l]) ~ k/(2*l+1) * (l*F[l-1] - (l+1)*F[l+1]) for l in 2:lmax-1]...
        D(F[lmax]) ~ k*F[lmax-1] - (lmax+1) * g.ℰ * F[lmax] # τ ≈ 1/ℰ
        δ ~ F0
        θ ~ 3*k*F[1]/4
        σ ~ F[2]/2
    ]
    ics = [
        δ ~ -2 * g.Ψ # adiabatic: δᵢ/(1+wᵢ) == δⱼ/(1+wⱼ) (https://cmb.wintherscoming.no/theory_initial.php#adiabatic)
        θ ~ 1//2 * (k^2/g.ℰ) * g.Ψ
        σ ~ 1//15 * (k/g.ℰ)^2 * g.Ψ
        [F[l] ~ +l//(2*l+1) * k/g.ℰ * F[l-1] for l in 3:lmax]...
    ]
    description = "Massless neutrinos"
    return extend(ν, System(eqs, τ, vars, pars; initialization_eqs=ics, name, kwargs...); description)
end

# TODO: use vector equations and simplify loops
"""
    massive_neutrinos(g; nx = 5, lmax = 4, name = :h, kwargs...)

Create a particle species for massive neutrinos in the spacetime with metric `g`.
"""
function massive_neutrinos(g; nx = 5, lmax = 4, name = :h, kwargs...)
    pars = @parameters begin
        m, [description = "Neutrino mass (in kg)"]
        m_eV, [description = "Neutrino mass (in eV/c^2)"] # TODO: only one m?
        T₀, [description = "Temperature today (in K)"]
        #x[1:nx], [description = "Dimensionless momentum bins"] # not working with MTKv10 # TODO: reintroduce?
        #W[1:nx], [description = "Gaussian momentum quadrature weights"] # not working with MTKv10 # TODO: reintroduce?
        Ω₀, [description = "Reduced background density today"]
        y₀, [description = "Temperature-reduced mass today"]
        Iρ₀, [description = "Density integral today"]
    end
    vars = @variables begin
        ρ(τ), [description = "Background density"]
        P(τ), [description = "Background pressure"]
        Ω(τ), [description = "Reduced background density"]
        T(τ), [description = "Temperature"]
        y(τ), [description = "Temperature-deuced mass"]
        w(τ), [description = "Equation of state"]
        cₛ²(τ, k), [description = "Speed of sound squared"]
        δ(τ, k), [description = "Overdensity"]
        σ(τ, k), [description = "Shear stress"]
        θ(τ, k), [description = "Velocity divergence"]
        u(τ, k), [description = "Velocity"]
        E(τ)[1:nx], [description = "Dimensionless energies"]
        ψ0(τ, k)[1:nx], [description = "Distribution function monopole"]
        ψ(τ, k)[1:nx,1:lmax], [description = "Distribution function multipoles"]
        In(τ), [description = "Number density integral"]
        Iρ(τ), [description = "Density integral"]
        IP(τ), [description = "Pressure integral"]
        Iδρ(τ, k), [description = "Overdensity integral"]
    end

    f₀(x) = 1 / (exp(x) + 1) # not exp(E); distribution function is "frozen in"; see e.g. Dodelson exercise 3.9
    dlnf₀_dlnx(x) = -x / (1 + exp(-x))
    ∫dx_x²_f₀(f) = sum(collect(f .* W)) # a function that approximates the weighted integral ∫dx*x^2*f(x)*f₀(x)

    # compute numerical reduced momenta x = q*c / (kB*T) and Gaussian quadrature weights, and map the symbolic parameters to them
    x, W = gauss(x -> x^2 * f₀(x), nx, 0.0, 1e3)
    #= # not working with MTKv10 # TODO: reintroduce _x and _W?
    defs = [
        collect(x .=> _x);
        collect(W .=> _W);
    ]
    =#

    eqs = [
        # parameter equations: compute Ω₀ parameter by duplicating time-dependent equations today # TODO: avoid
        m ~ m_eV * SymBoltz.eV/SymBoltz.c^2
        y₀ ~ m*c^2 / (kB*T₀)
        Iρ₀ ~ ∫dx_x²_f₀(@. √(x^2 + y₀^2)) # circumvent defining E₀[1:nx] because vector parameter dependencies doesn't work properly with setsym/remake
        Ω₀ ~ 8π/3 * 2/(2*π^2) * (kB*T₀)^4 / (ħ*c)^3 * Iρ₀ / ((H100*g.h*c)^2/GN)

        T ~ T₀ / g.a
        y ~ y₀ * g.a
        In ~ ∫dx_x²_f₀(1)
        Iρ ~ ∫dx_x²_f₀(E)
        IP ~ ∫dx_x²_f₀(x.^2 ./ E)
        ρ ~ 2/(2*π^2) * (kB*T)^4 / (ħ*c)^3 * Iρ / ((H100*g.h*c)^2/GN) # compute g/(2π²ħ³) * ∫dp p² √((pc)² + (mc²)²) / (exp(pc/(kT)) + 1) with dimensionless x = pc/(kT) and degeneracy factor g = 2
        P ~ 2/(6*π^2) * (kB*T)^4 / (ħ*c)^3 * IP / ((H100*g.h*c)^2/GN) # compute g/(6π²ħ³) * ∫dp p⁴ / √((pc)² + (mc²)²) / (exp(pc/(kT)) + 1) with dimensionless x = pc/(kT) and degeneracy factor g = 2
        w ~ P / ρ
        Ω ~ 8*Num(π)/3 * ρ

        Iδρ ~ ∫dx_x²_f₀(E .* ψ0)
        δ ~ Iδρ / Iρ
        u ~ ∫dx_x²_f₀(x .* ψ[:,1]) / (Iρ + IP/3)
        θ ~ u * k
        σ ~ (2//3) * ∫dx_x²_f₀(x.^2 ./ E .* ψ[:,2]) / (Iρ + IP/3)
        cₛ² ~ ∫dx_x²_f₀(x.^2 ./ E .* ψ0) / Iδρ # TODO: numerator ψ[:,0] or ψ[:,2]?
    ]
    ics = []
    for i in 1:nx
        push!(eqs, E[i] ~ √(x[i]^2 + y^2))
        append!(eqs, [
            D(ψ0[i]) ~ -k * x[i]/E[i] * ψ[i,1] - D(g.Φ) * dlnf₀_dlnx(x[i])
            D(ψ[i,1]) ~ k/3 * x[i]/E[i] * (ψ0[i] - 2*ψ[i,2]) - k/3 * E[i]/x[i] * g.Ψ * dlnf₀_dlnx(x[i])
            [D(ψ[i,l]) ~ k/(2*l+1) * x[i]/E[i] * (l*ψ[i,l-1] - (l+1) * ψ[i,l+1]) for l in 2:lmax-1]...
            D(ψ[i,lmax]) ~ k/(2*lmax+1) * x[i]/E[i] * (lmax*ψ[i,lmax-1] - (lmax+1) * ((2*lmax+1) * E[i]/x[i] * ψ[i,lmax] * g.ℰ/k - ψ[i,lmax-1])) # explicitly inserted ψ[lmax+1] to avoid array allocations in newer MTK (see example in https://github.com/SciML/ModelingToolkit.jl/issues/3708)
        ])
        append!(ics, [
            ψ0[i] ~ -1//4 * (-2*g.Ψ) * dlnf₀_dlnx(x[i])
            ψ[i,1] ~ -1//3 * E[i]/x[i] * (1/2*k/g.ℰ*g.Ψ) * dlnf₀_dlnx(x[i])
            ψ[i,2] ~ -1//2 * (1//15*(k/g.ℰ)^2*g.Ψ) * dlnf₀_dlnx(x[i])
            [ψ[i,l] ~ 0 for l in 3:lmax] # TODO: full ICs
        ])
    end

    description = "Massive neutrino"
    pars = [m, m_eV, T₀, #=x, W,=# Ω₀, y₀, T₀, Iρ₀] #  ModelingToolkit.scalarize(E₀)] # need every E₀ index
    return System(eqs, τ, vars, pars; initialization_eqs=ics, #=defaults=defs,=# name, description, kwargs...)
end

"""
    cold_dark_matter(g; name = :c, kwargs...)

Create a particle species for cold dark matter in the spacetime with metric `g`.
"""
function cold_dark_matter(g; name = :c, kwargs...)
    description = "Cold dark matter"
    return matter(g; adiabatic = true, name, description, kwargs...)
end

"""
    baryons(g; recombination=true, name = :b, kwargs...)

Create a particle species for baryons in the spacetime with metric `g`.
"""
function baryons(g; recombination = true, reionization = true, name = :b, kwargs...)
    description = "Baryons"
    b = matter(g; adiabatic = false, θinteract=true, name, description, kwargs...) |> complete
    if recombination # TODO: simply dont add recombination system when recombination = false
        @named rec = thermodynamics_recombination_recfast(g; reionization)
        eqs = Equation[]
    else
        vars = @variables begin
            κ(τ), [description = "Optical depth"]
            κ̇(τ), [description = "Optical depth derivative"]
            ρb(τ), [description = "Baryon background density"]
            Tγ(τ), [description = "Photon temperature"]
            cₛ²(τ), [description = "Speed of sound squared"]
        end
        eqs = [κ ~ 0, κ̇ ~ 0, cₛ² ~ 0]
        @named rec = System(eqs, τ, vars, [])
    end
    push!(eqs, b.cₛ² ~ rec.cₛ²)
    description = "Baryonic matter"
    b = extend(b, System(eqs, τ, [], []; name); description)
    b = compose(b, rec)
    return b
end

"""
    quintessence(g, v; name = :ϕ, kwargs...)

Create a species for a quintessence scalar field with potential `v` in the spacetime with metric `g`.
"""
function quintessence(g, v, v′, v′′; name = :Q, kwargs...)
    @variables begin
        ϕ(τ), [description = "Background scalar field"]
        ρ(τ), [description = "Effective background density"]
        P(τ), [description = "Effective background pressure"]
        w(τ), [description = "Equation of state"]
        δ(τ, k), [description = "Overdensity"]
        σ(τ, k), [description = "Shear stress"]
        V(τ), [description = "Potential of scalar field"]
        V′(τ), [description = "Potential derivative wrt. scalar field"]
        V′′(τ), [description = "Potential 2nd derivative wrt. scalar field"]
        K(τ), [description = "Effective kinetic energy"]
        m²(τ), [description = "Effective mass"]
        ϵs(τ), [description = "1st slow roll parameter"]
        ηs(τ), [description = "2nd slow roll parameter"]
        cₛ²(τ), [description = "Speed of sound squared"]
    end
    eqs = [
        V ~ v
        V′ ~ v′
        V′′ ~ v′′
        K ~ (D(ϕ)/g.a)^2 / 2 # ϕ̇²/2 = (ϕ′/a)²/2
        D(D(ϕ)) ~ -2 * g.ℰ * D(ϕ) - g.a^2 * V′ # with cosmic time: ϕ̈ + 3*E*ϕ̇ + V′ = 0
        ρ ~ K + V
        P ~ K - V
        w ~ P / ρ
        m² ~ V′′
        ϵs ~ (V′/V)^2 / (16*Num(π))
        ηs ~ (V′′/V) / (8*Num(π))

        δ ~ 0
        σ ~ 0
        cₛ² ~ 0
    ]
    description = "Quintessence dark energy"
    return System(eqs, τ; name, description, kwargs...)
end
function quintessence(g, v; name = :Q, kwargs...)
    @variables begin
        ϕ(τ), [description = "Background scalar field"]
    end
    ∂_∂ϕ = Differential(ϕ)
    v′ = ∂_∂ϕ(v(ϕ)) |> expand_derivatives |> simplify
    v′′ = ∂_∂ϕ(∂_∂ϕ(v(ϕ))) |> expand_derivatives |> simplify
    return quintessence(g, v(ϕ), v′, v′′; name, kwargs...)
end

# TODO: incorporate in metric/gravity; don't define as effective species
"""
    curvature(g; name = :K, kwargs...)

Create a species that effectively accounts for curvature in the spacetime with metric `g`.
"""
function curvature(g; name = :K, kwargs...)
    description = "Curvature"
    vars = @variables begin
        ρ(τ), [description = "Effective background density"]
        P(τ), [description = "Effective background pressure"]
        w(τ), [description = "Effective equation of state"]
        δ(τ, k),
        θ(τ, k),
        cₛ²(τ),
        σ(τ, k)
    end
    pars = @parameters begin
        Ω₀, [description = "Effective reduced background density today"] # dimless K is physical K*c²/H₀²
    end
    eqs = [
        w ~ -1//3
        ρ ~ 3/(8*Num(π))*Ω₀ / g.a^2
        P ~ w*ρ
        δ ~ 0
        θ ~ 0
        cₛ² ~ 0
        σ ~ 0
    ]
    return System(eqs, τ, vars, pars; name, description, kwargs...)
end
