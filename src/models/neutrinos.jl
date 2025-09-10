"""
    massless_neutrinos(g; lmax = 6, name = :ν, kwargs...)

Create a particle species for massless neutrinos in the spacetime with metric `g`.
"""
function massless_neutrinos(g; lmax = 6, name = :ν, kwargs...)
    description = "Massless neutrinos"
    ν = radiation(g; adiabatic = true, name, description, kwargs...) |> background |> complete

    vars = @variables begin
        F0(τ, k), [description = "Distribution function monopole"]
        F(τ, k)[1:lmax], [description = "Distribution function multipoles"]
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
        D(F[lmax]) ~ k*F[lmax-1] - (lmax+1) / τ * F[lmax]
        δ ~ F0
        θ ~ 3*k*F[1]/4
        σ ~ F[2]/2
    ]
    ics = [
        δ ~ -2 * g.Ψ # adiabatic: δᵢ/(1+wᵢ) == δⱼ/(1+wⱼ) (https://cmb.wintherscoming.no/theory_initial.php#adiabatic)
        θ ~ 1//2 * (k^2*τ) * g.Ψ
        σ ~ 1//15 * (k*τ)^2 * g.Ψ
        [F[l] ~ +l//(2*l+1) * k*τ * F[l-1] for l in 3:lmax]...
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
            D(ψ[i,lmax]) ~ k/(2*lmax+1) * x[i]/E[i] * (lmax*ψ[i,lmax-1] - (lmax+1) * ((2*lmax+1) * E[i]/x[i] * ψ[i,lmax] / (k*τ) - ψ[i,lmax-1])) # explicitly inserted ψ[lmax+1] to avoid array allocations in newer MTK (see example in https://github.com/SciML/ModelingToolkit.jl/issues/3708)
        ])
        append!(ics, [
            ψ0[i] ~ -1//4 * (-2*g.Ψ) * dlnf₀_dlnx(x[i])
            ψ[i,1] ~ -1//3 * E[i]/x[i] * (1/2*k*τ*g.Ψ) * dlnf₀_dlnx(x[i])
            ψ[i,2] ~ -1//2 * (1//15*(k*τ)^2*g.Ψ) * dlnf₀_dlnx(x[i])
            [ψ[i,l] ~ 0 for l in 3:lmax] # TODO: full ICs
        ])
    end

    description = "Massive neutrino"
    pars = [m, m_eV, T₀, #=x, W,=# Ω₀, y₀, T₀, Iρ₀] #  ModelingToolkit.scalarize(E₀)] # need every E₀ index
    return System(eqs, τ, vars, pars; initialization_eqs=ics, #=defaults=defs,=# name, description, kwargs...)
end
