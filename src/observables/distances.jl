# TODO: add formula
"""
    distance_luminosity(sol::CosmologySolution, ivs = sol.bg.t, τ0 = sol[sol.prob.M.τ0])

Compute luminosity distances
```math
d_L = \\frac{r}{a} = \\chi \\, \\mathrm{sinc} (\\sqrt{K} (τ₀-τ)),
```
at the independent variable values `ivs` relative to the (present) time `τ0`.
"""
function distance_luminosity(sol::CosmologySolution, ivs = sol.bg.t, τ0 = sol[sol.prob.M.τ0])
    M = sol.prob.M
    χ = sol(M.χ, ivs)
    Ωk0 = have(M, :K) ? sol[M.K.Ω₀] : 0.0
    r = sinc.(√(-Ωk0+0im)*χ/π) .* χ |> real # Julia's sinc(x) = sin(π*x) / (π*x)
    H0 = H100 * sol[M.g.h]
    a = sol(M.g.a, ivs)
    return @. r / a * c / H0 # to meters
end

# TODO: test @inferred
function distance_luminosity_function(M::System, pars_fixed, pars_varying, zs; bgopts = (alg = Tsit5(), reltol = 1e-5, maxiters = 1e3))
    isequal(ModelingToolkit.get_iv(M), M.g.a) || error("Independent variable must be $(M.g.a)")

    pars = merge(pars_fixed, Dict(pars_varying .=> NaN))
    as = @. 1 / (zs + 1)
    prob = CosmologyProblem(M, pars; pt = false, ivspan = (minimum(as), 1.0))
    probgen = parameter_updater(prob, pars_varying; build_initializeprob = Val{false})

    geta = getsym(prob, M.g.a)
    getτ = getsym(prob, M.τ)
    geth = getsym(prob, M.g.h)
    getΩk0 = getsym(prob, M.K.Ω₀)

    return p -> begin
        prob = probgen(p)
        sol = solve(prob; bgopts, saveat = as, save_end = true)
        a = geta(sol)
        τ = getτ(sol)
        h = geth(sol)
        Ωk0 = getΩk0(sol)
        τ0 = τ[end] # time today
        χ = τ0 .- τ
        r = @. real(sinc(√(-Ωk0+0im)*χ/π) * χ) # Julia's sinc(x) = sin(π*x) / (π*x)
        H0 = H100 * h
        return @. r / a * c / H0 # luminosity distance in meters
    end
end

@doc raw"""
    sound_horizon(sol::CosmologySolution)

Cumulatively integrate the sound horizon
```math
    rₛ(τ) = ∫_0^τ dτ cₛ = ∫_0^τ \frac{dτ}{√(3(1+3ρ_b/4ρ_γ))}
```
to the time steps of the solution `sol`.
"""
function sound_horizon(sol::CosmologySolution)
    M = sol.prob.M
    return integrate_cumulative(sol, 1 / √(3(1+3/4*M.b.ρ/M.γ.ρ))) # TODO: nonzero initial value?
end