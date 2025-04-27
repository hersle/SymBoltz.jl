using SymBoltz, Plots

M = SymBoltz.ΛCDM(; K = nothing, lmax = 6)
pars = SymBoltz.parameters_Planck18(M)

# 1) Source functions

lmax = 1000
ls = 10:10:lmax
k1s, k2s = SymBoltz.cmb_ks(ls[end])
prob = CosmologyProblem(M, pars)
sol1 = solve(prob, k1s; verbose = true, ptopts = (alg = SymBoltz.Rodas5P(), reltol = 1e-15))
sol2 = solve(prob, k2s; verbose = true)
ks = range(extrema(k2s)..., length=10)
τs = sol1[M.τ]

S1 = SymBoltz.source_temperature(sol1, ks, τs) #; sw=false, isw=false, dop=false, pol=true) # sol1(ks, τs, M.S)
S2 = SymBoltz.source_temperature(sol2, ks, τs) #; sw=false, isw=false, dop=false, pol=true) # sol2(ks, τs, M.S)
println("max(abs(S2-S1)) = ", maximum(abs.(S2 .- S1)))

p = plot(xlabel = "tanh(τ)", ylabel = "S", legend = nothing)
plot!(p, tanh.(τs), transpose(S1); linestyle = :solid, color = permutedims(1:length(ks)), alpha=0.5)
plot!(p, tanh.(τs), transpose(S2); linestyle = :dash,  color = permutedims(1:length(ks)), alpha=0.5)

# 2) LOS integration

ks = k2s #range(extrema(k2s)..., length=1*length(k2s))
τs, us, u′s = SymBoltz.los_substitution_range(sol1, (τ->tanh(τ)), (u->atanh(u)), (τ->1/cosh(τ)^2), length=500)
S = SymBoltz.source_temperature(sol1, ks, τs) # sol2(ks, exp.(lnτs), M.S)
Θls = SymBoltz.los_integrate(S, ls, ks, τs, us, u′s)

plot(ks, Θls[:,10], xlims=(0, 100))

# 3) Cl integrand

P0s = SymBoltz.spectrum_primordial(ks, sol1)
#P0s = sol1(ks, sol1[M.τ][begin], M.I.P)
Cls = similar(Θls, length(ls))
ks_with0 = [0.0; ks] # add dummy value with k=0 for integration
il = 1
dCl_dks_with0 = zeros(eltype(Θls), length(ks_with0))
@. dCl_dks_with0[2:end] = 2/π * ks^2 * P0s * Θls[:,il] * Θls[:,il]
#Cls[il] = integrate(ks_with0[2:end], dCl_dks_with0[2:end], integrator) # integrate over k (_with0 adds one additional point at (0,0))

plot(ks_with0, dCl_dks_with0; xlims=(0, 20))


# 4) CMB power spectrum
Cls = SymBoltz.spectrum_cmb(Θls, Θls, P0s, ls, ks)
Dls = SymBoltz.Dl(Cls, ls)

plot(log10.(ls), Dls)

# 4) With varying precision parameters

p = plot()
for Δk in [2π/16, 2π/32]
    Cls = SymBoltz.spectrum_cmb(:TT, M, pars, ls; Δk, verbose=true)
    Dls = SymBoltz.Dl(Cls, ls)
    plot!(p, ls, Dls; label = "Δk = $Δk")
end

p = plot()
for Δk_S in [10.0, 5.0, 2.5]
    Cls = SymBoltz.spectrum_cmb(:TT, M, pars, ls; Δk_S, verbose=true)
    Dls = SymBoltz.Dl(Cls, ls)
    plot!(p, ls, Dls; label = "Δk_S = $Δk_S")
end

p = plot()
for Δlnt in [0.09, 0.06, 0.03]
    Cls = SymBoltz.spectrum_cmb(:TT, M, pars, ls; Δlnt, verbose=true)
    Dls = SymBoltz.Dl(Cls, ls)
    plot!(p, ls, Dls; label = "Δlnt = $Δlnt")
end

p = plot()
for F in [1.0, 1.5, 2.0]
    Cls = SymBoltz.spectrum_cmb(:TT, M, pars, ls; kmax = F * ls[end], verbose=true)
    Dls = SymBoltz.Dl(Cls, ls)
    plot!(p, ls, Dls; label = "kmax/lmax = $F")
end

# E-mode source functions
ls = 10:5:1500
sol, ks = SymBoltz.solve_for_cmb(M, pars, ls[end])
τmin, τmax = extrema(sol[M.τ])
lnτs = range(log(τmin), log(τmax), step=0.01)
SPs = SymBoltz.source_polarization(sol, ks, exp.(lnτs)) # TODO: cut away early/late times? (e.g. S blows up today)
#plot(lnτs, Ss'; label=nothing)

ΘPls = SymBoltz.los_integrate(SPs, ls, ks, lnτs)
ΘPls .*= transpose(@. √((ls+2)*(ls+1)*(ls+0)*(ls-1)))

#plot(ks, Θls, xlims=(0,100))

P0s = SymBoltz.spectrum_primordial(ks)
Cls = SymBoltz.spectrum_cmb(ΘPls, ΘPls, P0s, ls, ks)
Dls = SymBoltz.Dl(Cls, ls)
plot(log10.(ls), Dls)
