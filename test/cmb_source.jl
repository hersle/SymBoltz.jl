using SymBoltz, Plots

M = SymBoltz.ΛCDM(; lmax = 6, h = nothing, Λanalytical = true)
pars = SymBoltz.parameters_Planck18(M)

# 1) Source functionz

lmax = 1000
Δk, Δk_S = 2π/24, 10.0 # TODO: factor in t0
sol1, k1s = SymBoltz.solve_for_cmb(M, pars, lmax; Δk, Δk_S = Δk_S, verbose = true)
sol2, k2s = SymBoltz.solve_for_cmb(M, pars, lmax; Δk, Δk_S = Δk,   verbose = true)
ks = range(extrema(k2s)..., length=50)
ts = exp.(range(-3, -2, length=250))

S1 = SymBoltz.source_temperature(sol1, ks, ts) # sol1(ks, ts, M.S)
S2 = SymBoltz.source_temperature(sol2, ks, ts) # sol2(ks, ts, M.S)
println("max(abs(S2-S1)) = ", maximum(abs.(S2 .- S1)))

p = plot(xlabel = "ln(t)", ylabel = "S")
plot!(p, log.(ts), transpose(S1); linestyle = :solid, color = permutedims(1:length(ks)), alpha=0.5)
plot!(p, log.(ts), transpose(S2); linestyle = :dash,  color = permutedims(1:length(ks)), alpha=0.5)

# 2) LOS integration

ls = 10:10:lmax
ks = range(extrema(sol2.ks)..., length=4*length(sol2.ks))
lnts = range(extrema(log.(sol2[M.t]))..., length=500)
S2 = SymBoltz.source_temperature(sol2, ks, exp.(lnts)) # sol2(ks, exp.(lnts), M.S)
Θl2s = SymBoltz.los_integrate(S2, ls, ks, lnts)

plot(ks, Θl2s[:,1], xlims=(0,100))

# 3) CMB power spectrum

P0s = SymBoltz.spectrum_primordial(ks)
Cl2s = SymBoltz.spectrum_cmb(Θl2s, Θl2s, P0s, ls, ks)
Dl2s = SymBoltz.Dl(Cl2s, ls)

plot(ls, Dl2s)

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
tmin, tmax = extrema(sol[M.t])
lnts = range(log(tmin), log(tmax), step=0.01)
SPs = SymBoltz.source_polarization(sol, ks, exp.(lnts)) # TODO: cut away early/late times? (e.g. S blows up today)
#plot(lnts, Ss'; label=nothing)

ΘPls = SymBoltz.los_integrate(SPs, ls, ks, lnts)
ΘPls .*= transpose(@. √((ls+2)*(ls+1)*(ls+0)*(ls-1)))

#plot(ks, Θls, xlims=(0,100))

P0s = SymBoltz.spectrum_primordial(ks)
Cls = SymBoltz.spectrum_cmb(ΘPls, ΘPls, P0s, ls, ks)
Dls = SymBoltz.Dl(Cls, ls)
plot(log10.(ls), Dls)
