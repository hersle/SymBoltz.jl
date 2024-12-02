# Sources and spectra

## Primordial power spectra

```@docs
SymBoltz.P0
```

#### Example

```@example
# TODO: create InflationModel or something
using SymBoltz, Unitful, UnitfulAstro, Plots
ks = 10 .^ range(-5, +1, length=100) / u"Mpc"
Ps = SymBoltz.P0(ks; As = 2.1e-9)
plot(log10.(ks*u"Mpc"), log10.(Ps/u"Mpc^3"); xlabel = "log10(k/Mpc⁻¹)", ylabel = "log10(P/Mpc³)")
```

## Matter power spectra

```@docs
SymBoltz.power_spectrum
```

#### Example

```@example
using SymBoltz, Unitful, UnitfulAstro, Plots
M = SymBoltz.ΛCDM()
pars = SymBoltz.parameters_Planck18(M)
ks = 10 .^ range(-5, +1, length=100) / u"Mpc"
Ps = SymBoltz.power_spectrum(M, pars, ks)
plot(log10.(ks*u"Mpc"), log10.(Ps/u"Mpc^3"); xlabel = "log10(k/Mpc⁻¹)", ylabel = "log10(P/Mpc³)")
```

## CMB power spectra

```@docs
SymBoltz.Cl
```

#### Example

```@example
# TODO: more generic, source functions, ... # hide
using SymBoltz, Unitful, UnitfulAstro, Plots
M = SymBoltz.ΛCDM()
pars = SymBoltz.parameters_Planck18(M)
ls = 10:10:1500
Cls = SymBoltz.Cl(M, pars, ls)
Dls = SymBoltz.Dl(Cls, ls)
plot(ls, Dls; xlabel = "l", ylabel = "Dₗ")
```
