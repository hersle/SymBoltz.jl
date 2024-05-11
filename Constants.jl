# TODO: avoid this; make recombination work naturally with Unitful units?
using PhysicalConstants, Unitful, UnitfulAstro
const c  = PhysicalConstants.CODATA2018.c_0 / u"m/s"
const h  = PhysicalConstants.CODATA2018.h / u"J*s"
const ħ  = PhysicalConstants.CODATA2018.ħ / u"J*s"
const kB = PhysicalConstants.CODATA2018.k_B / u"J/K"
const G  = PhysicalConstants.CODATA2018.G / u"m^3/kg/s^2"
const α  = PhysicalConstants.CODATA2018.α # fine structure constant, ≈ 1/137
const me = PhysicalConstants.CODATA2018.m_e / u"kg"
const mp = PhysicalConstants.CODATA2018.m_p / u"kg"
const σT = PhysicalConstants.CODATA2018.σ_e / u"m^2"
const km  = 1u"km/m"  |> NoUnits
const Mpc = 1u"Mpc/m" |> NoUnits
const EHion = 13.59844u"eV/J" |> NoUnits
const EHe1ion = 24.58738u"eV/J" |> NoUnits
const EHe2ion = 54.41776u"eV/J" |> NoUnits
k0 = 1 / 2997.92458 # h/Mpc

δkron(i, j) = (i == j ? 1 : 0)