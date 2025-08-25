# TODO: avoid this; make recombination work naturally with Unitful units?
using PhysicalConstants, Unitful, UnitfulAstro, PeriodicTable
const c  = PhysicalConstants.CODATA2018.c_0 / u"m/s"
const h  = PhysicalConstants.CODATA2018.h / u"J*s"
const ħ  = PhysicalConstants.CODATA2018.ħ / u"J*s"
const kB = PhysicalConstants.CODATA2018.k_B / u"J/K"
const GN = PhysicalConstants.CODATA2018.G / u"m^3/kg/s^2"
const α  = PhysicalConstants.CODATA2018.α # fine structure constant, ≈ 1/137

const σ  = PhysicalConstants.CODATA2018.σ / u"W*K^-4*m^-2" # Stefan-Boltzmann constant ("radiation constant"?)
const aR = 4/c*σ # "radiation constant" in https://arxiv.org/pdf/astro-ph/9909275
const σT = PhysicalConstants.CODATA2018.σ_e / u"m^2"
const km  = 1u"km/m"  |> NoUnits
const pc = 1u"pc/m" |> NoUnits
const kpc = 1u"kpc/m" |> NoUnits
const Mpc = 1u"Mpc/m" |> NoUnits
const Gpc = 1u"Gpc/m" |> NoUnits
const H100 = 100 * km/Mpc
const k0 = H100 * Mpc / c # h/Mpc
const eV = 1u"eV/J" |> NoUnits

const me = PhysicalConstants.CODATA2018.m_e / u"kg"
const mH = elements[:H].atomic_mass / u"kg" |> NoUnits
const mHe = elements[:He].atomic_mass / u"kg" |> NoUnits

δkron(i, j) = (i == j ? 1 : 0) # Kronecker delta

function k_dimensionless(k::Number, h)
    if unit(k) == NoUnits
        return k
    else
        H₀ = h * H100 # s⁻¹
        k0 = (H₀ / c) / u"m"
        return NoUnits(k / k0)
    end
end
