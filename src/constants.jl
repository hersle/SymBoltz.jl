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

# Hydrogen singlet transitions
const λH∞1s   =  91.17534e-9; const fH∞1s  = c/λH∞1s;  const EH∞1s  = h*fH∞1s # ∞ - 1s (read: wavelength Hydrogen ∞ to 1s)
const λH2s1s  = 121.56700e-9; const fH2s1s = c/λH2s1s; const EH2s1s = h*fH2s1s # 2s - 1s
const EH∞2s  = EH∞1s - EH2s1s # E(∞) - E(2s)

# Helium singlet transitions
const λHe∞1s  =  50.42590e-9; const fHe∞1s  = c/λHe∞1s;  const EHe∞1s  = h*fHe∞1s
const λHe2s1s =  60.14045e-9; const fHe2s1s = c/λHe2s1s; const EHe2s1s = h*fHe2s1s
const λHe2p1s =  58.43344e-9; const fHe2p1s = c/λHe2p1s; const EHe2p1s = h*fHe2p1s
const EHe2p2s = EHe2p1s - EHe2s1s
const EHe∞2s  = EHe∞1s - EHe2s1s
const EHe⁺∞1s = 54.4178 * eV

# Helium triplet transitions
const λHet∞2s = 260.0463e-9; const fHet∞2s = c/λHet∞2s; const EHet∞2s = h*fHet∞2s # ∞ - 2³s; ionization of lowest triplet state (4.77 or 4.8 eV) (read: "wavelength Helium triplet ∞ to 2s")
const λHet2p1s = 59.1411e-9; const fHet2p1s = c/λHet2p1s; const EHet2p1s = h*fHet2p1s
const λHet2s1s = 62.5563e-9; const fHet2s1s = c/λHet2s1s; const EHet2s1s = h*fHet2s1s
const EHet2p2s = EHet2p1s - EHet2s1s

δkron(i, j) = (i == j ? 1 : 0) # Kronecker delta

k_dimensionless(k::Number, h) = k
k_dimensionless(k::Quantity, h) = NoUnits(k / (h*H100 / c / u"m"))
k_dimensionless(k::Number, bgsol::ODESolution) = k
k_dimensionless(k::Quantity, bgsol::ODESolution) = k_dimensionless(k, getsym(bgsol, :h)(bgsol))
