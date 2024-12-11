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
const Mpc = 1u"Mpc/m" |> NoUnits
const Gpc = 1u"Gpc/m" |> NoUnits
const H100 = 100 * km/Mpc
const k0 = H100 * Mpc / c # h/Mpc
const eV = 1u"eV/J" |> NoUnits

me = PhysicalConstants.CODATA2018.m_e / u"kg"
mH = elements[:H].atomic_mass / u"kg" |> NoUnits
mHe = elements[:He].atomic_mass / u"kg" |> NoUnits

# Hydrogen transitions
λ_H_∞_1s   =  91.17534e-9; f_H_∞_1s  = c/λ_H_∞_1s;  E_H_∞_1s  = h*f_H_∞_1s # ∞ - 1s
λ_H_2s_1s  = 121.56700e-9; f_H_2s_1s = c/λ_H_2s_1s; E_H_2s_1s = h*f_H_2s_1s # 2s - 1s
                                                    E_H_∞_2s  = E_H_∞_1s - E_H_2s_1s # E_∞ - E_2s

# Helium singlet transitions
λ_He_∞_1s  =  50.42590e-9; f_He_∞_1s  = c/λ_He_∞_1s;  E_He_∞_1s  = h*f_He_∞_1s
λ_He_2s_1s =  60.14045e-9; f_He_2s_1s = c/λ_He_2s_1s; E_He_2s_1s = h*f_He_2s_1s
λ_He_2p_1s =  58.43344e-9; f_He_2p_1s = c/λ_He_2p_1s; E_He_2p_1s = h*f_He_2p_1s
                                                      E_He_2p_2s = E_He_2p_1s - E_He_2s_1s
                                                      E_He_∞_2s  = E_He_∞_1s - E_He_2s_1s
                                                      E_He⁺_∞_1s = 54.4178 * eV

# Helium triplet transitions # TODO: rename s,t to singlet,triplet?
λ_He_2p_1s_tri = 59.1411e-9; f_He_2p_1s_tri = c/λ_He_2p_1s_tri; E_He_2p_1s_tri = h*f_He_2p_1s_tri
λ_He_2s_1s_tri = 62.5563e-9; f_He_2s_1s_tri = c/λ_He_2s_1s_tri; E_He_2s_1s_tri = h*f_He_2s_1s_tri
                                                                E_He_2p_2s_tri = E_He_2p_1s_tri - E_He_2s_1s_tri

δkron(i, j) = (i == j ? 1 : 0)

function k_dimensionless(k::Number, h)
    if unit(k) == NoUnits
        return k
    else
        H₀ = h * H100 # s⁻¹
        k0 = (H₀ / c) / u"m"
        return NoUnits(k / k0)
    end
end
