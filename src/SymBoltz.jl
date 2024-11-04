module SymBoltz

# re-export commonly used ModelingToolkit functions
import ModelingToolkit: unknowns, observed, parameters, equations, initialization_equations, defaults, hierarchy
export unknowns, observed, parameters, equations, initialization_equations, defaults, hierarchy

using ModelingToolkit
using DifferentialEquations
using NonlinearSolve
using DataInterpolations

# TODO: non-linear: higher-order perturbations / halofit / NBodySystem?
# TODO: CosmoRec / HyRec recombination
# TODO: generate gravity equations
# TODO: modified gravity: coupled quintessence; DGP, parametrized framework, EFT of LSS, ...
# TODO: GPU-parallellized EnsembleProblem
# TODO: define components with @mtkmodel?
# TODO: try different AD sensitivity algorithms: https://docs.sciml.ai/SciMLSensitivity/stable/getting_started/
# TODO: connector systems for Compton scattering / recombination etc.

using ModelingToolkit: t_nounits as t, D_nounits as D # t is conformal time in units of 1/H0
k = only(GlobalScope.(@parameters k)) # perturbation wavenumber
ϵ = only(GlobalScope.(@parameters ϵ)) # perturbative expansion parameter

include("utils.jl")
include("constants.jl")
include("components/metric.jl")
include("components/gravity.jl")
include("components/species.jl")
include("components/thermodynamics.jl")
include("models.jl")
include("solve.jl")
include("spectra.jl")
include("parameters.jl")
include("plot.jl")

export RMΛ, ΛCDM, w0waCDM, QCDM, GRΛCDM, BDΛCDM
export CosmologyModel, CosmologySolution
export solve, shoot
export parameters_Planck18
export power_spectrum, Cl, Dl

end
