module SymBoltz

# re-export commonly used ModelingToolkit functions
import ModelingToolkit: unknowns, observed, parameters, equations, initialization_equations, defaults
export unknowns, observed, parameters, equations, initialization_equations, defaults

using ModelingToolkit
using DifferentialEquations
using NonlinearSolve
using DataInterpolations

# TODO: non-linear: higher-order perturbations vs halofit vs N-body?
# TODO: baryons: Recfast -> Recfast++ -> CosmoRec -> HyRec -> HyRec-2: call out, or integrate equations into my code to make use of my background calculation?
# TODO: composable models, generate equations
# TODO: modified gravity: quintessence -> Brans-Dicke -> coupled quintessence; DGP, parametrized framework, EFT of LSS, ...
# TODO: GPU-parallellized EnsembleProblem
# TODO: define components with @mtkmodel?
# TODO: try different AD sensitivity algorithms: https://docs.sciml.ai/SciMLSensitivity/stable/getting_started/
# TODO: define global spacetime structure here?
# TODO: solve BG in reverse, thermo forward, then perturbations forward?
# TODO: add ODESystems, ODEProblems, ... into one single CosmologyModel type

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
