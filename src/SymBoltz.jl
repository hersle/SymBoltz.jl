module SymBoltz

using ModelingToolkit
using DifferentialEquations
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
@parameters k # perturbation wavenumber
k = GlobalScope(k)

include("Utils.jl")
include("Cosmology.jl")
include("Constants.jl")
include("Components.jl")
include("Thermodynamics.jl")
include("Spectra.jl")
include("Plot.jl")

export ΛCDM
export CosmologyModel, CosmologySolution
export solve
export power_spectrum, Cl, Dl

# re-export commonly used ModelingToolkit functions
export unknowns, observed, parameters
export equations, initialization_equations, defaults

end