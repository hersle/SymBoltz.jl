module Symboltz

using ModelingToolkit
using SymbolicIndexingInterface: variable_index
using DifferentialEquations
import DifferentialEquations: solve
using DataInterpolations

# TODO: make simpler Cosmology interface
# TODO: non-linear: higher-order perturbations vs halofit vs N-body?
# TODO: baryons: Recfast -> Recfast++ -> CosmoRec -> HyRec -> HyRec-2: call out, or integrate equations into my code to make use of my background calculation?
# TODO: composable models, generate equations
# TODO: modified gravity: quintessence -> Brans-Dicke -> coupled quintessence; DGP, parametrized framework, EFT of LSS, ...
# TODO: GPU-parallellized EnsembleProblem
# TODO: define components with @mtkmodel?
# TODO: try different AD sensitivity algorithms: https://docs.sciml.ai/SciMLSensitivity/stable/getting_started/
# TODO: define global spacetime structure here?

@variables η # conformal time in units of 1/H0
Dη = Differential(η) # d/dη

@parameters k # perturbation wavenumber
k = GlobalScope(k)

include("Constants.jl")
include("Background.jl")
include("Thermodynamics.jl")
include("Perturbations.jl")
include("Spectra.jl")

export BackgroundSystem, ThermodynamicsSystem, PerturbationsSystem
export background_ΛCDM, thermodynamics_ΛCDM, perturbations_ΛCDM
export solve
export P0, P, Cl, Dl

end