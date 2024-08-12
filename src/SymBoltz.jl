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

export background_ΛCDM, thermodynamics_ΛCDM, perturbations_ΛCDM
export P0, P, Cl, Dl

end