module SymBoltz

using Reexport
import ModelingToolkit: parameters, unknowns # explicitly imported because they are extended
@reexport using ModelingToolkit # re-export commonly used ModelingToolkit functions
using OrdinaryDiffEq
using NonlinearSolve
using OhMyThreads
using Base.Threads
using Setfield
using StaticArrays

# TODO: generate gravity equations
# TODO: modified gravity: coupled quintessence; DGP, parametrized framework, EFT of LSS, ...
# TODO: define components with @mtkmodel?
# TODO: try different AD sensitivity algorithms: https://docs.sciml.ai/SciMLSensitivity/stable/getting_started/
# TODO: connector systems for Compton scattering / recombination etc.

# TODO: descriptions
@independent_variables τ [description = "Conformal time"] # conformal time in units of 1/H₀
D = Differential(τ)
k = only(GlobalScope.(@parameters k [description = "Perturbation mode wavenumber"])) # perturbation wavenumber

include("utils.jl")
include("constants.jl")
include("models/metric.jl")
include("models/gravity.jl")
include("models/generic_species.jl")
include("models/baryons.jl")
include("models/photons.jl")
include("models/neutrinos.jl")
include("models/dark_matter.jl")
include("models/dark_energy.jl")
include("models/curvature.jl")
include("models/inflation.jl")
include("models/cosmologies.jl")
include("solve.jl")
include("observables/distances.jl")
include("observables/fourier.jl")
include("observables/angular.jl")
include("parameters.jl")
include("plot.jl")

export RMΛ, ΛCDM, w0waCDM, QCDM, GRΛCDM, BDΛCDM
export CosmologyProblem, CosmologySolution
export solve, solvebg, solvept, remake, issuccess, parameter_updater
export parameters_Planck18
export spectrum_primordial, spectrum_matter, spectrum_matter_nonlinear, spectrum_cmb, correlation_function, variance_matter, stddev_matter, los_integrate, source_grid, sound_horizon, SphericalBesselCache
export express_derivatives

end
