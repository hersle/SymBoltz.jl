module SymBoltz

# re-export commonly used ModelingToolkit functions
import ModelingToolkit: unknowns, observed, parameters, equations, initialization_equations, defaults, hierarchy
export unknowns, observed, parameters, equations, initialization_equations, defaults, hierarchy

using ModelingToolkit
using OrdinaryDiffEq
using NonlinearSolve
using OhMyThreads
using Setfield

# TODO: generate gravity equations
# TODO: modified gravity: coupled quintessence; DGP, parametrized framework, EFT of LSS, ...
# TODO: define components with @mtkmodel?
# TODO: try different AD sensitivity algorithms: https://docs.sciml.ai/SciMLSensitivity/stable/getting_started/
# TODO: connector systems for Compton scattering / recombination etc.

using ModelingToolkit: t_nounits as t, D_nounits as D # t is conformal time in units of 1/H₀
k = only(GlobalScope.(@parameters k)) # perturbation wavenumber
ϵ = only(GlobalScope.(@parameters ϵ = 1)) # perturbative expansion parameter

include("utils.jl")
include("constants.jl")
include("components/metric.jl")
include("components/gravity.jl")
include("components/species.jl")
include("components/thermodynamics.jl")
include("components/inflation.jl")
include("models.jl")
include("solve.jl")
include("spectra.jl")
include("parameters.jl")
include("plot.jl")

export RMΛ, ΛCDM, w0waCDM, QCDM, GRΛCDM, BDΛCDM
export CosmologyProblem, CosmologySolution
export solve, shoot, remake
export parameters_Planck18
export spectrum_primordial, spectrum_matter, spectrum_matter_nonlinear, spectrum_cmb, correlation_function, variance_matter, stddev_matter

end
