module IncrementalInference

# @info "Multithreaded  convolutions possible, Threads.nthreads()=$(Threads.nthreads()).  See `addFactor!(.;threadmodel=MultiThreaded)`."

using Distributed
using Requires
using Reexport

@reexport using Distributions
@reexport using KernelDensityEstimate
@reexport using ApproxManifoldProducts
# @reexport using Graphs
@reexport using LinearAlgebra

using Manifolds

export ℝ, AbstractManifold
export ProductRepr
# common groups -- preferred defaults at this time.
export TranslationGroup, RealCircleGroup
# common non-groups -- TODO still teething problems to sort out in IIF v0.25-v0.26.
export Euclidean, Circle

import NLsolve
import NLSolversBase
import Optim

using
  Dates,
  TimeZones,
  DistributedFactorGraphs,
  DelimitedFiles,
  Statistics,
  Random,
  StatsBase,
  BSON,
  FileIO,
  ProgressMeter,
  DocStringExtensions,
  FunctionalStateMachine,
  JSON2,
  Combinatorics,
  UUIDs,
  TensorCast

using StaticArrays

using ManifoldsBase

# for BayesTree
using MetaGraphs

using Logging

# bringing in BSD 3-clause ccolamd
include("ccolamd.jl")
using SuiteSparse.CHOLMOD: SuiteSparse_long # For CCOLAMD constraints.
using .Ccolamd

# likely overloads or not exported by the upstream packages
import Base: convert, ==
import Distributions: sample
import Random: rand, rand!
import KernelDensityEstimate: getBW
import KernelDensityEstimate: getPoints
import ApproxManifoldProducts: kde!, manikde!
import ApproxManifoldProducts: getBW
import ApproxManifoldProducts: mmd
import ApproxManifoldProducts: isPartial
import DistributedFactorGraphs: reconstFactorData
import DistributedFactorGraphs: addVariable!, addFactor!, ls, lsf, isInitialized
import DistributedFactorGraphs: compare, compareAllSpecial
import DistributedFactorGraphs: rebuildFactorMetadata!
import DistributedFactorGraphs: getDimension, getManifold, getPointType, getPointIdentity
import DistributedFactorGraphs: getPPE, getPPEDict
import DistributedFactorGraphs: getFactorOperationalMemoryType
import DistributedFactorGraphs: getPoint, getCoordinates
import DistributedFactorGraphs: getVariableType


# will be deprecated in IIF
import DistributedFactorGraphs: isSolvable


# must be moved to their own repos
const KDE = KernelDensityEstimate
const MB = ManifoldsBase
const AMP = ApproxManifoldProducts
const FSM = FunctionalStateMachine
const IIF = IncrementalInference


const InstanceType{T} = Union{Type{<:T},T}
const NothingUnion{T} = Union{Nothing, T}
const BeliefArray{T} = Union{Array{T,2}, Adjoint{T, Array{T,2}} } # TBD deprecate?

## =============================
# API Exports

# Package aliases
export KDE, AMP, DFG, FSM, IIF

# TODO temporary for initial version of on-manifold products
KDE.setForceEvalDirect!(true)

include("ExportAPI.jl")


## =============================
# Source code

# FIXME, move up to DFG
# abstract type AbstractManifoldMinimize <: AbstractRelative end

# regular
include("entities/SolverParams.jl")

# needs SolverParams
const InMemDFGType = DFG.LightDFG{SolverParams}

include("entities/FactorOperationalMemory.jl")

include("Factors/GenericMarginal.jl")
# Special belief types for sampling as a distribution
include("entities/AliasScalarSampling.jl")
include("entities/OptionalDensities.jl")
include("entities/BeliefTypes.jl")
include("entities/FactorGradients.jl")


# Statistics helpers on manifolds
include("VariableStatistics.jl")

# factors needed for belief propagation on the tree
include("Factors/MsgPrior.jl")

include("entities/CliqueTypes.jl")
include("entities/JunctionTreeTypes.jl")

include("services/GraphInit.jl")
include("FactorGraph.jl")
include("services/BayesNet.jl")

# Serialization helpers
include("Serialization/entities/SerializingDistributions.jl")
include("Serialization/services/SerializingDistributions.jl")
include("Serialization/services/SerializationMKD.jl")
include("Serialization/services/DispatchPackedConversions.jl")

include("FGOSUtils.jl")
include("CompareUtils.jl")
include("NeedsResolution.jl")

# tree and init related functions
include("SubGraphFunctions.jl")
include("JunctionTree.jl")
include("TreeMessageAccessors.jl")
include("TreeMessageUtils.jl")
include("TreeBasedInitialization.jl")


# included variables of IIF, easy to extend in user's context
include("Variables/DefaultVariables.jl")
include("Variables/Circular.jl")

# included factors, see RoME.jl for more examples
include("Factors/GenericFunctions.jl")
include("Factors/Mixture.jl")
include("Factors/DefaultPrior.jl")
include("Factors/LinearRelative.jl")
include("Factors/EuclidDistance.jl")
include("Factors/Circular.jl")
include("Factors/PartialPrior.jl")
include("Factors/PartialPriorPassThrough.jl")
include("DefaultNodeTypes.jl") # older file

# Refactoring in progress
include("services/CalcFactor.jl")
# gradient tools
include("services/FactorGradients.jl")
include("services/CliqueTypes.jl")

# solving graphs
include("SolverUtilities.jl")
include("NumericalCalculations.jl")
include("DeconvUtils.jl")
include("ExplicitDiscreteMarginalizations.jl")
# include("InferDimensionUtils.jl")
include("services/EvalFactor.jl")
include("services/ApproxConv.jl")

include("ConsolidateParametricRelatives.jl") # FIXME CONSOLIDATE

include("GraphProductOperations.jl")
include("SolveTree.jl")
include("TetherUtils.jl")
include("TreeDebugTools.jl")
include("CliqStateMachineUtils.jl")

#EXPERIMENTAL parametric
include("ParametricCSMFunctions.jl")
include("ParametricUtils.jl")

#X-stroke
include("CliqueStateMachine.jl")

include("CanonicalGraphExamples.jl")

include("AdditionalUtils.jl")
include("SolverAPI.jl")

# Symbolic tree analysis files.
include("AnalysisTools.jl")

include("ManifoldSampling.jl")

# deprecation legacy support
include("Deprecated.jl")



exportimg(pl) = error("Please do `using Gadfly` to allow image export.")
function __init__()
  @require InteractiveUtils = "b77e0a4c-d291-57a0-90e8-8db25a27a240" include("RequireInteractiveUtils.jl")
  @require Gadfly="c91e804a-d5a3-530f-b6f0-dfbca275c004" include("EmbeddedPlottingUtils.jl")
  @require DifferentialEquations="0c46a032-eb83-5123-abaf-570d42b7fbaa" include("ODE/DERelative.jl")
  @require Interpolations="a98d9a8b-a2ab-59e6-89dd-64a1c18fca59" include("services/HeatmapSampler.jl")

  # combining neural networks natively into the non-Gaussian  factor graph object
  @require Flux="587475ba-b771-5e3f-ad9e-33799f191a9c" begin
    include("Flux/FluxModelsDistribution.jl")
    include("Serializatoin/services/FluxModelsSerialization.jl") # uses BSON
  end
end


export setSerializationNamespace!, getSerializationModule, getSerializationModules

end
