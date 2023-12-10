module IncrementalInference

# @info "Multithreaded  convolutions possible, Threads.nthreads()=$(Threads.nthreads()).  See `addFactor!(.;threadmodel=MultiThreaded)`."

using Distributed
using Reexport

@reexport using Distributions
@reexport using KernelDensityEstimate
@reexport using ApproxManifoldProducts
# @reexport using Graphs
@reexport using LinearAlgebra

using Manifolds
using RecursiveArrayTools: ArrayPartition
export ArrayPartition
using ManifoldDiff
using FiniteDifferences

using OrderedCollections: OrderedDict

import Optim

using Dates,
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
  JSON3,
  Combinatorics,
  UUIDs,
  TensorCast

using StructTypes

using StaticArrays

using ManifoldsBase
using ManifoldsBase: TypeParameter
# for BayesTree
using MetaGraphs
using Logging
using PrecompileTools

# JL 1.10 transition to IncrInfrApproxMinDegreeExt instead
# # bringing in BSD 3-clause ccolamd
# include("services/ccolamd.jl")
# using SuiteSparse.CHOLMOD: SuiteSparse_long # For CCOLAMD constraints.
# using .Ccolamd

# likely overloads or not exported by the upstream packages
import Base: convert, ==, getproperty
import Distributions: sample
import Random: rand, rand!
import KernelDensityEstimate: getBW
import KernelDensityEstimate: getPoints
import ApproxManifoldProducts: kde!, manikde!
import ApproxManifoldProducts: getBW
import ApproxManifoldProducts: mmd
import ApproxManifoldProducts: isPartial
import ApproxManifoldProducts: _update!
import DistributedFactorGraphs: reconstFactorData
import DistributedFactorGraphs: addVariable!, addFactor!, ls, lsf, isInitialized
import DistributedFactorGraphs: compare, compareAllSpecial
import DistributedFactorGraphs: rebuildFactorMetadata!
import DistributedFactorGraphs: getDimension, getManifold, getPointType, getPointIdentity
import DistributedFactorGraphs: getPPE, getPPEDict
import DistributedFactorGraphs: getFactorOperationalMemoryType
import DistributedFactorGraphs: getPoint, getCoordinates
import DistributedFactorGraphs: getVariableType
import DistributedFactorGraphs: AbstractPointParametricEst, loadDFG
import DistributedFactorGraphs: getFactorType
import DistributedFactorGraphs: solveGraph!, solveGraphParametric!

# will be deprecated in IIF
import DistributedFactorGraphs: isSolvable

# must be moved to their own repos
const KDE = KernelDensityEstimate
const MB = ManifoldsBase
const AMP = ApproxManifoldProducts
const FSM = FunctionalStateMachine
const IIF = IncrementalInference

const InstanceType{T} = Union{Type{<:T}, <:T}
const NothingUnion{T} = Union{Nothing, <:T}
const BeliefArray{T} = Union{<:AbstractMatrix{<:T}, <:Adjoint{<:T, AbstractMatrix{<:T}}} # TBD deprecate?

## =============================
# API Exports

# Package aliases
# FIXME, remove this and let the user do either import or const definitions
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

include("entities/HypoRecipe.jl")
include("entities/CalcFactor.jl")
include("entities/FactorOperationalMemory.jl")

include("Factors/GenericMarginal.jl")
# Special belief types for sampling as a distribution
include("entities/AliasScalarSampling.jl")
include("entities/ExtDensities.jl") # used in BeliefTypes.jl::SamplableBeliefs
include("entities/ExtFactors.jl")
include("entities/BeliefTypes.jl")

include("services/HypoRecipe.jl")

#
include("manifolds/services/ManifoldsExtentions.jl")
include("manifolds/services/ManifoldSampling.jl")

include("entities/FactorGradients.jl")

# Statistics helpers on manifolds
include("services/VariableStatistics.jl")

# factors needed for belief propagation on the tree
include("Factors/MsgPrior.jl")
include("Factors/MetaPrior.jl")

include("entities/CliqueTypes.jl")
include("entities/JunctionTreeTypes.jl")

include("services/JunctionTree.jl")
include("services/GraphInit.jl")
include("services/FactorGraph.jl")
include("services/BayesNet.jl")

# Serialization helpers
include("Serialization/entities/SerializingDistributions.jl")
include("Serialization/entities/AdditionalDensities.jl")
include("Serialization/services/SerializingDistributions.jl")
include("Serialization/services/SerializationMKD.jl")
include("Serialization/services/DispatchPackedConversions.jl")

include("services/FGOSUtils.jl")
include("services/CompareUtils.jl")

include("NeedsResolution.jl")

# tree and init related functions
include("services/SubGraphFunctions.jl")
include("services/JunctionTreeUtils.jl")
include("services/TreeMessageAccessors.jl")
include("services/TreeMessageUtils.jl")
include("services/TreeBasedInitialization.jl")

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

# older file
include("services/DefaultNodeTypes.jl")

# Refactoring in progress
include("services/CalcFactor.jl")
# gradient tools
include("services/FactorGradients.jl")
include("services/CliqueTypes.jl")

# solving graphs
include("services/SolverUtilities.jl")
include("services/NumericalCalculations.jl")
include("services/DeconvUtils.jl")
include("services/ExplicitDiscreteMarginalizations.jl")
# include("InferDimensionUtils.jl")
include("services/EvalFactor.jl")
include("services/ApproxConv.jl")


include("services/GraphProductOperations.jl")
include("services/SolveTree.jl")
include("services/TetherUtils.jl")
include("services/TreeDebugTools.jl")
include("CliqueStateMachine/services/CliqStateMachineUtils.jl")

# FIXME CONSOLIDATE
include("parametric/services/ConsolidateParametricRelatives.jl")
include("parametric/entities/MaxMixture.jl")
#EXPERIMENTAL parametric
include("parametric/services/ParametricCSMFunctions.jl")
include("parametric/services/ParametricUtils.jl")
include("parametric/services/ParametricOptim.jl")
include("parametric/services/ParametricManopt.jl")

#X-stroke
include("CliqueStateMachine/services/CliqueStateMachine.jl")

include("services/CanonicalGraphExamples.jl")

include("services/AdditionalUtils.jl")
include("services/SolverAPI.jl")

# Symbolic tree analysis files.
include("services/AnalysisTools.jl")

# extension densities on weakdeps
include("Serialization/entities/SerializingOptionalDensities.jl")
include("Serialization/services/SerializingOptionalDensities.jl")

include("../ext/WeakDepsPrototypes.jl")

# deprecation legacy support
include("Deprecated.jl")



@compile_workload begin
  # In here put "toy workloads" that exercise the code you want to precompile
  fg = generateGraph_Kaess()
  initAll!(fg)
  solveGraph!(fg)
  initParametricFrom!(fg, :default)
  solveGraphParametric!(fg)
end

export setSerializationNamespace!, getSerializationModule, getSerializationModules

end
