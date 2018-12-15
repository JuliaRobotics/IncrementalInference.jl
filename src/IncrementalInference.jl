module IncrementalInference

@info "Multithreaded  convolutions possible, Threads.nthreads()=$(Threads.nthreads()).  See `addFactor!(.;threadmodel=MultiThreaded)`."

using Distributed
using Reexport

@reexport using Distributions
@reexport using KernelDensityEstimate
@reexport using Graphs
@reexport using LinearAlgebra

using
  Statistics,
  Random,
  NLsolve,
  StatsBase,
  JLD2,
  FileIO,
  ProgressMeter,
  DocStringExtensions,
  Optim # might be deprecated in favor for only NLsolve dependency



const KDE = KernelDensityEstimate

import Base: convert
# import HDF5: root
import Distributions: sample
import Random: rand, rand!
import KernelDensityEstimate: kde!

export
  KDE,
  dlapi,  # data layer variables
  localapi,
  showcurrentdlapi,
  setdatalayerAPI!,
  DataLayerAPI,

  # general types for softtyping of variable nodes
  InferenceVariable,
  ContinuousScalar,
  ContinuousMultivariate,
  SamplableBelief,
  Prior,
  PackedPrior,
  LinearConditional,
  PackedLinearConditional,
  MixturePrior,
  PackedMixturePrior,
  MixtureLinearConditional,
  PackedMixtureLinearConditional,

  # using either dictionary or cloudgraphs
  VariableNodeData,
  PackedVariableNodeData,
  FactorMetadata,
  encodePackedType,
  FunctionNodeData,
  PackedFunctionNodeData,
  encodePackedType,
  decodePackedType,
  normalfromstring,
  categoricalfromstring,
  extractdistribution,

  FactorGraph,
  addNode!,
  addFactor!,
  doautoinit!,
  manualinit!,
  resetData!,
  getVert,
  getData,
  setData!,
  getVarNode,
  getVal,
  setVal!,
  getNumPts,
  getBWVal,
  setBW!,
  setValKDE!,
  isInitialized,
  ensureAllInitialized!,
  updateFullVert!,
  getOutNeighbors,
  BayesTree,
  EasyMessage,
  NBPMessage,
  ExploreTreeType,
  emptyFactorGraph,
  subgraphShortestPath,
  subgraphFromVerts,
  getEliminationOrder,
  buildBayesNet!,
  emptyBayesTree,
  buildTree!,
  prepBatchTree!,
  wipeBuildNewTree!,
  whichCliq,
  childCliqs,
  parentCliq,
  getKDE,
  getVertKDE,
  initializeNode!,
  batchSolve!,
  fifoFreeze!,
  getCurrentWorkspaceFactors,
  getCurrentWorkspaceVariables,

  #functors need
  getSample,
  freshSamples!,

  #Visualization
  writeGraphPdf,
  ls,
  lsf,
  ls2,
  lsRear,
  hasOrphans,
  printgraphmax,
  allnums,
  isnestednum,
  sortnestedperm,
  getfnctype,
  drawCopyFG,

  # Bayes (Junction) Tree
  evalPotential,
  evalFactor2,
  approxConv,
  approxConvBinary,

  # more debugging tools
  localProduct,  # TODO remove from Caesar docs make.jl import list v0.4.4+
  treeProductUp,

  # weiged sampling
  AliasingScalarSampler,
  rand!,
  rand,
  fastnorm,

  # new wrapper (experimental)
  CommonConvWrapper,


  # solve inference
  inferOverTree!,
  inferOverTreeR!,

  #development interface
  upMsgPassingRecursive,

  # Inference types
  InferenceType,
  PackedInferenceType,
  Singleton,
  Pairwise,
  # introduced for approximate convolution operations
  setThreadModel!,
  SingleThreaded,
  MultiThreaded,

  # functor abstracts
  FunctorInferenceType,
  FunctorPairwise,
  FunctorPairwiseMinimize,
  FunctorSingleton,
  # FunctorPartialSingleton,
  FunctorPairwiseNH,
  FunctorSingletonNH,

  # Solving utils
  findRelatedFromPotential,
  shuffleXAltD!,
  numericRoot,
  numericRootGenericRandomized,
  numericRootGenericRandomizedFnc,
  numericRootGenericRandomizedFnc!,

  # user functions
  proposalbeliefs,
  predictbelief,
  getCliqMat,
  getCliqMsgMat,

  # generic marginal used during elimitation game
  GenericMarginal,
  PackedGenericMarginal,

  uppA,
  convert,
  compare,
  extractdistribution,

  # factor graph operating system utils (fgos)
  convert2packedfunctionnode,
  encodefg,
  decodefg,
  savejld,
  loadjld,
  landmarks,

  # Temp placeholder for evaluating string types to real types
  _evalType,

  setUpMsg!,
  upMsg,
  setDwnMsg!,
  dwnMsg,

  compareField,
  compareFields,
  compareAll,

  # For 1D example,

  # TODO rename to ball radius
  Ranged,
  PackedRanged,

  # development
  # dev exports
  addGraphsVert!,
  makeAddEdge!,
  shuffleXAltD,
  reshapeVec2Mat # TODO deprecate




const NothingUnion{T} = Union{Nothing, T}

include("ccolamd.jl")

include("FactorGraphTypes.jl")
include("AliasScalarSampling.jl")
include("DefaultNodeTypes.jl")
include("DataLayerAPI.jl")
include("FactorGraph01.jl")
include("SerializingDistributions.jl")
include("DispatchPackedConversions.jl")
include("FGOSUtils.jl")

include("JunctionTree.jl")
include("GraphConstraintTypes.jl")
include("SolverUtilities.jl")
include("ExplicitDiscreteMarginalizations.jl")
include("ApproxConv.jl")
include("SolveTree01.jl")

# Hack for RoME module.
global serializationnamespace = Dict{String, Module}()

"""
    $(SIGNATURES)

De-serialization of IncrementalInference objects require discovery of foreign types.

Example:

Template to tunnel types from a user module:
```julia
# or more generic solution -- will always try Main if available
IIF.setSerializationNamespace!("Main" => Main)

# or a specific package such as RoME
using RoME
IIF.setSerializationNamespace!("RoME" => RoME)
```
"""
function setSerializationNamespace!(keyval::Pair{String, Module})
  serializationnamespace[keyval[1]] = keyval[2]
end


# global _romeModule = nothing
function setRoMEModule(romeModule)::Nothing
    error("Deprecated setRoMEModule, use IIF.setSerializationNamespace!(\"Main\" => Main) function instead")
    # global _romeModule
    _romeModule = romeModule
    return nothing
end

function getSerializationModule(mod::String="Main")::Module
  global serializationnamespace
  return serializationnamespace[mod]
end

# Old code that might be used again
# function getType(typestring::AS) where {AS <: AbstractString}
#  # eval(Meta.parse(typestring))()
#  # getfield(Main, Symbol(typestring))
#  getfield(@__MODULE__, Symbol(typestring))
# end
# function getRoMEModule()
#     return _romeModule
# end
function getRoMEModule()::Module
  @warn "Deprecated getRoMEModule, use getSerializationModule(...) instead"
  return getSerializationModule("RoME")
end

export setSerializationNamespace!
export setRoMEModule, getRoMEModule, getSerializationModule

end
