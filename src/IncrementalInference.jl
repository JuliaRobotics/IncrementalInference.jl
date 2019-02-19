module IncrementalInference

@info "Multithreaded  convolutions possible, Threads.nthreads()=$(Threads.nthreads()).  See `addFactor!(.;threadmodel=MultiThreaded)`."

using Distributed
using Reexport

@reexport using Distributions
@reexport using KernelDensityEstimate
@reexport using ApproxManifoldProducts
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
const AMP = ApproxManifoldProducts

import Base: convert
# import HDF5: root
import Distributions: sample
import Random: rand, rand!
import KernelDensityEstimate: getBW
import ApproxManifoldProducts: kde!

# TODO temporary for initial version of on-manifold products
KDE.setForceEvalDirect!(true)

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
  addVariable!,
  addFactor!,
  doautoinit!,
  manualinit!,
  resetData!,
  getVert,
  getData,
  setData!,
  getSofttype,
  getVarNode,
  getVal,
  getBW,
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
  FactorGraph,
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
  getChildren,
  parentCliq,
  getParent,
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
  getUpMsgs,
  setDwnMsg!,
  dwnMsg,
  getDwnMsg,

  # some utils
  compareField,
  compareFields,
  compareAll,
  getIdx,
  showFactor,

  # For 1D example,

  # TODO rename to L2 distance
  Ranged,
  PackedRanged,

  # development
  # dev exports
  addGraphsVert!,
  makeAddEdge!,
  shuffleXAltD,
  reshapeVec2Mat # TODO deprecate




const NothingUnion{T} = Union{Nothing, T}

# non-free, but not currently use!
include("ccolamd.jl")

# regular
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

include("Deprecated.jl")

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
  global serializationnamespace
  serializationnamespace[keyval[1]] = keyval[2]
end

function getSerializationModule(mod::String="Main")::Union{Module, Nothing}
  global serializationnamespace
  if haskey(serializationnamespace, mod)
    return serializationnamespace[mod]
  end
  return nothing
end

function getSerializationModules()::Dict{String, Module}
  global serializationnamespace
  return serializationnamespace
end

# Old code that might be used again
# function getType(typestring::AS) where {AS <: AbstractString}
#  # eval(Meta.parse(typestring))()
#  # getfield(Main, Symbol(typestring))
#  getfield(@__MODULE__, Symbol(typestring))
# end

export setSerializationNamespace!, getSerializationModule, getSerializationModules

end
