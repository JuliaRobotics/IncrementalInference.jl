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
  check_procs,

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
  convert, # for protobuf stuff
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


end
