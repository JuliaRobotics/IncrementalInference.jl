module IncrementalInference

info("Multithreaded  convolutions possible, Threads.nthreads()=$(Threads.nthreads()).  See `addFactor!(.;threadmodel=MultiThreaded)`.")

import Base: convert
import HDF5: root
import Distributions: sample
import Base: rand, rand!
import KernelDensityEstimate: kde!

using
  Graphs,
  NLsolve,
  Optim,
  Distributions,
  StatsBase,
  KernelDensityEstimate,
  HDF5,
  JLD,
  ProgressMeter,
  DocStringExtensions,
  Compat

export
  # pass through from Graphs.jl
  # plot,

  # added methods to functions from KernelDensityEstimate
  kde!,
  getPoints,

  # pass through functions commonly used lower down
  Npoints,
  Ndim,
  getBW,

  evalLikelihood,

  # data layer variables
  dlapi,
  localapi,
  showcurrentdlapi,
  # and callback setting function
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

  #functors need
  getSample,
  freshSamples!,

  #Visualization
  writeGraphPdf,
  ls,
  lsf,
  ls2,
  hasOrphans,
  allnums,
  isnestednum,
  sortnestedperm,
  getfnctype,
  drawCopyFG,

  # Tree stuff
  # spyCliqMat,
  evalPotential,
  evalFactor2,
  approxConv,
  approxConvBinary,

  # weiged sampling
  AliasingScalarSampler,
  rand!,
  rand,
  fastnorm,

  # dev
  CommonConvWrapper, # new wrapper (experimental) -- not ready for use

  # is deprecated
  FastGenericRoot,
  FastRootGenericWrapParam,
  GenericWrapParam,

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

  # For 1D example, should be refactored and renamed
  ## TODO will be deprecated
  Odo,
  odoAdd,
  PackedOdo,
  Obsv2,
  PackedObsv2,

  # TODO rename to ball radius
  Ranged,
  PackedRanged,

  # dev exports
  addGraphsVert!,
  makeAddEdge!,

  # define evalPotential functions outside IIF
  registerCallback!,

  # development
  shuffleXAltD,
  reshapeVec2Mat




const VoidUnion{T} = Union{Void, T}


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
include("TreePotentials01.jl")
include("ExplicitDiscreteMarginalizations.jl")
include("ApproxConv.jl")
include("SolveTree01.jl")


include("deprecated.jl")

# function plot(fg::FactorGraph)
#   Graphs.plot(fg.g)
# end

end
