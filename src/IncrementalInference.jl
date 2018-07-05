module IncrementalInference

import Base: convert
import HDF5: root
# import KernelDensityEstimate: root
import Distributions: sample
# import KernelDensityEstimate: sample
import KernelDensityEstimate: kde!
# import Graphs: plot

using
  Graphs,
  NLsolve,
  Optim,
  Distributions,
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

  #functors need
  getSample,

  #Visualization
  writeGraphPdf,
  ls,
  lsf,
  ls2,
  getfnctype,
  drawCopyFG,

  # Tree stuff
  # spyCliqMat,
  evalPotential,
  evalFactor2,
  approxConv,

  # dev
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

  # functor abstracts
  FunctorInferenceType,
  FunctorPairwise,
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
  FastGenericRoot,
  FastRootGenericWrapParam,

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
  Odo,
  odoAdd,
  PackedOdo,
  Obsv2,
  PackedObsv2,
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

  # deprecated
  # VNDencoder,
  # VNDdecoder,
  # FNDencode,
  # FNDdecode




const VoidUnion{T} = Union{Void, T}


include("FactorGraphTypes.jl")
include("DataLayerAPI.jl")
include("FactorGraph01.jl")
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
