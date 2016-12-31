module IncrementalInference

using
  Graphs,
  GraphViz,
  Gadfly,
  Colors,
  NLsolve,
  Distributions,
  KernelDensityEstimate,
  TransformUtils

export
  # using either dictionary or cloudgraphs
  VariableNodeData,
  PackedVariableNodeData,
  VNDencoder,
  VNDdecoder,
  FNDencode,
  FNDdecode,
  FunctionNodeData,
  PackedFunctionNodeData,
  FactorGraph,
  addNode!,
  addFactor!,
  resetData!,
  getVert,
  getData,
  getVarNode,
  getVal,
  setVal!,
  getBWVal,
  setBW!,
  setValKDE!,
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

  #functors need
  getSample,

  #Visualization
  investigateMultidimKDE,
  writeGraphPdf,
  ls,
  drawHorDens,
  drawHorBeliefsList,

  # Tree stuff
  spyCliqMat,
  evalPotential,
  evalFactor2,

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

  # Solving utils
  numericRoot,
  numericRootGenericRandomized,
  numericRootGenericRandomizedFnc,
  numericRootGenericRandomizedFnc!,
  FastGenericRoot,
  FastRootGenericWrapParam,

  # generic marginal used during elimitation game
  GenericMarginal,
  PackedGenericMarginal,

  uppA,
  convert, # for protobuf stuff
  compare,

  # Going to move to RoME.jl in future
  # For 1D example
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

  # CloudGraphs integration callback setting function
  setdatalayerAPI!,

  # development
  shuffleXAltD




include("FactorGraphTypes.jl")
include("DataLayerAPI.jl")
include("FactorGraph01.jl")
include("JunctionTree.jl")
include("GraphConstraintTypes.jl")
include("SolverUtilities.jl")
include("TreePotentials01.jl")
include("SolveTree01.jl")
include("SolverVisualization.jl")

end
