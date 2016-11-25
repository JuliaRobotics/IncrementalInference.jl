module IncrementalInference

blas_set_num_threads(2)

using
  Graphs,
  GraphViz,
  Gadfly,
  Colors,
  NLsolve,
  Distributions,
  KernelDensityEstimate,
  TransformUtils,
  CloudGraphs, Neo4j

export
  # actual CloudGraphs integration experimental code
  setCloudDataLayerAPI!,

  # using either dictionary or cloudgraphs
  VariableNodeData,
  PackedVariableNodeData,
  VNDencoder,
  VNDdecoder,
  FNDencode,
  FNDdecode,
  FunctionNodeData,
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

  #Visualization
  investigateMultidimKDE,
  writeGraphPdf,
  ls,
  drawHorDens,
  drawHorBeliefsList,
  # vstackedDensities, # global scope naming issue with msgPlots

  # Tree stuff
  spyCliqMat,
  evalPotential,
  evalFactor2,

  # solve inference
  inferOverTree!,
  inferOverTreeR!,
    #development interface
    upMsgPassingRecursive,

  Singleton,
  Pairwise,
  numericRoot,
  numericRootGenericRandomized,
  GenericMarginal,
  #Robot stuff
  PriorPose2,
  PackedPriorPose2,
  Pose2Pose2,
  PackedPose2Pose2,
  addPose2Pose2,
  Pose2DPoint2DBearingRange,
  Pose2DPoint2DRange,
  Point2DPoint2DRange,
  PriorPoint2D,
  PackedPose2DPoint2DBearingRange,
  solveLandm,
  solvePose2,
  solveSetSeps,
  addPose2Pose2!,
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

  PriorPose3,
  Pose3Pose3,
  projectParticles,

  # dev exports
  addGraphsVert!,
  makeAddEdge!,

  # CloudGraph stuff
  registerGeneralVariableTypes!,
  fullLocalGraphCopy!,
  removeGenericMarginals!,
  setBackendWorkingSet!,
  setDBAllReady!,
  getExVertFromCloud,
  getAllExVertexNeoIDs,
  getPoseExVertexNeoIDs,
  copyAllNodes!,
  copyAllEdges!,

  # DIDSON sonar model
LinearRangeBearingElevation,
project!,
project,
backprojectRandomized!,
residual!,
ominus,
evalPotential,
getSample



include("FactorGraphTypes.jl")
include("CloudGraphIntegration.jl") # Work in progress code
include("DataLayerAPI.jl")
include("FactorGraph01.jl")
include("JunctionTree.jl")
include("GraphConstraintTypes.jl")
include("TreePotentials01.jl")
include("TreePotentials02.jl")
include("TreePotentials03.jl")
include("SensorModels.jl")
include("SolveTree01.jl")
include("SolverVisualization.jl")

end
