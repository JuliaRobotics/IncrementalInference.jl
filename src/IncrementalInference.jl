module IncrementalInference

blas_set_num_threads(2)

using
  Graphs,
  GraphViz,
  Gadfly,
  Colors,
  NLsolve,
  Distributions,
  KernelDensityEstimate

export
  DataLayerAPI,
  setDataLayerAPI,
  VariableNodeData,
  PackedVariableNodeData,
  FunctionNodeData,
  FactorGraph,
  addNode!,
  addFactor!,
  getVal,
  setVal!,
  setBW!,
  setValKDE!,
  BayesTree,
  EasyMessage,
  NBPMessage,
  ExploreTreeType,
  emptyFactorGraph,
  getEliminationOrder,
  prepBatchTree!,
  wipeBuildNewTree!,
  whichCliq,
  getKDE,
  getVertKDE,
  # get2DSamples,
  # getAll2D,
  # get2DSampleMeans,
  # getAll2DMeans,
  # getAll2DPoses,
  # get2DPoseSamples,
  # get2DPoseMeans,
  # get2DPoseMax,
  # getAll2DLandmarks,
  # get2DLandmSamples,
  # get2DLandmMeans,
  # get2DLandmMax,

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

  #Robot stuff
  PriorPose2,
  PackedPriorPose2,
  Pose2Pose2,
  PackedPose2Pose2,
  addPose2Pose2,
  Pose2DPoint2DBearingRange,
  PackedPose2DPoint2DBearingRange,
  solveLandm,
  solvePose2,
  solveSetSeps,
  addPose2Pose2!,
  uppA,
  convert, # for magic protobuf stuff
  compare,

  # For 1D example
  Odo,
  odoAdd,
  PackedOdo,
  Obsv2,
  PackedObsv2,
  Ranged,

  # should improve abstraction
  R,
  se2vee,
  SE2,
  wrapRad


include("DataLayerAPI.jl")
include("FactorGraph01.jl")
include("GraphConstraintTypes.jl")
include("TreePotentials01.jl")
include("TreePotentials02.jl")
include("SolveTree01.jl")
include("SolverVisualization.jl")

end
