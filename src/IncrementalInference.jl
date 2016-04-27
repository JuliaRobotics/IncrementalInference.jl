module IncrementalInference

using
  Graphs,
  GraphViz,
  Gadfly,
  Colors,
  NLsolve,
  Distributions,
  KernelDensityEstimate

export
  FactorGraph,
  addNode!,
  addFactor!,
  getVal,
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
  Pose2Pose2,
  addPose2Pose2,
  Pose2DPoint2DBearingRange,
  solveLandm,
  solvePose2,
  solveSetSeps,
  addPose2Pose2!,
  uppA,

  # should improve abstraction
  R,
  se2vee,
  SE2,
  wrapRad



include("FactorGraph01.jl")
include("TreePotentials01.jl")
include("TreePotentials02.jl")
include("SolveTree01.jl")
include("SolverVisualization.jl")

end
