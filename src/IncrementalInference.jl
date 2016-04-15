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
  NBPMessage,
  ExploreTreeType,
  addNode!,
  addFactor!,
  getVal,
  BayesTree,
  emptyFactorGraph,
  getEliminationOrder,
  prepBatchTree!,
  wipeBuildNewTree!,
  whichCliq,
  get2DSamples,
  getAll2D,
  get2DSampleMeans,
  getAll2DMeans,
  getAll2DPoses,
  get2DPoseSamples,
  get2DPoseMeans,
  getKDE,
  getVertKDE,
  get2DPoseMax,
  getAll2DLandmarks,
  get2DLandmSamples,
  get2DLandmMeans,
  get2DLandmMax,
  writeGraphPdf,
  ls,
  # Tree stuff
  spyCliqMat,
  evalPotential,
  evalFactor2,
  inferOverTree!,
  inferOverTreeR!,

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

  #development interface
  upMsgPassingRecursive,

  #Visualization stuff for robots should be moved to RoME
  drawPosesLandms,
  getKDEMax,
  spyCliqMat



include("FactorGraph01.jl")
include("TreePotentials01.jl")
include("TreePotentials02.jl")
include("SolveTree01.jl")
include("SolverVisualization.jl")

end
