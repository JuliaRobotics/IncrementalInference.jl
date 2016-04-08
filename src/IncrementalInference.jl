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

  # should improve abstraction
  R



include("FactorGraph01.jl")
include("TreePotentials01.jl")
include("TreePotentials02.jl")
include("SolveTree01.jl")
include("SolverVisualization.jl")

end
