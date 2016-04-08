module IncrementalInference

export
  FactorGraph,
  addNode!,
  addFactor!,
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
  spyCliqMat


include("FactorGraph.jl")
include("TreePotentials01.jl")
include("TreePotentials02.jl")
include("SolveTree01.jl")
include("SolverVisualization.jl")

end
