module IncrementalInference

# @info "Multithreaded  convolutions possible, Threads.nthreads()=$(Threads.nthreads()).  See `addFactor!(.;threadmodel=MultiThreaded)`."

using Distributed
using Requires
using Reexport

@reexport using Distributions
@reexport using KernelDensityEstimate
@reexport using ApproxManifoldProducts
# @reexport using Graphs
@reexport using LinearAlgebra

using Manifolds

export ℝ, AbstractManifold, Euclidean, Circle

import NLsolve
import NLSolversBase
import Optim

using
  Dates,
  TimeZones,
  DistributedFactorGraphs,
  DelimitedFiles,
  Statistics,
  Random,
  StatsBase,
  JLD2,
  BSON,
  FileIO,
  ProgressMeter,
  DocStringExtensions,
  FunctionalStateMachine,
  JSON2,
  Combinatorics,
  UUIDs,
  TensorCast

using StaticArrays

using ManifoldsBase

# for BayesTree
using MetaGraphs

using Logging

# bringing in BSD 3-clause ccolamd
include("ccolamd.jl")
using SuiteSparse.CHOLMOD: SuiteSparse_long # For CCOLAMD constraints.
using .Ccolamd

# likely overloads or not exported by the upstream packages
import Base: convert, ==
import Distributions: sample
import Random: rand, rand!
import KernelDensityEstimate: getBW
import KernelDensityEstimate: getPoints
import ApproxManifoldProducts: kde!, manikde!
import ApproxManifoldProducts: getBW
import ApproxManifoldProducts: mmd
import ApproxManifoldProducts: isPartial
import DistributedFactorGraphs: addVariable!, addFactor!, ls, lsf, isInitialized
import DistributedFactorGraphs: compare, compareAllSpecial
import DistributedFactorGraphs: rebuildFactorMetadata!
import DistributedFactorGraphs: getDimension, getManifold, getPointType, getPointIdentity
import DistributedFactorGraphs: getPPE, getPPEDict
import DistributedFactorGraphs: getFactorOperationalMemoryType

# will be deprecated in IIF
import DistributedFactorGraphs: isSolvable


# must be moved to their own repos
const KDE = KernelDensityEstimate
const MB = ManifoldsBase
const AMP = ApproxManifoldProducts
const FSM = FunctionalStateMachine
const IIF = IncrementalInference

# Package aliases
export KDE, AMP, DFG, FSM, IIF


# TODO temporary for initial version of on-manifold products
KDE.setForceEvalDirect!(true)

const InstanceType{T} = Union{Type{<:T},T}

# DFG SpecialDefinitions
export AbstractDFG,
  InMemDFGType,
  getSolverParams,
  LightDFG,
  findShortestPathDijkstra, isPathFactorsHomogeneous,
  getSolvedCount, isSolved, setSolvedCount!,
  listSupersolves, listSolveKeys,
  deepcopySolvekeys!, deepcopySupersolve!,
  diagm,
  listDataEntries,
  FolderStore,
  addBlobStore!,
  addData!,
  getData,
  DFGVariable,
  DFGVariableSummary, 
  DFGFactor,
  DFGFactorSummary,
  deleteVariableSolverData!
  # listDataBlobs  # ERROR: LightDFG{} doesn't override 'listDataBlobs'.

# Inference types
export FunctorInferenceType, PackedInferenceType
export AbstractPrior, AbstractRelative, AbstractRelativeRoots, AbstractRelativeMinimize

# not sure if this is necessary
export convert

export *,
  CSMHistory,
  # getTreeCliqsSolverHistories,

  AbstractBayesTree,
  BayesTreeNodeData,
  PackedBayesTreeNodeData,

  # state machine methods
  StateMachine,
  exitStateMachine,
  printGraphSummary,
  printSummary,
  print,
  getGraphFromHistory,
  getCliqSubgraphFromHistory,
  sandboxStateMachineStep,
  # draw and animate state machine
  getStateLabel,
  histStateMachineTransitions,
  histGraphStateMachineTransitions,
  drawStateTransitionStep,
  drawStateMachineHistory,
  animateStateMachineHistoryByTime,
  animateStateMachineHistoryByTimeCompound,
  animateCliqStateMachines,
  makeCsmMovie,
  areSiblingsRemaingNeedDownOnly,

  # general types for softtyping of variable nodes
  BeliefArray,
  InferenceVariable,
  ContinuousScalar,
  SamplableBelief,
  PackedSamplableBelief,
  Prior,
  PackedPrior,
  MsgPrior,
  PackedMsgPrior,
  PartialPrior,
  PackedPartialPrior,

  ls2,
  # lsRear,
  # from DFG
  ls,
  lsf,
  listVariables,
  listFactors,
  exists,
  sortDFG,
  getLabel,
  getVariables,
  getVariableOrder,
  getPPE,
  getPPEDict,
  getVariablePPE,
  isVariable,
  isFactor,
  getFactorType,
  getSofttype,
  getVariableType,
  getLogPath,
  joinLogPath,
  lsfPriors,
  isPrior,
  lsTypes,
  lsfTypes,
  findClosestTimestamp,
  printVariable,
  printFactor,
  getTimestamp,
  deepcopyGraph,
  deepcopyGraph!,
  copyGraph!,
  getSolverData,

  # using either dictionary or cloudgraphs
  FactorMetadata,
  FunctionNodeData,
  PackedFunctionNodeData, # moved to DFG
  normalfromstring,
  categoricalfromstring,
  # extractdistribution,

  SolverParams,
  getSolvable,
  setSolvable!,
  addVariable!,
  deleteVariable!,
  addFactor!,
  deleteFactor!,
  addMsgFactors!,
  deleteMsgFactors!,
  factorCanInitFromOtherVars,
  doautoinit!,
  initManual!,
  initVariableManual!,
  resetInitialValues!,
  resetInitValues!,
  # asyncTreeInferUp!,
  # initInferTreeUp!,
  solveCliqWithStateMachine!,
  resetData!,
  resetTreeCliquesForUpSolve!,
  resetFactorGraphNewTree!,
  setVariableInitialized!,
  setVariableInferDim!,
  resetVariable!,
  getFactor,
  getFactorDim,
  getVariableDim,
  getVariableInferredDim,
  getVariableInferredDimFraction,
  getVariableSolvableDim,
  getFactorSolvableDim,
  getFactorInferFraction,
  isCliqFullDim,
  getVariable,
  getCliqueData,
  setCliqueData!,
  getManifold,  # new Manifolds.jl based operations
  getVal,
  getBW,
  setVal!,
  getNumPts,
  getBWVal,
  setBW!,
  setValKDE!,
  buildCliqSubgraph,

  #
  isPartial,
  isInitialized,
  isTreeSolved,
  isUpInferenceComplete,
  isCliqInitialized,
  isCliqUpSolved,
  areCliqVariablesAllInitialized,
  ensureSolvable!,
  initAll!,
  cycleInitByVarOrder!,
  getOutNeighbors,
  BayesTree, MetaBayesTree,
  TreeBelief,
  LikelihoodMessage,
  initfg,
  buildSubgraph,
  buildCliqSubgraph!,
  transferUpdateSubGraph!,
  getEliminationOrder,
  buildBayesNet!,
  buildTree!,
  buildTreeReset!,
  buildCliquePotentials,

  getCliqDepth,
  getTreeAllFrontalSyms,
  getTreeCliqUpMsgsAll,
  childCliqs,
  getChildren,
  parentCliq,
  getParent,
  getCliqSiblings,
  getNumCliqs,
  getBelief, getKDE,
  CliqStateMachineContainer,

  solveCliqUp!,
  solveCliqDown!,
  fifoFreeze!,

  # temp const types TODO
  TempUpMsgPlotting,

  #functors need
  getSample,
  sampleFactor!,
  sampleFactor,

  #Visualization
  drawGraph,
  drawGraphCliq,
  drawCliqSubgraphUpMocking,
  drawTree,
  drawTreeAsyncLoop,

  # Bayes (Junction) Tree
  evalFactor,
  approxConvBelief,
  approxConv,
  approxConvBinary,

  # more debugging tools
  localProduct,
  treeProductUp,
  approxCliqMarginalUp!,
  dontMarginalizeVariablesAll!,
  unfreezeVariablesAll!,
  resetVariableAllInitializations!,
  isMarginalized,
  setMarginalized!,
  isMultihypo,
  getMultihypoDistribution,
  getHypothesesVectors,

  # weiged sampling
  AliasingScalarSampler,
  rand!,
  rand,
  randToPoints,
  fastnorm,

  # new wrapper (experimental)
  CommonConvWrapper,
  
  getCliqVarInitOrderUp,
  getCliqNumAssocFactorsPerVar,

  # introduced for approximate convolution operations
  setThreadModel!,
  SingleThreaded,
  MultiThreaded,

  # user functions
  predictbelief,
  getCliqMat,
  getCliqAssocMat,
  getCliqMsgMat,
  getCliqFrontalVarIds,
  getFrontals,
  getCliqSeparatorVarIds,
  getCliqAllVarIds,
  getCliqVarIdsAll,
  getCliqAllVarSyms,
  getCliqVarIdsPriors,
  getCliqVarSingletons,
  getCliqAllFactIds,
  getCliqFactorIdsAll,
  getCliqFactors,
  areCliqVariablesAllMarginalized,

  # generic marginal used during elimitation game
  GenericMarginal,
  PackedGenericMarginal,

  # factor graph operating system utils (fgos)
  saveTree,
  loadTree,

  # Temp placeholder for evaluating string types to real types
  saveDFG,
  loadDFG!,  loadDFG,
  rebuildFactorMetadata!,

  getCliqVarSolveOrderUp,
  getFactorsAmongVariablesOnly,
  setfreeze!,

  #internal dev functions for recycling cliques on tree
  attemptTreeSimilarClique,

  # some utils
  compare,
  compareAllSpecial,
  getMeasurements,
  findFactorsBetweenFrom,
  addDownVariableFactors!,
  getDimension,
  getPointType, 
  getPointIdentity,
  setVariableRefence!,
  reshapeVec2Mat,
  accumulateFactorChain



export  buildCliqSubgraph_StateMachine


const NothingUnion{T} = Union{Nothing, T}
const BeliefArray{T} = Union{Array{T,2}, Adjoint{T, Array{T,2}} }


# regular
include("entities/SolverParams.jl")

# JT TODO move to somewhere more fitting? (DF, perhaps not remember its IIF.SolverParams)
const InMemDFGType = DFG.LightDFG{SolverParams}


include("entities/FactorOperationalMemory.jl")


include("Factors/GenericMarginal.jl")
include("entities/OptionalDensities.jl")
include("entities/FactorGradients.jl")

# Special belief types for sampling as a distribution
include("AliasScalarSampling.jl")
include("HeatmapSampler.jl")

include("entities/BeliefTypes.jl")

# factors needed for belief propagation on the tree
include("Factors/MsgPrior.jl")

include("entities/CliqueTypes.jl")
include("entities/JunctionTreeTypes.jl")

include("FactorGraph.jl")
include("SerializingDistributions.jl")
include("SerializationMKD.jl")
include("DispatchPackedConversions.jl")

include("FGOSUtils.jl")
include("CompareUtils.jl")
include("NeedsResolution.jl")

# tree and init related functions
include("SubGraphFunctions.jl")
include("JunctionTree.jl")
include("TreeMessageAccessors.jl")
include("TreeMessageUtils.jl")
include("TreeBasedInitialization.jl")



# included variables of IIF, easy to extend in user's context
include("Variables/DefaultVariables.jl")
include("Variables/Circular.jl")

# included factors, see RoME.jl for more examples
include("Factors/GenericFunctions.jl")
include("Factors/MsgLikelihoods.jl")
include("Factors/Mixture.jl")
include("Factors/DefaultPrior.jl")
include("Factors/LinearRelative.jl")
include("Factors/EuclidDistance.jl")
include("Factors/Circular.jl")
include("Factors/PartialPrior.jl")
include("Factors/PartialPriorPassThrough.jl")
include("DefaultNodeTypes.jl") # older file


# Refactoring in progress
include("services/CalcFactor.jl")
# gradient tools
include("services/FactorGradients.jl")
include("services/CliqueTypes.jl")


# solving graphs
include("SolverUtilities.jl")
include("NumericalCalculations.jl")
include("DeconvUtils.jl")
include("ExplicitDiscreteMarginalizations.jl")
include("InferDimensionUtils.jl")
include("services/EvalFactor.jl")
include("services/ApproxConv.jl")
include("GraphProductOperations.jl")
include("SolveTree.jl")
include("TetherUtils.jl")
include("TreeDebugTools.jl")
include("CliqStateMachineUtils.jl")
include("CSMOccuranceUtils.jl")

#EXPERIMENTAL parametric
include("ParametricCSMFunctions.jl")
include("ParametricUtils.jl")

#X-stroke
include("CliqueStateMachine.jl")

include("CanonicalGraphExamples.jl")

include("AdditionalUtils.jl")
include("SolverAPI.jl")

# Symbolic tree analysis files.
include("AnalysisTools.jl")

# deprecation legacy support
include("Deprecated.jl")



exportimg(pl) = error("Please do `using Gadfly` to allow image export.")
function __init__()
  @require InteractiveUtils = "b77e0a4c-d291-57a0-90e8-8db25a27a240" include("RequireInteractiveUtils.jl")
  @require Gadfly="c91e804a-d5a3-530f-b6f0-dfbca275c004" include("EmbeddedPlottingUtils.jl")
  @require DifferentialEquations="0c46a032-eb83-5123-abaf-570d42b7fbaa" include("ODE/DERelative.jl")
  @require Interpolations="a98d9a8b-a2ab-59e6-89dd-64a1c18fca59" include("HeatmapSampler.jl")

  # combining neural networks natively into the non-Gaussian  factor graph object
  @require Flux="587475ba-b771-5e3f-ad9e-33799f191a9c" begin
    include("Flux/FluxModelsDistribution.jl")
    include("Flux/FluxModelsSerialization.jl") # uses BSON
  end
end


export setSerializationNamespace!, getSerializationModule, getSerializationModules

end
