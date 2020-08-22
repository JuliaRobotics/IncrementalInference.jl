module IncrementalInference

@info "Multithreaded  convolutions possible, Threads.nthreads()=$(Threads.nthreads()).  See `addFactor!(.;threadmodel=MultiThreaded)`."

using Distributed
using Requires
using Reexport

@reexport using Distributions
@reexport using KernelDensityEstimate
@reexport using ApproxManifoldProducts
@reexport using Graphs
@reexport using LinearAlgebra

using
  Dates,
  TimeZones,
  DistributedFactorGraphs,
  DelimitedFiles,
  Statistics,
  Random,
  NLsolve,
  NLSolversBase,
  Optim,
  StatsBase,
  JLD2,
  FileIO,
  ProgressMeter,
  DocStringExtensions,
  FunctionalStateMachine,
  JSON2,
  Combinatorics

# experimental for replacing BayesTree on Graphs.jl
using MetaGraphs

using Logging

# bringing in BSD 3-clause ccolamd
include("ccolamd.jl")
using SuiteSparse.CHOLMOD: SuiteSparse_long # For CCOLAMD constraints.
using .Ccolamd


import Base: convert, ==
import Distributions: sample
import Random: rand, rand!
import KernelDensityEstimate: getBW
import ApproxManifoldProducts: kde!, manikde!
import ApproxManifoldProducts: mmd
import DistributedFactorGraphs: addVariable!, addFactor!, ls, lsf, isInitialized, compare, compareAllSpecial
import DistributedFactorGraphs: rebuildFactorMetadata!
import DistributedFactorGraphs: getDimension, getManifolds

# will be deprecated in IIF
import DistributedFactorGraphs: isSolvable


# must be moved to their own repos
const KDE = KernelDensityEstimate
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
  getSolvedCount, isSolved, setSolvedCount!,
  listSupersolves, listSolveKeys,
  deepcopySolvekeys!, deepcopySupersolve!,
  diagm,
  listDataEntries,
  FolderStore,
  addBlobStore!,
  getData
  # listDataBlobs  # ERROR: LightDFG{} doesn't override 'listDataBlobs'.

# Inference types
export FunctorInferenceType, PackedInferenceType
export AbstractPrior, AbstractRelativeFactor, AbstractRelativeFactorMinimize

export *,
  notifyCSMCondition,
  CSMHistory,
  # getTreeCliqsSolverHistories,

  updateCliqSolvableDims!,
  fetchCliqSolvableDims,
  AbstractBayesTree,
  BayesTreeNodeData,
  PackedBayesTreeNodeData,

  # state machine methods
  StateMachine,
  exitStateMachine,
  filterHistAllToArray,
  cliqHistFilterTransitions,
  printCliqSummary,
  printHistoryLine,
  printCliqHistorySummary,
  printCliqHistorySequential,
  printGraphSummary,
  printSummary,
  print,
  getGraphFromHistory,
  getCliqSubgraphFromHistory,
  sandboxStateMachineStep,
  sandboxCliqResolveStep,
  # draw and animate state machine
  getStateLabel,
  histStateMachineTransitions,
  histGraphStateMachineTransitions,
  drawStateTransitionStep,
  drawStateMachineHistory,
  animateStateMachineHistoryByTime,
  animateStateMachineHistoryByTimeCompound,
  animateCliqStateMachines,
  csmAnimate,
  makeCsmMovie,
  lockUpStatus!,
  unlockUpStatus!,
  lockDwnStatus!,
  unlockDwnStatus!,
  getSiblingsDelayOrder,
  areSiblingsRemaingNeedDownOnly,

  # general types for softtyping of variable nodes
  BeliefArray,
  InferenceVariable,
  ContinuousScalar,
  ContinuousMultivariate,
  SamplableBelief,
  Prior,
  PackedPrior,
  MsgPrior,
  PackedMsgPrior,
  PartialPrior,
  PackedPartialPrior,
  LinearConditional,
  PackedLinearConditional,
  MixturePrior,
  PackedMixturePrior,
  MixtureLinearConditional,
  PackedMixtureLinearConditional,

  ls2,
  # lsRear,
  # from DFG
  ls,
  lsf,
  listVariables,
  listFactors,
  exists,
  sortDFG,
  # getVariableIds,
  getVariableOrder,
  calcVariablePPE,
  getPPE,
  getPPEDict,
  getVariablePPE,
  isVariable,
  isFactor,
  # from dfg
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

  # using either dictionary or cloudgraphs
  FactorMetadata,
  FunctionNodeData,
  PackedFunctionNodeData, # moved to DFG
  normalfromstring,
  categoricalfromstring,
  extractdistribution,

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
  initVariable!,
  resetInitialValues!,
  resetInitValues!,
  asyncTreeInferUp!,
  initInferTreeUp!,
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
  getCliqSiblingsPriorityInitOrder,
  isCliqFullDim,
  getVariable,
  getCliqueData,
  setCliqueData!,
  getManifolds,
  getVal,
  getBW,
  setVal!,
  getNumPts,
  getBWVal,
  setBW!,
  setValKDE!,
  buildCliqSubgraph,
  cliqInitSolveUpByStateMachine!,

  # state machine functions
  finishCliqSolveCheck_StateMachine,
  determineCliqNeedDownMsg_StateMachine,
  blockUntilSiblingsStatus_StateMachine,
  doesCliqNeeddownmsg_StateMachine,
  slowCliqIfChildrenNotUpsolved_StateMachine,
  buildCliqSubgraph_StateMachine,
  isCliqUpSolved_StateMachine,
  testCliqCanRecycled_StateMachine,
  buildCliqSubgraphForDown_StateMachine,

  #
  isPartial,
  isInitialized,
  isTreeSolved,
  isUpInferenceComplete,
  isCliqInitialized,
  isCliqUpSolved,
  areCliqVariablesAllInitialized,
  areCliqChildrenNeedDownMsg,
  areCliqChildrenAllUpSolved,
  ensureSolvable!,
  ensureAllInitialized!,
  doCliqInitDown!,
  cycleInitByVarOrder!,
  prepCliqInitMsgsUp,
  prepCliqInitMsgsDown!,
  getOutNeighbors,
  BayesTree,
  TreeBelief,
  LikelihoodMessage,
  FullExploreTreeType,
  ExploreTreeType,
  initfg,
  buildSubgraph,
  buildCliqSubgraph!,
  transferUpdateSubGraph!,
  transferUpdateSubGraph!,
  getEliminationOrder,
  buildBayesNet!,
  emptyBayesTree,
  buildTree!,
  buildTreeFromOrdering!,
  resetBuildTreeFromOrder!,
  prepBatchTree!,
  wipeBuildNewTree!,
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
  solveTree!,
  solveGraph!,
  solveCliq!,
  fifoFreeze!,

  # temp const types TODO
  TempUpMsgPlotting,

  #functors need
  getSample,
  freshSamples!,
  freshSamples,

  #Visualization
  drawGraph,
  drawGraphCliq,
  drawCliqSubgraphUpMocking,
  drawTree,
  drawTreeAsyncLoop,
  # printgraphmax,

  # Bayes (Junction) Tree
  evalFactor2,
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
  isMultihypo,
  getMultihypoDistribution,
  getHypothesesVectors,
  isCliqMarginalizedFromVars,
  isCliqParentNeedDownMsg,
  setCliqAsMarginalized!,
  updateTreeCliquesAsMarginalizedFromVars!,

  # weiged sampling
  AliasingScalarSampler,
  rand!,
  rand,
  fastnorm,

  # new wrapper (experimental)
  CommonConvWrapper,

  getCliqVarInitOrderUp,
  getCliqInitVarOrderDown,
  blockCliqUntilChildrenHaveUpStatus,
  blockCliqSiblingsParentNeedDown,
  getCliqNumAssocFactorsPerVar,

  # introduced for approximate convolution operations
  setThreadModel!,
  SingleThreaded,
  MultiThreaded,

  # Solving utils
  findRelatedFromPotential,
  shuffleXAltD!,
  numericRoot,
  numericRootGenericRandomizedFnc!,

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

  uppA,
  convert,
  extractdistribution,

  # factor graph operating system utils (fgos)
  saveTree,
  loadTree,
  # landmarks,
  setCliqDrawColor,

  # Temp placeholder for evaluating string types to real types
  saveDFG,
  loadDFG!,  loadDFG,
  rebuildFactorMetadata!,

  getCliqVarSolveOrderUp,
  getFactorsAmongVariablesOnly,
  setfreeze!,

  #internal dev functions for recycling cliques on tree
  attemptTreeSimilarClique,
  getCliqSiblingsPartialNeeds,

  # some utils
  compare,
  compareAllSpecial,
  getIdx,
  getMeasurements,
  findFactorsBetweenFrom,
  addDownVariableFactors!,
  getDimension,
  setVariableRefence!,
  shuffleXAltD,
  reshapeVec2Mat,
  accumulateFactorChain,

  # For 1D example,
  # TODO rename to L2 distance
  Ranged,
  PackedRanged


const NothingUnion{T} = Union{Nothing, T}

# regular
include("FactorGraphTypes.jl")

# JT TODO move to somewhere more fitting? (DF, perhaps not remember its IIF.SolverParams)
const InMemDFGType = DFG.LightDFG{SolverParams}

import DistributedFactorGraphs: getFactorOperationalMemoryType

getFactorOperationalMemoryType(dfg::SolverParams) = CommonConvWrapper


include("AliasScalarSampling.jl")
include("DefaultNodeTypes.jl")
include("CliqueTypes.jl")
include("BeliefTypes.jl")
include("JunctionTreeTypes.jl")
include("FactorGraph.jl")
include("SerializingDistributions.jl")
include("DispatchPackedConversions.jl")
include("FGOSUtils.jl")
include("CompareUtils.jl")

# tree and init related functions
include("SubGraphFunctions.jl")
include("JunctionTree.jl")
include("TreeMessageAccessors.jl")
include("TreeMessageUtils.jl")
include("TreeMsgDwnConsolidation.jl")
include("TreeBasedInitialization.jl")

# solving graphs
include("SolverUtilities.jl")
include("DeconvUtils.jl")
include("ExplicitDiscreteMarginalizations.jl")
include("InferDimensionUtils.jl")
include("ApproxConv.jl")
include("SolveTree.jl")
include("TetherUtils.jl")
include("CliqStateMachine.jl")
include("CliqStateMachineUtils.jl")

#EXPERIMENTAL parametric
include("SolveTree_Parametric.jl")
include("CliqStateMachine_Parametric.jl")
include("ParametricUtils.jl")

# special variables and factors, see RoME.jl for more examples
include("GraphConstraintTypes.jl")
include("Variables/Sphere1D.jl")
include("Factors/Sphere1D.jl")
include("CanonicalGraphExamples.jl")

include("AdditionalUtils.jl")
include("SolverAPI.jl")
# Symbolic tree analysis files.
include("AnalysisTools.jl")

# deprecation legacy support
include("Deprecated.jl")



exportimg(pl) = error("Please do `using Gadfly` before IncrementalInference is used to allow image export.")
function __init__()
    @require InteractiveUtils = "b77e0a4c-d291-57a0-90e8-8db25a27a240" include("RequireInteractiveUtils.jl")

    @require Gadfly="c91e804a-d5a3-530f-b6f0-dfbca275c004" include("EmbeddedPlottingUtils.jl")
end

# Old code that might be used again
# function getType(typestring::AS) where {AS <: AbstractString}
#  # eval(Meta.parse(typestring))()
#  # getfield(Main, Symbol(typestring))
#  getfield(@__MODULE__, Symbol(typestring))
# end

export setSerializationNamespace!, getSerializationModule, getSerializationModules

end
