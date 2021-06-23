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

export ℝ, AbsstractManifold, Euclidean, Circle

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


import Base: convert, ==
import Distributions: sample
import Random: rand, rand!
import KernelDensityEstimate: getBW
import KernelDensityEstimate: getPoints
import ApproxManifoldProducts: kde!, manikde!
import ApproxManifoldProducts: mmd
# import ApproxManifoldProducts: getManifolds # Deprecated
# import ApproxManifoldProducts: getManifold # might be used again later
import DistributedFactorGraphs: addVariable!, addFactor!, ls, lsf, isInitialized
import DistributedFactorGraphs: compare, compareAllSpecial
import DistributedFactorGraphs: rebuildFactorMetadata!
import DistributedFactorGraphs: getDimension, getManifold, getPointType, getPointIdentity
import DistributedFactorGraphs: getPPE, getPPEDict

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
  getCliqSiblingsPriorityInitOrder,
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

# regular
include("FactorGraphTypes.jl")

# JT TODO move to somewhere more fitting? (DF, perhaps not remember its IIF.SolverParams)
const InMemDFGType = DFG.LightDFG{SolverParams}

import DistributedFactorGraphs: getFactorOperationalMemoryType

getFactorOperationalMemoryType(dfg::SolverParams) = CommonConvWrapper


include("AliasScalarSampling.jl")
include("Flux/entities.jl")
include("BeliefTypes.jl")
include("CalcFactor.jl")

# Refactoring in progress
include("Factors/MsgLikelihoods.jl")

include("CliqueTypes.jl")

include("JunctionTreeTypes.jl")
include("FactorGraph.jl")
include("SerializingDistributions.jl")
include("DispatchPackedConversions.jl")

include("Variables/DefaultVariables.jl")

include("FGOSUtils.jl")
include("CompareUtils.jl")
include("NeedsResolution.jl")

# tree and init related functions
include("SubGraphFunctions.jl")
include("JunctionTree.jl")
include("TreeMessageAccessors.jl")
include("TreeMessageUtils.jl")
include("TreeBasedInitialization.jl")

# special variables and factors, see RoME.jl for more examples
include("GraphConstraintTypes.jl")
include("Factors/Mixture.jl")
include("Factors/DefaultPrior.jl")
include("Factors/LinearRelative.jl")
include("Factors/EuclidDistance.jl")
include("Factors/Circular.jl")
include("Variables/Circular.jl")
include("Factors/PartialPrior.jl")
include("DefaultNodeTypes.jl") # older file

# solving graphs
include("SolverUtilities.jl")
include("NumericalCalculations.jl")
include("DeconvUtils.jl")
include("ExplicitDiscreteMarginalizations.jl")
include("InferDimensionUtils.jl")
include("ApproxConv.jl")
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



exportimg(pl) = error("Please do `using Gadfly` before IncrementalInference is used to allow image export.")
function __init__()
  @require InteractiveUtils = "b77e0a4c-d291-57a0-90e8-8db25a27a240" include("RequireInteractiveUtils.jl")

  @require Gadfly="c91e804a-d5a3-530f-b6f0-dfbca275c004" include("EmbeddedPlottingUtils.jl")

  @require DifferentialEquations="0c46a032-eb83-5123-abaf-570d42b7fbaa" include("ODE/DERelative.jl")

  # combining neural networks natively into the non-Gaussian  factor graph object
  @require Flux="587475ba-b771-5e3f-ad9e-33799f191a9c" begin
    include("Flux/FluxModelsDistribution.jl")
    include("Flux/FluxModelsSerialization.jl") # uses BSON
  end
end


export setSerializationNamespace!, getSerializationModule, getSerializationModules

end
