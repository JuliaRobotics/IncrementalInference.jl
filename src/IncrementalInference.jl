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
  hasVariable,
  getSolverParams,
  LightDFG,
  getSolvedCount, isSolved, setSolvedCount!,
  listSupersolves, listSolvekeys,
  deepcopySolvekeys!, deepcopySupersolve!,
  diagm

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
  # getCliqSolveHistory,
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
  exists,
  sortDFG,
  # getVariableIds,
  getVariableOrder,
  calcVariablePPE,
  getPPE,
  # getPPEs,
  getVariablePPE,
  # getVariablePPEs,
  # sortVarNested,
  # drawCopyFG,
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
  # VariableNodeData,
  # PackedVariableNodeData,
  FactorMetadata,
  # encodePackedType,
  FunctionNodeData,
  PackedFunctionNodeData, # moved to DFG
  # encodePackedType,
  # decodePackedType,
  normalfromstring,
  categoricalfromstring,
  extractdistribution,

  # FactorGraph,
  SolverParams,
  getSolvable,
  setSolvable!,
  # addNode!,
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
  # getVariablePotentialDims,
  getVariableSolvableDim,
  getFactorSolvableDim,
  getFactorInferFraction,
  getCliqSiblingsPriorityInitOrder,
  isCliqFullDim,
  getVariable,
  getCliqueData,
  setCliqueData!,
  getManifolds,
  # getVarNode,
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
  # doCliqInferAttempt_StateMachine,
  determineCliqNeedDownMsg_StateMachine,
  # blockUntilChildrenStatus_StateMachine,
  blockUntilSiblingsStatus_StateMachine,
  doesCliqNeeddownmsg_StateMachine,
  slowCliqIfChildrenNotUpsolved_StateMachine,
  # whileCliqNotSolved_StateMachine,
  buildCliqSubgraph_StateMachine,
  isCliqUpSolved_StateMachine,
  # determineAllChildrenNeedDownMsg_StateMachine,
  testCliqCanRecycled_StateMachine,
  buildCliqSubgraphForDown_StateMachine,

  #
  isPartial,
  isInitialized,
  isTreeSolved,
  isUpInferenceComplete,
  isCliqInitialized,
  # isCliqReadyInferenceUp,
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
  # updateFullVert!,
  getOutNeighbors,
  BayesTree,
  TreeBelief,
  # NBPMessage,
  LikelihoodMessage,
  FullExploreTreeType,
  ExploreTreeType,
  # FactorGraph,
  initfg,
  buildSubgraph,
  # copyGraph,
  buildCliqSubgraph!,
  # subgraphShortestPath,
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
  # isReadyCliqInferenceUp,
  childCliqs,
  getChildren,
  parentCliq,
  getParent,
  getCliqSiblings,
  getNumCliqs,
  getBelief, getKDE,
  # initializeNode!,
  CliqStateMachineContainer,
  solveTree!,
  solveGraph!,
  solveCliq!,
  fifoFreeze!,

  # temp const types TODO
  # TempBeliefMsg,
  TempUpMsgPlotting,

  #functors need
  getSample,
  freshSamples!,
  freshSamples,

  #Visualization
  # writeGraphPdf, # deprecated
  drawGraph,
  drawGraphCliq,
  drawCliqSubgraphUpMocking,
  drawTree,
  drawTreeAsyncLoop,
  # printgraphmax,

  # Bayes (Junction) Tree
  # evalPotential,
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

  # solve inference
  # inferOverTree!,
  # inferOverTreeR!,
  # inferOverTreeIterative!,

  #development interface
  # getTreeCliqSolveOrderUp,
  # getCliqOrderUpSolve,
  getCliqVarInitOrderUp,
  getCliqInitVarOrderDown,
  # getCliqueStatusUp,
  blockCliqUntilChildrenHaveUpStatus,
  blockCliqSiblingsParentNeedDown,
  getCliqNumAssocFactorsPerVar,
  # upMsgPassingRecursive,
  # downMsgPassingRecursive,

  # upMsgPassingIterative!,
  # downMsgPassingIterative!,

  # introduced for approximate convolution operations
  setThreadModel!,
  SingleThreaded,
  MultiThreaded,

  # Solving utils
  findRelatedFromPotential,
  shuffleXAltD!,
  numericRoot,
  # numericRootGenericRandomized,
  # numericRootGenericRandomizedFnc,
  numericRootGenericRandomizedFnc!,

  # user functions
  # proposalbeliefs,
  predictbelief,
  getCliqMat,
  getCliqAssocMat,
  getCliqMsgMat,
  getCliqFrontalVarIds,
  getFrontals,
  getCliqSeparatorVarIds,
  getCliqAllVarIds,
  getCliqVarIdsAll,
  # getCliqVars,
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

  # getSym,
  # doCliqInferenceUp!,
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

  # development
  # dev exports
# TODO deprecate
  # addGraphsVert!,
  # makeAddEdge!




# TODO should be deprecated
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

    @require Gadfly="c91e804a-d5a3-530f-b6f0-dfbca275c004" begin
      @info "Defining spyCliqMat(..) for visualizing association matrix of a clique in the Bayes (Junction) tree"

      exportimg(pl) = Gadfly.PNG(pl)

      export spyCliqMat

      """
          $SIGNATURES

      Draw the clique association matrix, with keyword arguments for more or less console print outs.

      Notes
      * Columns are variables, rows are factors.
      * Drawn from up message passing perspective.
      * Blue color implies no factor association.
      * Frontal, separator, and upmessages are all drawn at different intensity of red.
      * Downward messages not shown, as they would just be singletons of the full separator set.
      """
      function spyCliqMat(cliq::TreeClique; showmsg=true, suppressprint::Bool=false)
        mat = deepcopy(getCliqMat(cliq, showmsg=showmsg))
        # TODO -- add improved visualization here, iter vs skip
        mat = map(Float64, mat)*2.0.-1.0
        numlcl = size(getCliqAssocMat(cliq),1)
        mat[(numlcl+1):end,:] *= 0.9
        mat[(numlcl+1):end,:] .-= 0.1
        numfrtl1 = floor(Int,length(getCliqueData(cliq).frontalIDs) + 1)
        mat[:,numfrtl1:end] *= 0.9
        mat[:,numfrtl1:end] .-= 0.1
        if !suppressprint
          @show getCliqueData(cliq).itervarIDs
          @show getCliqueData(cliq).directvarIDs
          @show getCliqueData(cliq).msgskipIDs
          @show getCliqueData(cliq).directFrtlMsgIDs
          @show getCliqueData(cliq).directPriorMsgIDs
        end
        if size(mat,1) == 1
          mat = [mat; -ones(size(mat,2))']
        end
        sp = Gadfly.spy(mat)
        push!(sp.guides, Gadfly.Guide.title("$(getLabel(cliq)) || $(getCliqueData(cliq).frontalIDs) :$(getCliqueData(cliq).separatorIDs)"))
        push!(sp.guides, Gadfly.Guide.xlabel("fmcmcs $(getCliqueData(cliq).itervarIDs)"))
        push!(sp.guides, Gadfly.Guide.ylabel("lcl=$(numlcl) || msg=$(size(getCliqMsgMat(cliq),1))" ))
        return sp
      end
      function spyCliqMat(bt::AbstractBayesTree, lbl::Symbol; showmsg=true, suppressprint::Bool=false)
        spyCliqMat(getClique(bt,lbl), showmsg=showmsg, suppressprint=suppressprint)
      end
    end
end

# Old code that might be used again
# function getType(typestring::AS) where {AS <: AbstractString}
#  # eval(Meta.parse(typestring))()
#  # getfield(Main, Symbol(typestring))
#  getfield(@__MODULE__, Symbol(typestring))
# end

export setSerializationNamespace!, getSerializationModule, getSerializationModules

end
