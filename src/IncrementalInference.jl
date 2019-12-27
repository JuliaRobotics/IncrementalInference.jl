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
  DistributedFactorGraphs,
  DelimitedFiles,
  Statistics,
  Random,
  NLsolve,
  StatsBase,
  JLD2,
  FileIO,
  ProgressMeter,
  DocStringExtensions,
  FunctionalStateMachine,
  Optim, # might be deprecated in favor for only NLsolve dependency
  JSON2

using Logging

const KDE = KernelDensityEstimate
const AMP = ApproxManifoldProducts
const DFG = DistributedFactorGraphs
const FSM = FunctionalStateMachine

import Base: convert
# import HDF5: root
import Distributions: sample
import Random: rand, rand!
import KernelDensityEstimate: getBW
import ApproxManifoldProducts: kde!, manikde!
import DistributedFactorGraphs: addVariable!, addFactor!, ls, lsf, isInitialized, hasOrphans, compare, compareAllSpecial

# missing exports
import DistributedFactorGraphs: PackedFunctionNodeData, FunctionNodeData

# will be deprecated in IIF
import DistributedFactorGraphs: isSolvable


# TODO temporary for initial version of on-manifold products
KDE.setForceEvalDirect!(true)


export
  KDE,
  AMP,
  DFG,
  FSM,

  # DFG SpecialDefinitions
  AbstractDFG,
  hasVariable,
  getSolverParams,

  *,
  notifyCSMCondition,
  CSMHistory,
  getTreeCliqsSolverHistories,
  assignTreeHistory!,

  updateCliqSolvableDims!,
  fetchCliqSolvableDims,
  BayesTreeNodeData,
  PackedBayesTreeNodeData,

  # state machine methods
  StateMachine,
  exitStateMachine,
  getCliqSolveHistory,
  filterHistAllToArray,
  cliqHistFilterTransitions,
  printCliqSummary,
  printCliqHistorySummary,
  printGraphSummary,
  printSummary,
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
  lsRear,
  # from DFG
  ls,
  lsf,
  getVariableIds,
  getVariableOrder,
  calcVariablePPE,
  sortVarNested,
  hasOrphans,
  drawCopyFG,
  isVariable,
  isFactor,
  # from dfg
  getfnctype,
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

  # using either dictionary or cloudgraphs
  # VariableNodeData,
  # PackedVariableNodeData,
  FactorMetadata,
  getpackedtype,
  encodePackedType,
  FunctionNodeData,
  PackedFunctionNodeData, # moved to DFG
  encodePackedType,
  decodePackedType,
  normalfromstring,
  categoricalfromstring,
  extractdistribution,

  FactorGraph,
  SolverParams,
  addNode!,
  addVariable!,
  deleteVariable!,
  addFactor!,
  deleteFactor!,
  addMsgFactors!,
  deleteMsgFactors!,
  factorCanInitFromOtherVars,
  doautoinit!,
  manualinit!,
  initVariable!,
  asyncTreeInferUp!,
  initInferTreeUp!,
  solveCliqWithStateMachine!,
  resetData!,
  resetTreeCliquesForUpSolve!,
  resetFactorGraphNewTree!,
  setVariableInitialized!,
  setVariableInferDim!,
  resetVariable!,
  resetCliqSolve!,
  getFactor,
  getFactorDim,
  getVariableDim,
  getVariableInferredDim,
  getVariableInferredDimFraction,
  getVariablePotentialDims,
  getVariableSolvableDim,
  getFactorSolvableDim,
  getFactorInferFraction,
  getCliqSiblingsPriorityInitOrder,
  isCliqFullDim,
  getVariable,
  # getVert, # deprecated use DFG.getVariable getFactor instead
  getData,
  setData!,
  getManifolds,
  getVarNode,
  getVal,
  getBW,
  setVal!,
  getNumPts,
  getBWVal,
  setBW!,
  setValKDE!,
  buildCliqSubgraph,
  # buildCliqSubgraphUp,
  # buildCliqSubgraphDown,
  setCliqUpInitMsgs!,
  cliqInitSolveUpByStateMachine!,

  # state machine functions
  finishCliqSolveCheck_StateMachine,
  doCliqInferAttempt_StateMachine,
  determineCliqNeedDownMsg_StateMachine,
  blockUntilChildrenStatus_StateMachine,
  blockUntilSiblingsStatus_StateMachine,
  doesCliqNeeddownmsg_StateMachine,
  slowCliqIfChildrenNotUpsolved_StateMachine,
  whileCliqNotSolved_StateMachine,
  buildCliqSubgraph_StateMachine,
  isCliqUpSolved_StateMachine,
  determineAllChildrenNeedDownMsg_StateMachine,
  testCliqCanRecycled_StateMachine,
  buildCliqSubgraphForDown_StateMachine,

  #
  isPartial,
  isInitialized,
  isTreeSolved,
  isUpInferenceComplete,
  isCliqInitialized,
  isCliqReadyInferenceUp,
  isCliqUpSolved,
  areCliqVariablesAllInitialized,
  areCliqChildrenNeedDownMsg,
  areCliqChildrenAllUpSolved,
  ensureAllInitialized!,
  doCliqAutoInitUpPart1!,
  doCliqAutoInitUpPart2!,
  doCliqInitDown!,
  cycleInitByVarOrder!,
  doCliqUpSolve!,
  getCliqInitUpMsgs,
  getCliqStatus,
  setCliqStatus!,
  getSolveCondition,
  prepCliqInitMsgsUp,
  prepCliqInitMsgsDown!,
  updateFullVert!,
  getOutNeighbors,
  BayesTree,
  EasyMessage,
  NBPMessage,
  FullExploreTreeType,
  ExploreTreeType,
  FactorGraph,
  initfg,
  buildSubgraphFromLabels,
  subgraphShortestPath,
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
  hasCliq,
  getCliq,
  whichCliq,
  getCliqDepth,
  getTreeAllFrontalSyms,
  getTreeCliqUpMsgsAll,
  stackCliqUpMsgsByVariable,
  getCliqChildMsgsUp,
  getCliqParentMsgDown,
  getCliqDownMsgsAfterDownSolve,
  isReadyCliqInferenceUp,
  childCliqs,
  getChildren,
  parentCliq,
  getParent,
  getCliqSiblings,
  getKDE,
  getVertKDE,
  initializeNode!,
  CliqStateMachineContainer,
  batchSolve!,
  solveTree!,
  solveCliq!,
  fifoFreeze!,
  getCurrentWorkspaceFactors,
  getCurrentWorkspaceVariables,

  # temp const types TODO
  TempBeliefMsg,
  TempUpMsgPlotting,

  #functors need
  getSample,
  freshSamples!,

  #Visualization
  writeGraphPdf, # deprecated, but first move code to drawGraph before deleting
  drawGraph,
  drawGraphCliq,
  drawCliqSubgraphUpMocking,
  drawTree,
  printgraphmax,
  # allnums,
  # isnestednum,
  # sortnestedperm,

  # Bayes (Junction) Tree
  evalPotential,
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
  inferOverTree!,
  inferOverTreeR!,
  inferOverTreeIterative!,

  #development interface
  getTreeCliqSolveOrderUp,
  getCliqOrderUpSolve,
  getCliqVarInitOrderUp,
  getCliqInitVarOrderDown,
  getCliqStatusUp,
  blockCliqUntilChildrenHaveUpStatus,
  blockCliqSiblingsParentNeedDown,
  getCliqNumAssocFactorsPerVar,
  upMsgPassingRecursive,
  downMsgPassingRecursive,

  upMsgPassingIterative!,
  downMsgPassingIterative!,

  # Inference types
  InferenceType,
  PackedInferenceType,
  Singleton,
  Pairwise,
  # introduced for approximate convolution operations
  setThreadModel!,
  SingleThreaded,
  MultiThreaded,

  # functor abstracts
  FunctorInferenceType,
  FunctorPairwise,
  FunctorPairwiseMinimize,
  FunctorSingleton,
  # FunctorPartialSingleton,
  FunctorPairwiseNH,
  FunctorSingletonNH,

  # Solving utils
  findRelatedFromPotential,
  shuffleXAltD!,
  numericRoot,
  numericRootGenericRandomized,
  numericRootGenericRandomizedFnc,
  numericRootGenericRandomizedFnc!,
  solveFactorMeasurements,

  # user functions
  proposalbeliefs,
  predictbelief,
  getCliq,
  getCliqMat,
  getCliqAssocMat,
  getCliqMsgMat,
  getCliqFrontalVarIds,
  getFrontals,                     # duplicate
  getCliqSeparatorVarIds,
  getCliqAllVarIds,
  getCliqVarIdsAll,
  getCliqVars,
  getCliqAllVarSyms,
  getCliqVarIdsPriors,
  getCliqVarSingletons,
  getCliqAllFactIds,
  getCliqFactorIdsAll,
  getCliqFactors,
  areCliqVariablesAllMarginalized,
  setTreeCliquesMarginalized!,

  # generic marginal used during elimitation game
  GenericMarginal,
  PackedGenericMarginal,

  uppA,
  convert,
  extractdistribution,

  # factor graph operating system utils (fgos)
  saveTree,
  loadTree,
  # convert2packedfunctionnode,
  # encodefg,
  # decodefg,
  # savejld,
  # loadjld,
  landmarks,
  setCliqDrawColor,

  # csm utils
  fetchCliqTaskHistoryAll!,
  fetchAssignTaskHistoryAll!,

  # Temp placeholder for evaluating string types to real types
  _evalType,
  saveDFG,
  loadDFG,
  rebuildFactorMetadata!,

  setUpMsg!,
  upMsg,
  getDwnMsgs,
  getUpMsgs,
  setDwnMsg!,
  dwnMsg,
  getDwnMsgs,
  getCliqMsgsUp,
  getCliqMsgsDown,
  getCliqVarSolveOrderUp,

  getSym,
  doCliqInferenceUp!,
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

  # For 1D example,

  # TODO rename to L2 distance
  Ranged,
  PackedRanged,

  # development
  # dev exports
  addGraphsVert!,
  makeAddEdge!,
  shuffleXAltD,
  reshapeVec2Mat, # TODO deprecate

  # OBSOLETE TODO REMOVE #TODO TODO
  subGraphFromVerts,
  getMaxVertId,
  dlapi,
  localapi,
  showcurrentdlapi,
  setdatalayerAPI!,
  DataLayerAPI



# TODO should be deprecated
const NothingUnion{T} = Union{Nothing, T}

# non-free, but not currently use!
include("ccolamd.jl")

# regular
include("FactorGraphTypes.jl")

# const InMemDFGType = DFG.LightDFG{SolverParams} #swap out default in v0.8.0/v0.9.0?
const InMemDFGType = DFG.GraphsDFG{SolverParams} # JT TODO move to somewhere more fitting?

include("BeliefTypes.jl")
include("AliasScalarSampling.jl")
include("DefaultNodeTypes.jl")
include("FactorGraph01.jl")
include("SerializingDistributions.jl")
include("DispatchPackedConversions.jl")
include("FGOSUtils.jl")
include("CompareUtils.jl")

include("JunctionTreeTypes.jl")
include("SubGraphFunctions.jl")
include("JunctionTree.jl")
include("TreeBasedInitialization.jl")
include("GraphConstraintTypes.jl")
include("SolverUtilities.jl")
include("ExplicitDiscreteMarginalizations.jl")
include("InferDimensionUtils.jl")
include("ApproxConv.jl")
include("SolveTree01.jl")
include("TetherUtils.jl")
include("CliqStateMachine.jl")
include("CliqStateMachineUtils.jl")

# special variables and factors, see RoME.jl for more examples
include("Variables/Sphere1D.jl")
include("Factors/Sphere1D.jl")

include("AdditionalUtils.jl")
include("SolverAPI.jl")

include("CanonicalGraphExamples.jl")
include("Deprecated.jl")

exportimg(pl) = error("Please do `using Gadfly` before IncrementalInference is used to allow image export.")
function __init__()
    @require Gadfly="c91e804a-d5a3-530f-b6f0-dfbca275c004" begin
      @info "Defining spyCliqMat(..) for visualizing association matrix of a clique in the Bayes (Juntion) tree"

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
      function spyCliqMat(cliq::Graphs.ExVertex; showmsg=true, suppressprint::Bool=false)
        mat = deepcopy(getCliqMat(cliq, showmsg=showmsg))
        # TODO -- add improved visualization here, iter vs skip
        mat = map(Float64, mat)*2.0.-1.0
        numlcl = size(getCliqAssocMat(cliq),1)
        mat[(numlcl+1):end,:] *= 0.9
        mat[(numlcl+1):end,:] .-= 0.1
        numfrtl1 = floor(Int,length(getData(cliq).frontalIDs) + 1)
        mat[:,numfrtl1:end] *= 0.9
        mat[:,numfrtl1:end] .-= 0.1
        if !suppressprint
          @show getData(cliq).itervarIDs
          @show getData(cliq).directvarIDs
          @show getData(cliq).msgskipIDs
          @show getData(cliq).directFrtlMsgIDs
          @show getData(cliq).directPriorMsgIDs
        end
        sp = Gadfly.spy(mat)
        push!(sp.guides, Gadfly.Guide.title("$(cliq.attributes["label"]) || $(cliq.attributes["data"].frontalIDs) :$(cliq.attributes["data"].separatorIDs)"))
        push!(sp.guides, Gadfly.Guide.xlabel("fmcmcs $(cliq.attributes["data"].itervarIDs)"))
        push!(sp.guides, Gadfly.Guide.ylabel("lcl=$(numlcl) || msg=$(size(getCliqMsgMat(cliq),1))" ))
        return sp
      end
      function spyCliqMat(bt::BayesTree, lbl::Symbol; showmsg=true, suppressprint::Bool=false)
        spyCliqMat(whichCliq(bt,lbl), showmsg=showmsg, suppressprint=suppressprint)
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
