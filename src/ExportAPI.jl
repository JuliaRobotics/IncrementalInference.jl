# the IncrementalInference API


# DFG SpecialDefinitions
export AbstractDFG,
  getSolverParams,
  GraphsDFG,
  LocalDFG,
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
export AbstractPackedFactor, AbstractFactor
export AbstractPrior, AbstractRelative
export AbstractRelativeRoots, AbstractRelativeMinimize, AbstractManifoldMinimize

# not sure if this is necessary
export convert, *

export
  CSMHistory,
  # getTreeCliqsSolverHistories,

  AbstractBayesTree,
  BayesTreeNodeData,
  PackedBayesTreeNodeData,

  # state machine methods
  StateMachine,
  exitStateMachine,
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
  initVariable!,
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
  getBelief,
  CliqStateMachineContainer,

  solveCliqUp!,
  solveCliqDown!,
  fifoFreeze!,

  #functors need
  preambleCache,
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
  calcProposalBelief,
  approxConvBelief,
  approxConv,

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

  # Factor operational memory
  CommonConvWrapper,
  CalcFactor,
  
  getCliqVarInitOrderUp,
  getCliqNumAssocFactorsPerVar,

  # introduced for approximate convolution operations
  setThreadModel!,
  SingleThreaded,
  MultiThreaded,

  # user functions
  predictbelief,
  propagateBelief,
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
  reshapeVec2Mat


export incrSuffix

export calcPPE, calcVariablePPE
export setPPE!, setVariablePosteriorEstimates!
export getPPEDict
export getPPESuggested, getPPEMean, getPPEMax
export getPPESuggestedAll
export loadDFG
export findVariablesNear, defaultFixedLagOnTree!
export fetchDataJSON


export Position, Position1, Position2, Position3, Position4
export ContinuousScalar, ContinuousEuclid # TODO figure out if this will be deprecated, Caesar.jl #807
export Circular, Circle

# serializing distributions
export packDistribution, unpackDistribution
export PackedCategorical #, PackedDiscreteNonParametric
export PackedUniform, PackedNormal
export PackedZeroMeanDiagNormal, PackedZeroMeanFullNormal, PackedDiagNormal, PackedFullNormal
export PackedManifoldKernelDensity
export PackedAliasingScalarSampler, PackedHeatmapGridDensity, PackedLevelSetGridNormal
export PackedRayleigh

export Mixture, PackedMixture

export sampleTangent
export samplePoint

export  buildCliqSubgraph_StateMachine

